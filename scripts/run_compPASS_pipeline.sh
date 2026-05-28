#!/bin/bash
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

# Main script for evaluating profile of chromatin reads from Pol2 ChIP/HiChIP data
# Yaw Asante | yxa181@case.edu | August 2025 | Gryder Lab, CWRU


# STEP 0: Reads in configuration file, includes path to BAM files and conda environments
loc=$(pwd)
source $loc/scripts/compass_config_file.sh


# STEP 1: Handles user inputs
function usage {
    echo -e "usage: "
    echo -e "  run_compPASS_pipeline.sh \\"
    echo -e "    -A sample_a_input  \\" # the control BAM file (prefix)
    echo -e "    -B sample_b_input  \\" # the perturbed BAM file (prefix)
    echo -e "   [ -g GENE_REF     ]  \\" #
    echo -e "   [ -n sample_names ]  \\" # sampleA_name,sampleB_name (comma separated)
    echo -e "   [ -m not_spikein  ]  \\" # 
    echo -e "   [ -r resegment    ]  \\" # TRUE or FALSE (uses updated compass_config_file.sh to make new segments)
    echo -e "   [ -o output_dir   ]  \\" # default: creates a folder named "results"
    echo -e "   [ -v version      ]  \\" #
    echo -e "   [ -h]"
    echo -e "Use option -h|--help for more information"
    exit;
}

# tip 1: to use direct paths to input, set DATASOURCE="" in compass_config_file.sh
# tip 2: to use direct paths to output, set OUTPUTDIR="" in compass_config_file.sh
# tip 3: CAPITAL options are required, lowercase options are optional

function help {
    echo "Starting from RNA Pol2 reads (bam) for control and case samples, calculate and visualize Pol2 defect modes"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -A|--sample_a       : Prefix of sample A bam file in your designated input folder. Used to name bigwig"
    echo "   -B|--sample_b       : Prefix of sample B bam file in your designated input folder. Used to name bigwig"
    echo " [ -g|--gene_ref      ] : BED file of gene chr, start, end and strand information. Defaults to hg38 RefSeq intervals for genes of length > 2kb."
    echo " [ -n|--sample_names  ] : Names of samples A and B to be used in plots (format: name_a,name_b) Defaults to control,case."
    echo " [ -m|--not_use_spike ] : Do not use spike in. If T, uses RPM normalization. Defaults to F."
    echo " [ -r|--resegment     ] : Set to T to resegment gene reference into regions of interest. Defaults to F."
    echo " [ -o|--output_dir    ] : Specify the name of the folder where results go in your designated output folder. Defaults to results."
    echo " [ -v|--version       ] : Displays the version."
    echo " [ -h|--help          ] : Help message."
    exit;
}

function version {
    echo "version 0.0.5, Feb 4th 2026"
    exit;
}


# Transforms long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--sample_a")      set -- "$@" "-A" ;;
      "--sample_b")      set -- "$@" "-B" ;;
      "--gene_ref")      set -- "$@" "-g" ;;
      "--sample_names")  set -- "$@" "-n" ;;
      "--not_use_spike") set -- "$@" "-m" ;;
      "--resegment")     set -- "$@" "-r" ;;
      "--output_dir")    set -- "$@" "-o" ;;
      "--version")       set -- "$@" "-v" ;;
      "--help")          set -- "$@" "-h" ;;
       *)                set -- "$@" "$arg"
  esac
done

if [ $# -lt 1 ];then
    usage
    exit 1;
fi

g=$GENEREF
n="sample_a,sample_b"
r="F"
m="F"
o="results"


while getopts ":A:B:g:n:m:r:o:v:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  B) B=$OPTARG;;
  g) g=$OPTARG;;
  n) n=$OPTARG;;
  m) m=$OPTARG;;
  r) r=$OPTARG;;
  o) o=$OPTARG;;
  v) version  ;;
  h) help     ;;
  \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
  :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
    esac
done


m=$(echo "$m" | tr '[:lower:]' '[:upper:]')
if [[ $m == "TRUE" ]];then
    m="T";
fi

if [[ $m == "FALSE" ]];then
    m="F";
fi

if [[ $m != "T" && $m != "F" ]];then
    echo "Normalization variable not recognized, using RPM";
    m="F";
fi

n1=$(echo "${n}" | cut -d "," -f 1)
n2=$(echo "${n}" | cut -d "," -f 2)
remake_segments="$r"
outdir="$o"
src_bam1=$(realpath "$DATASOURCE/$A/$A.bam")
src_bam2=$(realpath "$DATASOURCE/$B/$B.bam")

if [[ ! -f "${src_bam1}" ]]; then
    echo "sample_a not found. Exiting."
    usage
    exit 1
fi

if [[ ! -f "${src_bam2}" ]]; then
    echo "sample_b not found. Exiting."
    usage
    exit 1
fi

# Activates conda environment
echo ""
echo "Step 1: Activate conda environment"
source activate $CONDAVAL
# used on cluster
#module load SAMtools
#module load BEDTools/2.30.0-GCC-11.2.0
echo " - Env active."

# Makes a new directory in output
if [[ -d $OUTPUTDIR"/"$outdir ]];then
    echo " - Warning, directory $o already exists, content may be overwritten ..."
else   
    mkdir "$OUTPUTDIR/$outdir"
fi
shortout=$OUTPUTDIR"/"$outdir

# STEP 2: If not present or requested, (re)generates segmented gene BED3+3 from gene reference (GENEREF, BED3+2)
echo ""
echo "Step 2: Check for segmented BED reference"

#SEGMENTED=$(echo "$GENEREF" | sed s/".tsv"/"_segments.bed"/g)
SEGMENTED_BASENAME=$(basename "$GENEREF" .tsv)_segments.bed

#if [[ ! -f "$REFSOURCE/${SEGMENTED_BASENAME}" || ${remake_segments} == "T" ]];then
#    bash ${SCRIPTSOURCE}/make_segmented_bed.sh $GENEREF $REFSOURCE/${SEGMENTED_BASENAME} $SCRIPTSOURCE
#    #echo "bash ${SCRIPTSOURCE}/make_segmented_bed.sh $GENEREF $shortout/$SEGMENTED $SCRIPTSOURCE"
#else
#    echo " - Segmented BED already made."
#fi
#
SEGMENTED="$REFSOURCE/${SEGMENTED_BASENAME}"
# NOTE: file is created in sorted order; downstream conversions rely on all 4 segments being consecutive

# STEP 3: Uses input to generate Pol2 array for samples A and B
echo ""
echo "Step 3: Make annotated coverage map from segmented gene reference."

coverage_file=${n1}"_vs_"${n2}"_gene_coverage.bed"
if [[ ! -s "$shortout/${coverage_file}" ]]; then
        echo " - Reading from input BAMs, please be patient ..."
        bedtools multicov -bams "${src_bam1}" "${src_bam2}" -bed $SEGMENTED > $shortout"/"${coverage_file}
        echo " - Coverage BED made"
else
     # check if full annotated file made
    linesGoal=$(wc -l $SEGMENTED | cut -d " " -f 1)
    linesActual=$(wc -l $shortout"/"${coverage_file} | cut -d " " -f 1)

    if [[ $linesActual != $linesGoal ]]; then
        echo " - Incomplete coverage BED found, remaking ... "
        bedtools multicov -bams ${src_bam1} ${src_bam2} -bed $SEGMENTED > $shortout"/"${coverage_file}
        echo " - Coverage BED made"
        #echo "bedtools multicov -bams ${src_bam1} ${src_bam2} -bed $SEGMENTED > ${coverage_file}"
    else
        echo " - Coverage for ${n1}_${n2} comparison already made. Continuing ..."
    fi
fi


# STEP 4a: Get normalization read counts
echo ""
echo "Step 4a: Get read counts for normalization"

counts1=0
counts2=0

if [[ $m == "F" ]]; then
    echo "Using RRPM normalization."
    # using spike-in, grab from SpikeIn folder
    spikein1=${DATASOURCE}"/"$A"/SpikeIn/spike_map_summary"
    spikein2=${DATASOURCE}"/"$B"/SpikeIn/spike_map_summary"

    if [ ! -f "${spikein1}" ]; then
        echo "Spike-in values for ${n1} not found. Exiting ..."
        exit
    else
        counts1=$(python3 $SCRIPTSOURCE"/get_spikein_count.py" $spikein1)
    fi

    if [ ! -f "${spikein2}" ]; then
        echo "Spike-in values for ${n2} not found. Exiting ..."
        exit
    else
        counts2=$(python3 $SCRIPTSOURCE"/get_spikein_count.py" $spikein2)
    fi

    read_counts_csp="${counts1},${counts2}"

else
    echo "Using RPM normalization."
    # not using spike-in, grab from flagstat.txt
    flagstat1="${DATASOURCE}/$A/$A.flagstat.txt"
    flagstat2="${DATASOURCE}/$B/$B.flagstat.txt"

    if [ ! -f "${flagstat1}" ]; then
        echo "Flagstat/readcounts for ${n1} not found. Exiting ..."
        exit
    else
        counts1=$(python3 $SCRIPTSOURCE"/get_primary_count.py" $flagstat1)
    fi

    if [ ! -f "${flagstat2}" ]; then
        echo "Flagstat/readcounts for ${n2} not found. Exiting ..."
        exit
    else
        counts2=$(python3 $SCRIPTSOURCE"/get_primary_count.py" $flagstat2)
    fi

    read_counts_csp="${counts1},${counts2}"
fi

echo " - Counts for sample a - "${n1}": "${counts1}
echo " - Counts for sample b - "${n2}": "${counts2}

# STEP 4b: Uses pol2 read lists (BED3+4) from samples to generate Log2FC table (Raw data output )
echo ""
echo "Step 4b: Make Pol2 ratio table"

comparison_file=${n1}"_vs_"${n2}"_comp_pol2.tsv"
if [[ ! -s "$shortout/${comparison_file}" ]]; then
    Rscript $SCRIPTSOURCE"/generate_comparison_table_from_multicov.R" $shortout ${coverage_file} \
     ${read_counts_csp} ${pro_start} ${tssr_start} ${gene_start} ${tesr_end} \
     ${SCRIPTSOURCE} ${EXCLUDE_CHR_Y} $shortout"/"${comparison_file}
     #echo "Pol2 table made"
    #echo "Rscript $SCRIPTSOURCE/generate_comparison_table_from_multicov.R $shortout ${coverage_file} ${read_counts_csp} ${pro_start} ${tssr_start} ${gene_start} ${tesr_end} > ${comparison_file}"
else
    echo " - Pol2 table for ${n1}_vs_${n2} already made. Continuing ..."
fi


# STEP 5: Uses Log2FC table to classify genes and draw plots (Raw data output + Visual outputs)
echo ""
echo "Step 5: Draw comparison plots from ratio table"
mkdir -p $shortout/${n1}_vs_${n2}_genelists
output_table=${n1}"_vs_"${n2}"_pol2_ratio_tb.tsv"

Rscript $SCRIPTSOURCE/classify_genes_from_pol2_states.R $shortout ${comparison_file} \
 ${SCRIPTSOURCE} ${GCT} $n1 $n2 $shortout"/"${output_table}
#echo "Gene classifications made"


# STEP 6: Generates profile plots across gene body regions for classified genes (Visual outputs)
echo ""
echo "Step 6: Generate gene profile plots"

mkdir -p $shortout"/metagene_plots" 
# copy genelists over to heatmap folder
cp $shortout/${n1}_vs_${n2}_genelists/*.txt $shortout"/metagene_plots" 

# create config file
genes_to_plot=$(ls $shortout"/"${n1}_vs_${n2}"_genelists")
namelist=""
bedlist=""
for list_to_plot in ${genes_to_plot};do
    curname=$(echo "${list_to_plot}" | cut -d "." -f 1)
    setname=$(echo "${list_to_plot}" | rev | cut -d "_" -f 1,2 | rev)

    namelist=$namelist","${setname}
    bedlist=$bedlist","${curname}".txt"
done 

if [[ $namelist != "" ]];then
    namelist="${namelist:1}"
    bedlist="${bedlist:1}"

    # copy meta-gene config template, replace sample names, add BAM absolute paths and get read lengths
    exampleread=$(samtools view ${src_bam1} | head -1 | cut -f 10)
    readlen=${#exampleread}
    echo " - Expected read length from BAMs: "$readlen

    maxreads=$(( counts1 > counts2 ? counts1 : counts2 ))
    readcap=$(( (maxreads + 20000000) / 20000000))

    if [[ $m == "F" ]]; then
        readcap=$((readcap * 2));
        echo " - Max RRPM value for profile plots: "$readcap
    else
		readcap=$((readcap / 2));
        echo " - Max RPM value for profile plots: "$readcap
    fi
    #cp $SCRIPTSOURCE"/heatmap_config_template.txt" $shortout"/heatmap.config"

    cat $SCRIPTSOURCE"/heatmap_config_template.txt" |\
            sed  "s#%SAMPLE_A_BAM%#$src_bam1#" |\
            sed  "s#%SAMPLE_B_BAM%#$src_bam2#" |\
            sed  "s#%BAM_READ_LENGTH%#$readlen#" |\
            sed  "s#%MAX_VALUE_PLOT%#$readcap#" |\
            sed  "s#%SAMPLE_A_NAME%#$n1#" |\
            sed  "s#%SAMPLE_B_NAME%#$n2#" > $shortout"/heatmap.config"


    cp $SCRIPTSOURCE"/genes_refseq.hg38fix.bed" $shortout/metagene_plots
    cd $shortout/metagene_plots

    echo " - Generating profile plots"
    #echo "python3 $SCRIPTSOURCE/plotHeatmap.py -b $bedlist -c $shortout/heatmap.config \
    #                        -l $namelist -o $shortout/metagene_plots -n 10 -s 30 \
    #                        -m $m -g 1 -t $THREADCOUNT -r $r"

    if [[ $m == "T" || m == "True" || m == "TRUE" ]]; then
    python3 $SCRIPTSOURCE/plotHeatmap.py -b $bedlist -c $shortout/heatmap.config \
                            -l $namelist -o $shortout/metagene_plots -n 10 -s 30 \
                            -e $leftgenebp -y $genebodybp -k $rightgenebp \
                            -m -g 1 -t $THREADCOUNT
    else
    python3 $SCRIPTSOURCE/plotHeatmap.py -b $bedlist -c $shortout/heatmap.config \
                            -l $namelist -o $shortout/metagene_plots -n 10 -s 30 \
                            -e $leftgenebp -y $genebodybp -k $rightgenebp \
                            -g 1 -t $THREADCOUNT
    fi
    
else
    echo "No major gene changes found, no profile plots made"
fi

# STEP 7: Cleans up and removes temporary files
#echo "Removing temporary files ..."
mv $LOCATION/slurm-${SLURM_JOB_ID}.out $shortout/run_log.txt
chmod 774 $shortout/*
#chgrp beg33 $shortout/* 

#bash cleaning_up.sh

echo ""
echo "Process complete. Thank you for using compPASS!"
