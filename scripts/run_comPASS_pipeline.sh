#!/bin/bash
# Main script for evaluating profile of chromatin reads from Pol2 ChIP/HiChIP data
# Yaw Asante | yxa181@case.edu | January 7th, 2025 | Gryder Lab, CWRU


# STEP 0: Reads in configuration file, includes path to BAM files and conda environments
source config_file.sh

# STEP 1: Handles user inputs
function usage {
    echo -e "usage: "
    echo -e "  run_full_comPAS_pipeline.sh \\"
    echo -e "    -A sample_a_input  \\"
    echo -e "    -B sample_b_input  \\"
    echo -e "   [-g GENE_REF ]      \\"
    echo -e "   [-n sample_names ]  \\"
    echo -e "   [-f scale_factors ] \\"
    echo -e "   [-r resegment ]     \\"
    echo -e "   [-o output_dir ]    \\"
    echo -e "   [-v version    ]    \\"
    echo -e "   [-h]"
    echo -e "Use option -h|--help for more information"
}

# tip 1: to use direct paths to input, set DATASOURCE="" in system_config.sh
# tip 2: to use direct paths to output, set OUTPUTDIR="" in system_config.sh
# tip 3: CAPITAL options are required, lowercase options are optional

function help {
    echo 
    echo "Starting from RNA Pol2 reads (bam) for two samples, calculate and visualize Pol2 defect modes"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -A|--sample_a       : Name of sample A bam file in your designated input folder."
    echo "   -B|--sample_b       : Name of sample B bam file in your designated input folder."
    echo " [ -g|--gene_ref     ] : BED file of gene chr, start, end and strand information. Defaults to hg38 RefSeq intervals for genes of length > 330."
    echo " [ -n|--sample_names ] : Names of samples A and B to be used in plots (format: name_a,name_b) Defaults to sample_a,sample_b."
    echo " [ -f|--scalefactors ] : Scale factors (format: 123,123) for each sample, mimicking spike-in. Defaults to 1,1."
    echo " [ -r|--resegment    ] : Set to T to resegment gene reference into regions of interest. Defaults to F."
    echo " [ -o|--output_dir   ] : Specify the name of the folder where results go in your designated output folder. Defaults to results."
    echo " [ -v|--version      ] : Displays the version."
    echo " [ -h|--help         ] : Help message."
    exit;
}

function version {
    echo "version 1.0, January 8th 2025"
    exit;
}

# Transforms long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--sample_a")     set -- "$@" "-A" ;;
      "--sample_b")     set -- "$@" "-B" ;;
      "--gene_ref")     set -- "$@" "-g" ;;
      "--sample_names") set -- "$@" "-n" ;;
      "--scalefactors") set -- "$@" "-f" ;;
      "--resegment")    set -- "$@" "-r" ;;
      "--output_dir")   set -- "$@" "-o" ;;
      "--version")      set -- "$@" "-v" ;;
      "--help")         set -- "$@" "-h" ;;
       *)               set -- "$@" "$arg"
  esac
done

g=$GENEREF
n="sample_a,sample_b"
f="1,1"
r="F"
o="results"


while getopts ":A:B:g:n:f:v:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  B) B=$OPTARG;;
  g) g=$OPTARG;;
  n) n=$OPTARG;;
  f) f=$OPTARG;;
  r) r=$OPTARG;;
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

if [ -z "$A" ]; then
    echo "sample_a not found. Exiting."
    usage
    exit 1
fi

if [ -z "$B" ]; then
    echo "sample_b not found. Exiting."
    usage
    exit 1
fi

# Activates conda environment
source activate $CONDAPATH"/comPAS"

# Makes a new directory in output
if [ -d $OUTPUTDIR"/"$o ];then
    echo "Warning, directory $o already exists, content may be overwritten ..."
else   
    mkdir $OUTPUTDIR"/"$o
fi
shortout=$OUTPUTDIR"/"$o

# STEP 2: If not present or requested, (re)generates segmented gene BED3+3 from gene reference (GENEREF, BED3+2)
SEGMENTED=$(echo "$GENEREF" | sed s/".tsv"/"_segments.bed"/g)

if [ -z "$SEGMENTED" || "$r" == "T" ];then
    bash $SCRIPTSOURCE"/"make_segmented_bed.sh $GENREF $SEGMENTED
fi
# note, file is created in sorted order; downstream conversions rely on all 4 segments being consecutive

# STEP 3: Uses input to generate Pol2 array for samples A and B
echo "Reading from input bam files, please be patient ..."
#if [ "$T" == "bigwig"]; then
    #python3 $SCRIPTSOURCE"/get_table_from_bigwig.py" -a $DATASOURCE"/"$A -b $DATASOURCE"/"$B -s $SEGMENTED -o $shortout
#    multiBigwigSummary BED-file -b $DATASOURCE"/"$A $DATASOURCE"/"$B --BED $SEGMENTED --outRawCounts $shortout"/gene_coverage.bed" 
#else
    bedtools multicov -bams $DATASOURCE"/"$A $DATASOURCE"/"$B -bed $SEGMENTED > $shortout"/gene_coverage.bed"
#fi

# STEP 4: Uses pol2 read lists (BED3+4) from samples to generate Log2FC table (Raw data output )
Rscript $SCRIPTSOURCE"/generate_comparison_table_from_multicov.R" $shortout "gene_coverage.bed" $f ${pro_start} ${tssr_start} ${gene_start} ${tesr_end}

# STEP 5: Uses Log2FC table to classify genes and draw plots(Raw data output + Visual outputs)
Rscript $SCRIPTSOURCE"/classify_genes_from_pol2_states.R" $shortout "comparison_tb.tsv" 

# STEP 6: Uses provided gene_categories table to further classify genes (Raw data output + Visual outputs)
#echo "Step 6: Classifying genes of interest according to ${GCT} ... "


# STEP 7: Generates profile plots across gene body regions for classified genes (Visual outputs)
echo "Step 7: Generating gene profile plots"

# reads sense, left to right
leftgenebp="2000"  # extends past promoter region
genebodybp="5000"  # gene body regions scale, used for relative length
rightgenebp="4000" # matches tesr region

python3 $SCRIPTSOURCE"/plot_geneprofiles.py" $shortout "total_comp_pol2_class.tsv" $leftgenebp $rightgenebp $PROCESSORS


# STEP 8: Cleans up and removes temporary files
#echo "Removing temporary files ..."
#bash cleaning_up.sh

echo "Process complete. Thank you for using comPAS!"