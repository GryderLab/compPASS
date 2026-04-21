#!/bin/bash
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Use to replot heatmaps at a selected RPM value, reuses an existing config file


source="/mnt/vstor/SOM_GENE_BEG33/ChIP_seq/hg38/projects/Pol2_ratios/compPASS_outs/"
n1=$1
n2=$2
shortout=${source}/$3
MAX=$4
not_use_spike=$5
output=$6


m=${not_use_spike}
m=$(echo "$m" | tr '[:lower:]' '[:upper:]')
if [[ m == "TRUE" ]];then
    m = "T";
fi

if [[ m == "FALSE" ]];then
    m = "F";
fi

loc=$(pwd)
source $loc/scripts/compass_config_file.sh
module load BEDTools/2.30.0-GCC-11.2.0
source activate $CONDAVAL

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

    cp $SCRIPTSOURCE"/genes_refseq.hg38fix.bed" $shortout/metagene_plots/
    cd $shortout/metagene_plots/

    echo " - Generating profile plots"
    #echo "python3 $SCRIPTSOURCE/plotHeatmap.py -b $bedlist -c $shortout/heatmap.config \
    #                        -l $namelist -o $shortout/metagene_plots -n 10 -s 30 \
    #                        -m $m -g 1 -t $THREADCOUNT -r $r"

    if [[ $m == "T" ]]; then
    python3 $SCRIPTSOURCE/plotHeatmap.py -b $bedlist -c $shortout/heatmap.config \
                            -l $namelist -o $shortout/metagene_plots/ -n 10 -s 30 \
                            -e $leftgenebp -y $genebodybp -k $rightgenebp \
                            -m -g 1 -t $THREADCOUNT -x $MAX -a $output #-r
    else
    python3 $SCRIPTSOURCE/plotHeatmap.py -b $bedlist -c $shortout/heatmap.config \
                            -l $namelist -o $shortout/metagene_plots/ -n 10 -s 30 \
                            -e $leftgenebp -y $genebodybp -k $rightgenebp \
                            -g 1 -t $THREADCOUNT -x $MAX -a $output #-r
    fi
    echo " - Profile plots made"
else
    echo "No major gene changes found, no profile plots made"
fi