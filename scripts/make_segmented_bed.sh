#!/bin/bash
# Takes a BED format gene reference and converts into a segmented BED file
# Yaw Asante | yxa181@case.edu | January 7th, 2025 | Gryder Lab, CWRU
source ./config_file.sh

# Define the function to generate the BED file
make_bed_file() {
    local ref_file=$1
    local bed_file=$2
    
    # Open reference file and output BED file
    while IFS=$'\t' read -r chr start end strand gene overlap; do
        # Skip the header line in the reference file
        if [[ "$chr" == "chr"* ]]; then
            # Calculate positions based on strand
            if [[ "$strand" == "+" ]]; then
                pro_start_pos=$((start + pro_start))
                pro_end_pos=$((start + pro_end))
                tssr_start_pos=$((start + tssr_start))
                tssr_end_pos=$((start + tssr_end))
                gene_start_pos=$((start + gene_start))
                gene_end_pos=$((end + gene_end))
                tesr_start_pos=$((end + tesr_start))
                tesr_end_pos=$((end + tesr_end))
            else
                pro_start_pos=$((end - pro_end))
                pro_end_pos=$((end - pro_start))
                tssr_start_pos=$((end - tssr_end))
                tssr_end_pos=$((end - tssr_start))
                gene_start_pos=$((start - gene_end))
                gene_end_pos=$((start - gene_start))
                tesr_start_pos=$((start - tesr_end))
                tesr_end_pos=$((start - tesr_start))
            fi

            # Skip if gene end position is less than or equal to gene start position
            if [[ $gene_end_pos -le $gene_start_pos ]]; then
                continue
            fi

            # Write the results to the BED file
            echo -e "$chr\t$pro_start_pos\t$pro_end_pos\t$strand\t$gene\tPromoter\t$overlap" >> "$bed_file"
            echo -e "$chr\t$tssr_start_pos\t$tssr_end_pos\t$strand\t$gene\tTSSR\t$overlap" >> "$bed_file"
            echo -e "$chr\t$gene_start_pos\t$gene_end_pos\t$strand\t$gene\tGene Body\t$overlap" >> "$bed_file"
            echo -e "$chr\t$tesr_start_pos\t$tesr_end_pos\t$strand\t$gene\tTESR\t$overlap" >> "$bed_file"
        fi
    done < "$ref_file"
}

# Usage: ./make_bed_file.sh reference_file.bed output_bed_file.bed
make_bed_file "$1" "$2"
