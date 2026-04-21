#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
from pathlib import Path

# Default values
output_dir = "heatmap_out"
left_ext_bp = 3000
body_ext_bp = 5000
right_ext_bp = 4000
recalculate = False
not_spikein = False
num_processors = 1
bin_size = 25
smooth_length = 75
dpi = 600
rpm_max = -1
type_ = "gene"
sort_sample = 1
data_path = "/mnt/vstor/SOM_GENE_BEG33/ChIP_seq/hg38/DATA/"
tss_file = "/mnt/vstor/SOM_GENE_BEG33/software/compPASS/scripts/genes_refseq.hg38fix.bed"
plot_name="pol2_state_shift_profiles"

def run_command(cmd):
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)


def format_dir(directory):
    return directory if directory.endswith("/") else directory + "/"


def read_config_file(config_file):
    bws, bwlabels, colors, min_values, max_values = [], [], [], [], []
    has_max = False

    with open(config_file, 'r') as f:
        for line in f:
            tokens = line.strip().split("\t")
            if len(tokens) < 4:
                continue
            bw_file, bw_label, ext_len, color = tokens[:4]
            min_value, max_value = 0, 0
            if len(tokens) >= 5:
                has_max = True
                values = tokens[4].split(',')
                if len(values) == 1:
                    max_value = tokens[4]
                else:
                    min_value = values[0]
                    max_value = values[1]

            main, ext = Path(bw_file).stem, Path(bw_file).suffix[1:]

            if ext == "bam":
                bam_file = bw_file
                bw_file = f"{data_path}/{main}/{main}.{bin_size}.RPM.bw"
                if not Path(bw_file).exists():
                    cmd = f"bamCoverage -e {ext_len} -b {bam_file} -o {bw_file} --outFileFormat bigwig --smoothLength {smooth_length} --binSize {bin_size} --normalizeUsing CPM -p {num_processors}"
                    run_command(cmd)
                    subprocess.run(f"chgrp beg33 {bw_file}", shell=True)

                spike_in_file = f"{Path(bw_file).parent}/SpikeIn/spike_map_summary"
                if Path(spike_in_file).exists() and not not_spikein:
                    with open(spike_in_file, 'r') as spike_in_f:
                        spike_in_f.readline()  # skip header
                        summary_line = spike_in_f.readline()
                        summary_fields = summary_line.split('\t')
                        scale_factor_raw = 1000000 / float(summary_fields[2])
                        scale_factor = f"{scale_factor_raw:.4f}"
                        bw_file = f"{data_path}/{main}/{main}.{bin_size}.scaled.bw"
                        if not Path(bw_file).exists():
                            cmd = f"bamCoverage -e {ext_len} -b {bam_file} -o {bw_file} --outFileFormat bigwig --smoothLength {smooth_length} --binSize {bin_size} --scaleFactor {scale_factor} -p {num_processors}"
                            run_command(cmd)
                            subprocess.run(f"chgrp beg33 {bw_file}", shell=True)

            bws.append(bw_file)
            bwlabels.append(bw_label)
            colors.append(color)
            min_values.append(min_value)
            max_values.append(max_value)

    return bws, bwlabels, colors, min_values, max_values, has_max


def parse_args():
    parser = argparse.ArgumentParser(description="Generate profile plots")
    parser.add_argument('-b', required=True, help="BED/Gene list (comma-separated)")
    parser.add_argument('-c', required=True, help="Config file")
    parser.add_argument('-r', action='store_true', help="Recalculate the matrix file")
    parser.add_argument('-m', action='store_true', help="Do not use spikeIn data")
    parser.add_argument('-l', help="BED label list (comma-separated)")
    parser.add_argument('-o', default=output_dir, help="Output directory, default: 'heatmap_out'")
    parser.add_argument('-t', type=int, default=num_processors, help="Threads to be used, default: THREADCOUNT")
    parser.add_argument('-e', type=int, default=left_ext_bp, help="Extension length, default: 3000")
    parser.add_argument('-y', type=int, default=body_ext_bp, help="Body length, default: 5000")
    parser.add_argument('-k', type=int, default=right_ext_bp, help="Extension length, default: 4000")
    parser.add_argument('-n', type=int, default=bin_size, help="Bin size of bigwig file, default: 10")
    parser.add_argument('-s', type=int, default=smooth_length, help="Smooth length of bigwig file, default: 30")
    parser.add_argument('-g', type=int, default=sort_sample, help="Sort by sample, default: 1")
    parser.add_argument('-x', type=float, default=rpm_max, help="Y-max in RPM as a float value, overrides config max if set")
    parser.add_argument('-a', type=str, default=plot_name, help="heatmap name")
    return parser.parse_args()


def main():
    args = parse_args()

    bed_list = args.b
    config_file = args.c
    output_dir = format_dir(args.o)
    num_processors = args.t
    left_ext_bp = args.e
    body_ext_bp = args.y
    right_ext_bp = args.k
    bin_size = args.n # ignore, hard lock for now 
    smooth_length = args.s # ignore, hard lock for now 
    #dpi = args.d # ignore, hard lock for now 
    recalculate = args.r # var
    not_spikein = args.m # include as scalefactor input
    sort_sample = args.g 
    plot_max = args.x
    plot_name=args.a
    type_ = "gene"
    bed_labels = args.l # input in bash
    tss_file = f"{output_dir}/genes_refseq.hg38fix.bed"

    # Validate type
    #if type_ not in ["bed", "gene"]:
    #    sys.exit(f"Unknown type: {type_}. Please use 'bed' or 'gene'")

    # Default label generation if not provided
    if not bed_labels:
        bed_labels = ','.join([Path(f).stem for f in bed_list.split(',')])

    # Prepare output directory
    os.makedirs(output_dir, exist_ok=True)

    # force type gene
    type_="gene"

    # Handle gene type BED list
 
    gene_lists = bed_list.split(',')
    out_list = []
    for gene_list in gene_lists:
        out_name = f"{output_dir}/{Path(gene_list).stem}_tss.bed"
        out_list.append(out_name)
        subprocess.run(f"chgrp beg33 {gene_list}", shell=True)
        subprocess.run(f"dos2unix {gene_list}", shell=True)

        genes = {}
        with open(gene_list, 'r') as gene_file:
            for line in gene_file:
                genes[line.strip()] = ""

        with open(out_name, 'w') as tss_out_file, open(tss_file, 'r') as tss_file_obj:
            for line in tss_file_obj:
                fields = line.strip().split("\t")
                if fields[3] in genes:
                    tss_out_file.write(line)
    
    bed_list = ','.join(out_list)

    bed_list = bed_list.replace(",", " ")
    bed_labels = bed_labels.replace(",", " ")

    # Read configuration file
    bws, bwlabels, colors, min_values, max_values, has_max = read_config_file(config_file)

    bw_list = " ".join(bws)
    bwlabel_list = " ".join(bwlabels)
    color_list = " ".join(colors)

    ymax_option = ""
    zmax_option = ""
    if has_max:
        ymax_option = f"--yMin 0 --yMax {' '.join(map(str, max_values))}"
        zmax_option = f"--zMin {' '.join(map(str, min_values))} --zMax {' '.join(map(str, max_values))}"
    
    # treated as no default
    if plot_max > 0:
        ymax_option = f"--yMin 0 --yMax {plot_max}"

    matrix_zfile = f"{output_dir}matrix.gz"
    matrix_file = f"{output_dir}matrix.tab"
    heatmap_file = f"{output_dir}heatmap.pdf"
    heatmap_bed_file = f"{output_dir}heatmap.bed"
    profile_file = f"{output_dir}{plot_name}.pdf"

    if not Path(matrix_zfile).exists() or recalculate:
        ref_point = "TSS" if type_ == "gene" else "center"
        cmd = f"computeMatrix scale-regions -R {bed_list} -S {bw_list} --beforeRegionStartLength {left_ext_bp} --regionBodyLength {body_ext_bp} --afterRegionStartLength {right_ext_bp} --samplesLabel {bwlabel_list} -p {num_processors} -o {matrix_zfile} --outFileNameMatrix {matrix_file}"
        run_command(cmd)

    ref_label = "TSS"# if type_ == "gene" else "Peak"
    #cmd = f"plotHeatmap -m {matrix_zfile} --missingDataColor white --colorList {color_list} --regionsLabel {bed_labels} -o {heatmap_file} --refPointLabel {ref_label} --dpi {dpi} --sortUsingSamples {sort_sample} --outFileSortedRegions {heatmap_bed_file} --interpolationMethod nearest {ymax_option} {zmax_option}"
    #run_command(cmd)

    cmd = f"plotProfile --perGroup -m {matrix_zfile} --numPlotsPerRow 2 --regionsLabel {bed_labels} -o {profile_file} --refPointLabel {ref_label} {ymax_option}"
    run_command(cmd)

    subprocess.run(f"chgrp beg33 {output_dir}*", shell=True)


if __name__ == "__main__":
    main()
