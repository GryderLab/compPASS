# compPASS: Comparison of Pol2 Activity State Shifts

The **compPASS** pipeline captures the landscape of major changes in RNA Polymerase II (Pol2) profile across the genome. By moving beyond the traditional pausing ratio and quantifying multiple region-based ratios, compPASS classifies genes into distinct activity states and helps characterize the nature of a perturbation applied to cells with respect to Pol2 behavior.

<img width="522" height="393" alt="compPASS overview" src="https://github.com/user-attachments/assets/b350ef1d-cf94-439a-92eb-1e692c07b35d" />

---

## Workflow

<img width="1016" height="263" alt="compPASS workflow" src="https://github.com/user-attachments/assets/ee3cb890-558e-4c67-97dd-d14cdf7c385a" />

---

## Inputs

compPASS requires sequencing data in BAM format from two ChIP-seq or HiChIP-style experiments — ideally a control and a case — with total RNA Pol2 as the target.

For each sample, the input folder must follow this structure:

```
SAMPLE_NAME/
├── SAMPLE_NAME.bam
├── SAMPLE_NAME.bam.bai
└── SAMPLE_NAME.flagstat.txt
```

If spike-in reads are used, include an additional file at `SAMPLE_NAME/SpikeIn/spike_map_summary` containing one header line and one data line in the following tab-separated format:

```
total_reads    human_reads    orthogonal_reads
```

---

## Installation

Create the conda environment from the provided YAML file:

```bash
conda env create -f compPASS_env.yaml
```

Once installation is complete, organize your input data into the folder structure described above under the `input/` directory (or any input folder of your choice).

---

## Outputs

Each successful run of the compPASS pipeline produces the following outputs:

- **Results table** — A tab-separated values (TSV) file containing RPM values for sample A's regions of interest along with calculated ratios, the equivalent RPM values and ratios for sample B, and the log2 fold change (sample B RPM / sample A RPM) for each region.
- **Scatterplots** — Plots of L2FC values for paired regions (promoter vs. TSS, TSS vs. gene body, and gene body vs. TES). Genes are highlighted in red or blue according to their classification as stalling/entering, pausing/releasing, or clogging/unloading.
- **Radar plots** — Visualizations showing the number of genes assigned to each of the eight activity state classifications.
- **Gene lists** — Text files listing the genes assigned to each class.
- **Average signal plots** — Plots of Pol2-associated RPM signal averaged across the loci of genes in each class.
- **Enrichment dot plots** — Enrichment visualizations for the user-provided gene sets across the compPASS gene categories.
