# comparison of Pol2 Activity State Shifts (compPASS)

The compPASS pipeline allows investigators to capture the landscape of major changes in RNA Polymerase II profile across the genome. The goal is to capture ratios beyond just the pausing ratio and use the identified genes to characterize the nature of the perturbation (with respect to Pol2) applied to the cells.

#<img width="522" height="393" alt="image" src="https://github.com/user-attachments/assets/b350ef1d-cf94-439a-92eb-1e692c07b35d" />

# Workflow
#<img width="510" height="386" alt="image" src="https://github.com/user-attachments/assets/501666e8-8568-4762-a4f6-39cc26d33f73" />

<img width="1100" height="663" alt="image" src="https://github.com/user-attachments/assets/5adf6ac5-2129-40e5-81c4-66d6c339c1fb" />

Inputs:
- Sequencing data in BAM format from two ChIP-seq or HiChIP style experiments (ideally control and case) with total RNA Pol2 as the target.

Outputs of compPASS
For each successful run of the compPASS pipeline, the user should receive:
•	A tab-separated values table which captures the RPM values for sample A’s regions of interest and calculated ratios, RPM values for sample B and calculated ratios and the log2 of sample B’s RPM / sample A’s values (representing the L2FC)
•	Scatterplots showing the L2FC of the following paired regions – promoter region and TSS region, TSS region and gene body region, and gene body region and TES region. These plots also highlight in red or blue respectively, the genes classified as stalling or entering, pausing or releasing and clogging or unloading.
•	Radar plots showing the number of genes assigned to each of the eight classifications
•	Text files listing genes in each class
•	Average plots of Pol2-associated RPM signal across the loci of genes in each class
•	Enrichment dot plots for the provided sets in the gene categories tables
