The comparison of Pol2 Activity State Shifts (compPASS) pipeline allows investigators to capture the landscape of major changes in RNA Polymerase II profile across the genome. The goal is to capture ratios beyond just the pausing ratio and use the identified genes to characterize the nature of the perturbation (with respect to Pol2) applied to the cells.

<img width="522" height="393" alt="image" src="https://github.com/user-attachments/assets/b350ef1d-cf94-439a-92eb-1e692c07b35d" />


<img width="521" height="387" alt="image" src="https://github.com/user-attachments/assets/344df6dc-3cca-4900-913f-bb5caa191e44" />

Inputs:
- Sequencing data in BAM format from two ChIP-seq or HiChIP style experiments (ideally control and case) with total RNA Pol2 as the target.

Outputs: 
- A comparison table of Pol2 reads per million across gene proximal sites for both case in control
- Plots showcasing the difference in read density as well as Pol2 ratios (cumulative distribution, enrichment for gene categories and meta-gene average plots)
- Gene lists for all experiencing major Pol2 activity state shift
