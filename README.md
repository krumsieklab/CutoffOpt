# CutoffOpt

This repository contains the R code to replicate the findings reported in the paper Benedetti et al. "A strategy to incorporate prior knowledge into correlation network cutoff selection", _Nature Communications_ (2020).

This code was created with R version 4.0.1 and Rstudio version 1.3.959.

## Glycomics Results

Sourcing the script _Main.R_ will reproduce the main glycomics results. The [preprocessed IgG glycomics data](https://doi.org/10.6084/m9.figshare.5335861) will be automatically downloaded from the figshare repository and the following files will be saved in the working directory:

- **Figure3B.pdf** -> corresponding to the GeneNet curve in **Figure 3B** in the paper 
- **Figure3C.pdf** -> corresponding to **Figure 3C** in the paper 
- **Figure4A.pdf** -> corresponding to **Figure 4A** in the paper 
- **Figure4B.pdf** -> corresponding to **Figure 4B** in the paper 
- **Figure4C.pdf** -> corresponding to **Figure 4C** in the paper 

<ins>**Note:**</ins>
The current version of the code performs _nboot=10_ bootstrapping to compute the confidence intervals of all all plots, while in the paper we used _nboot=1000_. Increasing the number of bootstrapping will substantially increase the runtime, which for _nboot=10_ is roughly 3 hours on a MacBook Pro (macOS version 10.15.1) with a 2.3 GHz Quad-Core Intel Core i5 processor and 16GB of RAM. Since in the paper figures we typically report the average across the bootstrapping results, results obtained with the default _nboot_ value will not be as smooth as the ones reported in the paper due to the low bootstrapping samples, but will be qualitatively equivalent.

## Metabolomics Results

The metabolomics datasets used in the paper is not publicly available due to patient privacy policies. Data access can be obtained upon reasonable request (see paper for details).

## Transcriptomics Results

Sourcing the file _TranscriptomicsResults.R_ will automatically generate the following files in the working directory:

- **Figure6.pdf** -> corresponding to **Figure 6** in the paper 
- **Figure7A.pdf** -> corresponding to **Figure 7A** in the paper 
- **Figure7B.html** -> corresponding to **Figure 7B** in the paper 
- **Figure7C.html** -> corresponding to **Figure 7C** in the paper 
- **Figure7D.html** -> corresponding to **Figure 7D** in the paper 
- **Figure7E.html** -> corresponding to **Figure 7E** in the paper 

- **SupplementaryData1.html** -> corresponding to **Supplementary Data 1** in the paper (interactive version of Figure 6)
- **SupplementaryFigure8.pdf** -> corresponding to **Supplementary Figure 8** in the paper 

- **TranscriptomicsNetworkData.Rds** -> this data file is needed for the R-Markdown file _TranscriptomicsNetworks.Rmd_, which will create an interactive shiny app to explore the optimization and network results for all significant transcriptomics pathways.

<ins>**Disclaimer:**</ins>
In order to spare the user the long computation time necessary to regenerate from scratch the transcriptomics results, the current version of the _TranscriptomicsResults.R_ script loads precomputed versions of the biological references (STRING, CORUM and pathway annotations from Reactome), of the [TCGA PANCAN12 RNA-seq data](https://xenabrowser.net/datapages/?cohort=TCGA%20PANCAN12%20(PANCAN12)&addHub=https%3A%2F%2Flegacy.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) corrected for covariates, as well as of the bootstrapping results. All the precomputed file are automatically downloaded from figshare upon code sourcing. This allows to generate the figures listed above in roughly XXX minutes on a MacBook Pro (macOS version 10.15.1) with a 2.3 GHz Quad-Core Intel Core i5 processor and 16GB of RAM. In order to generate all data from scratch with _nboot=100_ as in the paper, the user needs to set _use_precomputed=F_ at line XX and source the script. This will take an estimated XX hours on the above-mentioned machine.