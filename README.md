# CutoffOpt

This repository contains R code to replicate the findings from: Benedetti et al. "A strategy to incorporate prior knowledge into correlation network cutoff selection", _Nature Communications_ (2020).

This code was created with R version 4.0.1 and Rstudio version 1.3.959.

## Glycomics Results

Sourcing the script _Main.R_ will reproduce the main glycomics results. The [preprocessed IgG glycomics data](https://doi.org/10.6084/m9.figshare.5335861) will automatically be downloaded from the figshare repository and the following files will be saved in the working directory (each file corresponds to the respective paper figure panel):

- **Figure3B.pdf** 
- **Figure3C.pdf** 
- **Figure4A.pdf** 
- **Figure4B.pdf** 
- **Figure4C.pdf** 

<ins>**Note:**</ins>
The current version of the code performs _nboot=10_ bootstraps to compute the confidence intervals of all plots, while in the paper we used _nboot=1000_. Increasing the number will substantially increase the runtime, which for _nboot=10_ is roughly 3 hours on a MacBook Pro (macOS version 10.15.1) with a 2.3 GHz Quad-Core Intel Core i5 processor and 16GB of RAM. Since in the paper figures we typically report the average across the bootstrapping results, results obtained with the default _nboot_ value given here will not be as smooth as the ones reported in the paper, but will be qualitatively the same.

## Metabolomics Results

The metabolomics datasets used in the paper is not publicly available due to study participant privacy policies. Data can be obtained upon request (see paper for details).

## Transcriptomics Results

Sourcing the file _TranscriptomicsResults.R_ will automatically generate the following files in the working directory (each file corresponds to the respective paper figure panel):

- **Figure6.pdf** 
- **Figure7A.pdf** 
- **Figure7B.html** 
- **Figure7C.html** 
- **Figure7D.html** 
- **Figure7E.html**  

- **SupplementaryData1.html** -> corresponding to **Supplementary Data 1**  (interactive version of Figure 6)
- **SupplementaryFigure8.pdf** -> corresponding to **Supplementary Figure 8** 

- **TranscriptomicsNetworkData.Rds** -> This data file is needed for the R-Markdown file _TranscriptomicsNetworks.Rmd_, which will create an interactive Shiny app to explore the optimization and network results for all significant transcriptomics pathways.

<ins>**Notes on precomputed data:**</ins>
In order to circumvent the long computation time necessary to regenerate the transcriptomics results from scratch, the current version of the _TranscriptomicsResults.R_ script loads precomputed versions of the biological references (STRING, CORUM and pathway annotations from Reactome), of the [TCGA PANCAN12 RNA-seq data](https://xenabrowser.net/datapages/?cohort=TCGA%20PANCAN12%20(PANCAN12)&addHub=https%3A%2F%2Flegacy.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) corrected for covariates, as well as precomputed versions of the bootstrapping results. All precomputed files are automatically downloaded from figshare upon code sourcing. This allows to generate the figures listed above in roughly !!!!!!!!XXX!!!!!!!!!! minutes on a MacBook Pro (macOS version 10.15.1) with a 2.3 GHz Quad-Core Intel Core i5 processor and 16GB of RAM. To generate all data from scratch with _nboot=100_ as in the paper, the user needs to set _use_precomputed=F_ at line !!!!!XX!!!!!!!! and source the script. This will take an estimated !!!!!!!XX!!!!!!!!! hours on the above-mentioned machine.
