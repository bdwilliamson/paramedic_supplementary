## Supplementary materials for the `paramedic` paper

This repository contains code to reproduce the analyses in "Estimation of microbial abundances from 16s data and qPCR abundances" by Williamson, Willis, and Hughes. All analyses were implemented in the freely available R programming language.

This README file provides an overview of the code available in the repository.  

-----

## Data directory

-----

## Code directory

We have separated our code further into two sub-directories based on the two main objectives of the manuscript:

1. Numerical experiments to evaluate the operating characteristics of our proposed method under varying data-generating mechanisms.
2. An analysis of 433 taxa sampled from the vaginal microbiome of 1213 women.

Within each sub-directory, we further subdivide the code into an R directory (hosting all of the R code for the analysis) and a shell directory (hosting all of the code for batch submission to a high-performance cluster computing environment). All analyses were performed on a Linux cluster using the Slurm batch scheduling system. If you use a difference batch scheduling system, the individual code files are flagged with the line where you can change batch variables. If you prefer to run the analyses locally, you may -- however, these analyses will then take a large amount of time.

-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/bdwilliamson/paramedic_supplementary/issues).