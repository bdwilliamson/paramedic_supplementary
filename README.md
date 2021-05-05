# Supplementary materials for the `paramedic` paper

This repository contains code to reproduce the analyses in "A multi-view model for relative and absolute microbial abundances" by Williamson, Hughes, and Willis. All analyses were implemented in the freely available R programming language; specifically, version 3.4.3. This may cause difficulty with random number seed generation if you are using a different R version on your home computer. To get the same random number streams that we used, please include the code `RNGkind(sample.kind = "Rounding")` prior to any invocation of ` set.seed()`.

This README file provides an overview of the code available in the repository.  

## Code directory

This directory contains all R and bash scripts necessary to run the numerical experiments and replicate the data analyses.

## Stan directory

This directory contains all Stan code necessary to run the numerical experiments and replicate the data analyses.

-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/bdwilliamson/paramedic_supplementary/issues).
