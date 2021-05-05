# The `code` directory

This directory contains all of the R and bash scripts necessary to replicate the numerical experiments and data analyses in "A multi-view model for relative and absolute microbial abundances" by Williamson, Hughes, and Willis.

We have separated our code further into two sub-directories based on the two main objectives of the manuscript:

1. Numerical experiments to evaluate the operating characteristics of our proposed method under varying data-generating mechanisms.
2. An analysis of 127 taxa sampled from the vaginal microbiome of 55 women.

Within each sub-directory, we further subdivide the code into an R directory (hosting all of the R code for the analysis) and a shell directory (hosting all of the code for batch submission to a high-performance cluster computing environment). All analyses were performed on a Linux cluster using the Slurm batch scheduling system; each individual analysis was run on a cluster node with at least 4 cores and 16GB of memory (though each analysis may have been allocated less memory at run-time by the batch scheduler). If you use a difference batch scheduling system, the individual code files are flagged with the line where you can change batch variables. If you prefer to run the analyses locally, you may -- however, these analyses will then take a large amount of time.
