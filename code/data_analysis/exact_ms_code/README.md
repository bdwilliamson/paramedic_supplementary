# The exact code used in the `paramedic` paper data analyses

This repository contains the exact code used to analyze the data presented in "A multi-view model for relative and absolute microbial abundances" by Williamson, Hughes, and Willis. As in all other cases, the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system. Additionally, the code will not run unless you have access to the study data.

## Primary results (main and jackknife)

These results can be reproduced by running
```{sh}
./shell/submit_all_ms_data_analyses.sh 1 4 3
```

Note that the command-line arguments correspond to whether or not to adjust for HIV-1 status, `alpha_sigma`, and `kappa_sigma`, respectively. Once you have the results of this analysis, you can run
```{sh}
./shell/load_qpcr_data_analysis.sh 127 13 55 1
```
to obtain the figures and tables.

## Sensitivity analyses

The first sensitivity analysis investigated not adjusting for HIV-1 status, and can be reproduced by running
```{sh}
./shell/submit_all_ms_data_analyses.sh 0 4 3
```
Once you have the results of this analysis, you can run
```{sh}
./shell/load_qpcr_data_analysis.sh 127 13 55 0
```
to obtain the figures and tables.

Finally, you can investigate the hyperparameters on the distribution of the efficiencies by running
```{sh}
./shell/submit_all_ms_data_analyses.sh 1 2 3
```
