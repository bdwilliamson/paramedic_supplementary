# The `stan` directory

This directory contains all of the Stan files necessary to replicate the numerical experiments and data analyses in "A multi-view model for relative and absolute microbial abundances" by Williamson, Hughes, and Willis. A brief description of each Stan file is as follows:
* `predict_qpcr_noncentered.stan`: The main, non-centered parameterization of the proposed efficiency-naive hierarchical model (i.e., the concentrations are modeled by a transformation of a standard normal random variable)
* `predict_qpcr_with_varying_efficiency_noncentered.stan`: The main, non-centered parameterization of the proposed varying-efficiency hierarchical model (i.e., the concentrations are modeled by a transformation of a standard normal random variable)
* `predict_qpcr.stan`: The main, centered parameterization of the proposed efficiency-naive hierarchical model (i.e., the concentrations are modeled by a normal random variable with non-zero mean and non-identity covariance matrix)
* `predict_qpcr_with_varying_efficiency.stan`: The main, centered parameterization of the proposed varying-efficiency hierarchical model (i.e., the concentrations are modeled by a normal random variable with non-zero mean and non-identity covariance matrix)
* `predict_qpcr_noncentered_vectorized.stan`: An attempt at vectorizing the model
* `predict_qpcr_noncentered_sigma_beta_*.stan`: Non-centered efficiency-naive models with differing, fixed values of `sigma_beta` (the default was `50`)
* `predict_qpcr_sigma_beta_*.stan`: Centered efficiency-naive models with differing, fixed values of `sigma_beta` (the default was `50`)
* `predict_qpcr_small_var.stan`, `predict_qpcr_wide_var.stan`: Centered efficiency-naive models with small (sqrt(10)) and large (sqrt(50)) variances on `beta`
* `predict_qpcr_with_varying_efficiency_noncentered_hard_e_constraint.stan`: the proposed varying-efficiency model with a hard constraint on the efficiencies
