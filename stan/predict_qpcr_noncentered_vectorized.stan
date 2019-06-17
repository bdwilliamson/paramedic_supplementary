data{
    int<lower=0> N;
    int<lower=0> q_obs;
    int<lower=0> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
}
parameters{
    matrix[N,q] log_mu_tilde;
    vector[q] beta;
    row_vector<lower=0>[q] Sigma;
}
transformed parameters{
    matrix<lower=0>[N,q] mu;
    matrix<lower=0>[q,q] Sigma_tilde;
    matrix[q,N] beta_tilde;

    Sigma_tilde = rep_matrix(Sigma, q);
    beta_tilde = rep_matrix(beta, N);
    mu = exp(beta_tilde + (Sigma_tilde*log_mu_tilde'))';
}
model {
    beta ~ normal(0, sqrt(50.0));
    Sigma ~ lognormal(0, sqrt(50.0));
    
    for (j in 1:N){
        log_mu_tilde[j,] ~ normal(0, 1);
        V[j] ~ poisson(mu[j,1:q_obs]);
        W[j] ~ multinomial(to_vector(mu[j]/sum(mu[j])));
    }
}

