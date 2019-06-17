data{
    int<lower=0> N;
    int<lower=0> q_obs;
    int<lower=0> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
}
parameters{
    vector[q] log_mu_tilde[N];
    vector[q] beta;
    vector<lower=0>[q] Sigma;
}
transformed parameters{
    vector<lower=0>[q] mu[N];
    for (j in 1:N){
        mu[j] = exp(beta + Sigma .* log_mu_tilde[j]);
    }
}
model {
    beta ~ normal(0, sqrt(196.0));
    Sigma ~ lognormal(0, sqrt(196.0));
    
    for (j in 1:N){
        log_mu_tilde[j] ~ normal(0, 1);
        V[j] ~ poisson(mu[j,1:q_obs]);
        W[j] ~ multinomial(mu[j]/sum(mu[j]));
    }
}

