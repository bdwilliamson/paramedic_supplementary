data{
    int<lower=0> N;
    int<lower=0> q_obs;
    int<lower=0> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
}
parameters{
    vector[q] log_mu[N];
    vector[q] beta;
    row_vector<lower=0>[q] Sigma;
}
transformed parameters{
    vector<lower=0>[q] mu[N];
    vector<lower=0>[q] p[N];
    mu = exp(log_mu);
    for (j in 1:N){
        p[j] = mu[j]/sum(mu[j]);
    }
}
model {
    beta ~ normal(0, sqrt(10.0));
    Sigma ~ lognormal(0, sqrt(10.0));
    
    for (j in 1:N){
        log_mu[j] ~ normal(beta, Sigma);
        V[j] ~ poisson(mu[j,1:q_obs]);
        W[j] ~ multinomial(p[j]);
    }
}

