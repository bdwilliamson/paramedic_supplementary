data{
    int<lower=0> N;
    int<lower=0> h;
    int<lower=0> k;
    int<lower=0> X[N,h];
    int<lower=0> Y[N,k];
}
parameters{
    vector[k] log_mu[N];
    vector[k] beta;
    cov_matrix[k] Sigma;
}
transformed parameters{
    vector<lower=0>[k] mu[N];
    mu = exp(log_mu);
}
model {
    vector[k] pY[N];
    for (i in 1:k){
        beta[i] ~ normal(1, sqrt(50.0));
        Sigma[i,i] ~ lognormal(0, sqrt(50.0));
    }
    for (j in 1:N){
        log_mu[j] ~ multi_normal(beta, Sigma);
        pY[j] = mu[j]/sum(mu[j]);
        X[j,] ~ poisson(mu[j,1:h]);
        Y[j,] ~ multinomial(pY[j]);
    }
}
