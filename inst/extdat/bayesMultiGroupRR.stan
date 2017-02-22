data {
  int<lower=0> n; // total number of observations
  int<lower=0> n_groups; // number of groups
  int<lower=0> n_loci; //
  int<lower=0> group[n]; // group membership
  vector[n] y; // phenotypes
  matrix[n, n_loci] X; // genotypes
}
parameters {
  real mu; // overall mean
  vector[n_groups] beta; // group-specific mean
  matrix[n_loci, n_groups] alpha_raw; // group-specific locus effects
  vector<lower=0>[n_loci] alpha_mu;
  vector<lower=0>[n_groups] sigma_alpha; // group-specific locus effect standard dev
  cholesky_factor_corr[n_groups] L; // Left Cholesky factor of correlaton matrix;
  real sigma_beta; // group mean standard dev
  vector<lower=0>[n_groups] sigma_e; // group-specific errors standard dev
  real<lower=0> sigma_alpha_mu; // overall locus effect standard dev
}
transformed parameters {
  matrix[n_loci, n_groups] alpha; // group-specific locus effects
  vector[n] theta; // The linear predictor
  //vector[n_groups] zero_vector;
  //Sigma = quad_form_diag(Omega, sigma_alpha);
  //for (k in 1:n_groups) zero_vector[k] = 0;

  for (i in 1:n_loci)
    alpha[i] = alpha_mu[i] + alpha_raw[i] * (L' * diag_matrix(sigma_alpha));

  //alpha = alpha_raw * (L' * diag_matrix(sigma_alpha));

  for (i in 1:n) {
    int grp = group[i];
    theta[i] = mu + beta[grp] + row(X, i) * col(alpha, grp);
  }

}
model {
  mu ~ normal(0, 1000);
  beta ~ normal(mu, sigma_beta);
  alpha_mu ~ normal(0, sigma_alpha_mu);
  sigma_beta ~ cauchy(0, 5);
  sigma_e ~ cauchy(0, 5);
  sigma_alpha ~ cauchy(0, 5);
  sigma_alpha_mu ~ cauchy(0, 5);
  L ~ lkj_corr_cholesky(1); 

  for (i in 1:n_loci)
    row(alpha_raw, i) ~ normal(0, 1);
    //row(alpha, i) ~ multi_normal(zero_vector, Sigma);

  for (i in 1:n)
    y[i] ~ normal(theta[i], sigma_e[group[i]]);
}
generated quantities {
  corr_matrix[n_groups] Omega; // Correlation matrix of locus effects
  cov_matrix[n_groups] Sigma; // Covariance matrix of locus effects
  vector<lower=0, upper=1>[n_groups] h2; // heritability
  Omega = tcrossprod(L);
  Sigma = quad_form_diag(Omega, sigma_alpha);
  for (i in 1:n_groups) h2[i] = sigma_alpha[i] / (sigma_alpha[i] + sigma_e[i]);
}