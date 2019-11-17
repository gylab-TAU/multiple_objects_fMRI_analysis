data {
  //int<lower=1> D; // number of predictors - 2
  int<lower=0> N; // number of samples  - ?
  int<lower=1> L; // number of subjects  - 15
  real y[N];
  int<lower=1,upper=L> ll[N];
  real x1[N];
  real x2[N];
}
parameters {
  real mu_1;
  real<lower=0> sigma_1;
  real mu_2;
  real<lower=0> sigma_2;
  real mu_interaction;
  real<lower=0> sigma_interaction;
  real beta_1_raw[L];
  real beta_2_raw[L];
  real beta_interaction_raw[L];
  real<lower=0> sigma_subj;
}
transformed parameters {
  real beta_1[L];
  real beta_2[L];
  real beta_interaction[L];
  real mu_sum;
  real mu_diff;
  for (l in 1:L){
    // implies: beta_1 ~ normal(mu_1, sigma_1)
    beta_1[l] = mu_1 + sigma_1 * beta_1_raw[l];
    beta_2[l] = mu_2 + sigma_2 * beta_2_raw[l];
    beta_interaction[l] = mu_interaction + sigma_interaction * beta_interaction_raw[l];
  }  
  
  mu_sum = mu_1+mu_2;
  mu_diff=(mu_1-mu_2)/(mu_1+mu_2);
}
model {
  // priors
  target += normal_lpdf(mu_1 | 0, 5);
  target += normal_lpdf(mu_2 | 0, 5);
  target += normal_lpdf(mu_interaction | 0, 5);
  target += cauchy_lpdf(sigma_1 | 0, 5) - cauchy_lccdf(0 | 0, 5);
  target += cauchy_lpdf(sigma_2 | 0, 5) - cauchy_lccdf(0 | 0, 5);
  target += cauchy_lpdf(sigma_interaction | 0, 5) - cauchy_lccdf(0 | 0, 5);
  target += cauchy_lpdf(sigma_subj | 0, 5) - cauchy_lccdf(0 | 0, 5);
  
  // calculating betas for each subjet:
  for (l in 1:L) {
    target += normal_lpdf (beta_1_raw[l] | 0, 1);
    target += normal_lpdf (beta_2_raw[l] | 0, 1);
    target += normal_lpdf (beta_interaction_raw[l] | 0, 1);
//    target += cauchy_lpdf(sigma_subj[l] | 0, 20) - cauchy_lccdf(0 | 0, 20);
  }
  for (n in 1:N)
    target += normal_lpdf (y[n] | (beta_1[ll[n]]*x1[n] + beta_2[ll[n]]*x2[n] + beta_interaction[ll[n]]*x1[n]*x2[n]), sigma_subj);
}
