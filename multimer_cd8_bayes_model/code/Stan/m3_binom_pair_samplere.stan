data {
  int<lower=1> N;                 // observations
  int<lower=1> S;                 // samples (here: patients at TP1)
  int<lower=1> J_pair;
  int<lower=1> J_pep;
  int<lower=1> J_allele;
  int<lower=1> K_age;

  array[N] int<lower=0> y;
  array[N] int<lower=0> n;

  array[N] int<lower=1, upper=S> sample_of_obs;
  array[N] int<lower=1, upper=J_pep> pep_of_obs;
  array[N] int<lower=1, upper=J_allele> allele_of_obs;
  array[N] int<lower=1, upper=J_pair> pair_of_obs;

  array[S] int<lower=0, upper=1> severity;
  array[S] int<lower=0, upper=1> sex_M;
  matrix[S, K_age] B_age;         // recommend mean-centered in R
}

parameters {
  // fixed
  real b0;
  real b_sexM;
  vector[K_age] b_age;
  real b_sev;

  // sample random intercept
  real<lower=0> sd_samp;
  vector[S] z_samp;

  // hierarchical SDs
  real<lower=0> sd_pep0;
  real<lower=0> sd_pepSev;
  real<lower=0> sd_hla0;
  real<lower=0> sd_hlaSev;
  real<lower=0> sd_pair0;

  // non-centered REs
  vector[J_pep] z_pep0;
  vector[J_pep] z_pepSev;
  vector[J_allele] z_hla0;
  vector[J_allele] z_hlaSev;
  vector[J_pair] z_pair0;
}

transformed parameters {
  vector[S] samp_re = sd_samp * z_samp;

  vector[J_pep]    re_pep0   = sd_pep0   * z_pep0;
  vector[J_pep]    re_pepSev = sd_pepSev * z_pepSev;
  vector[J_allele] re_hla0   = sd_hla0   * z_hla0;
  vector[J_allele] re_hlaSev = sd_hlaSev * z_hlaSev;
  vector[J_pair]   re_pair0  = sd_pair0  * z_pair0;
}

model {
  // fixed effects: weakly informative, a bit tighter than before
  b0     ~ normal(0, 1.5);
  b_sexM ~ normal(0, 0.7);
  b_sev  ~ normal(0, 0.7);
  b_age  ~ normal(0, 0.7);

  // SD priors: tighten to help geometry; adjust if you *know* effects are larger
  sd_samp   ~ normal(0, 0.7);

  sd_pep0   ~ normal(0, 0.5);
  sd_pepSev ~ normal(0, 0.35);
  sd_hla0   ~ normal(0, 0.6);
  sd_hlaSev ~ normal(0, 0.4);
  sd_pair0  ~ normal(0, 0.35);

  // non-centered
  z_samp   ~ std_normal();
  z_pep0   ~ std_normal();
  z_pepSev ~ std_normal();
  z_hla0   ~ std_normal();
  z_hlaSev ~ std_normal();
  z_pair0  ~ std_normal();

  // likelihood
  for (n_i in 1:N) {
    int s = sample_of_obs[n_i];
    int p = pep_of_obs[n_i];
    int a = allele_of_obs[n_i];
    int j = pair_of_obs[n_i];
    int sev = severity[s];

    real eta =
      b0
      + b_sexM * sex_M[s]
      + dot_product(B_age[s], b_age)
      + samp_re[s]
      + b_sev * sev
      + re_pep0[p] + re_pepSev[p] * sev
      + re_hla0[a] + re_hlaSev[a] * sev
      + re_pair0[j];

    target += binomial_logit_lpmf(y[n_i] | n[n_i], eta);
  }
}

generated quantities {
  array[N] real log_lik;
  for (n_i in 1:N) {
    int s = sample_of_obs[n_i];
    int p = pep_of_obs[n_i];
    int a = allele_of_obs[n_i];
    int j = pair_of_obs[n_i];
    int sev = severity[s];

    real eta =
      b0
      + b_sexM * sex_M[s]
      + dot_product(B_age[s], b_age)
      + samp_re[s]
      + b_sev * sev
      + re_pep0[p] + re_pepSev[p] * sev
      + re_hla0[a] + re_hlaSev[a] * sev
      + re_pair0[j];

    log_lik[n_i] = binomial_logit_lpmf(y[n_i] | n[n_i], eta);
  }
}
