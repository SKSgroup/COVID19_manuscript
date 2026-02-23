data {
  int<lower=1> N;                 // number of observations (testable sample × peptide × HLA combinations)
  int<lower=1> S;                 // number of samples/patients (here: TP1 samples)
  int<lower=1> J_pair;            // number of unique peptide-HLA pairs
  int<lower=1> J_pep;             // number of unique peptides
  int<lower=1> J_allele;          // number of unique HLA alleles
  int<lower=1> K_age;             // number of age-basis columns (spline basis dimension)

  array[N] int<lower=0> y;        // successes: pMHC+ (multimer+) CD8 events for each observation
  array[N] int<lower=0> n;        // trials: total CD8 count for each observation

  // Observation -> group indexing
  // each observation points to the index of its sample, peptide, allele, and pair)
  array[N] int<lower=1, upper=S> sample_of_obs;
  array[N] int<lower=1, upper=J_pep> pep_of_obs;
  array[N] int<lower=1, upper=J_allele> allele_of_obs;
  array[N] int<lower=1, upper=J_pair> pair_of_obs;

  array[S] int<lower=0, upper=1> severity;   // 0 = mild, 1 = severe
  array[S] int<lower=0, upper=1> sex_M;      // 0 = female, 1 = male
  matrix[S, K_age] B_age;                    // age spline basis (centered/scaled in R)
}

parameters {
  // Fixed effects (population-level)
  real b0;                    // global intercept (baseline log-odds)
  real b_sexM;                // male effect
  vector[K_age] b_age;        // spline coefficients for age
  real b_sev;                 // global severity shift (shared across all observations)

  // Sample/patient random intercept (non-centered)
  // Captures overall between-patient differences in responsiveness
  real<lower=0> sd_samp;      // SD of sample random intercepts
  vector[S] z_samp;           // standard-normal latent variables

  // Hierarchical SDs for epitope/HLA structure
  real<lower=0> sd_pep0;      // SD of peptide effects in mild (baseline)
  real<lower=0> sd_pepSev;    // SD of peptide-specific severity shifts
  real<lower=0> sd_hla0;      // SD of HLA effects in mild (baseline)
  real<lower=0> sd_hlaSev;    // SD of HLA-specific severity shifts
  real<lower=0> sd_pair0;     // SD of residual peptide-HLA pair effects (baseline only)

  // Non-centered random effects (standard normal latents)
  vector[J_pep] z_pep0;       // peptide baseline effects
  vector[J_pep] z_pepSev;     // peptide severity-shift effects
  vector[J_allele] z_hla0;    // HLA baseline effects
  vector[J_allele] z_hlaSev;  // HLA severity-shift effects
  vector[J_pair] z_pair0;     // residual pair effects (no severity interaction)
}

transformed parameters {
  // Convert non-centered parameterization to actual random effects
  // (improves sampling geometry vs directly sampling centered random effects)
  vector[S] samp_re = sd_samp * z_samp;

  vector[J_pep]    re_pep0   = sd_pep0   * z_pep0;      // peptide effect in mild
  vector[J_pep]    re_pepSev = sd_pepSev * z_pepSev;    // peptide-specific shift in severe
  vector[J_allele] re_hla0   = sd_hla0   * z_hla0;      // HLA effect in mild
  vector[J_allele] re_hlaSev = sd_hlaSev * z_hlaSev;    // HLA-specific shift in severe
  vector[J_pair]   re_pair0  = sd_pair0  * z_pair0;     // residual peptide-HLA pair effect (baseline only)
}

model {
  // Priors: fixed effects
  // Weakly informative priors on log-odds scale
  b0     ~ normal(0, 1.5);
  b_sexM ~ normal(0, 0.7);
  b_sev  ~ normal(0, 0.7);
  b_age  ~ normal(0, 0.7);

  // Priors: hierarchical SDs
  // Half-normal due to lower=0 constraints
  sd_samp   ~ normal(0, 0.7);

  sd_pep0   ~ normal(0, 0.5);
  sd_pepSev ~ normal(0, 0.35);
  sd_hla0   ~ normal(0, 0.6);
  sd_hlaSev ~ normal(0, 0.4);
  sd_pair0  ~ normal(0, 0.35);

  // Priors: non-centered latent variables
  z_samp   ~ std_normal();
  z_pep0   ~ std_normal();
  z_pepSev ~ std_normal();
  z_hla0   ~ std_normal();
  z_hlaSev ~ std_normal();
  z_pair0  ~ std_normal();

  // Likelihood
  // Binomial logistic regression for pMHC+ fraction among CD8 cells
  for (n_i in 1:N) {
    int s = sample_of_obs[n_i];   // sample/patient index for observation n_i
    int p = pep_of_obs[n_i];      // peptide index
    int a = allele_of_obs[n_i];   // HLA allele index
    int j = pair_of_obs[n_i];     // peptide-HLA pair index
    int sev = severity[s];        // sample-level severity (0 mild / 1 severe)

    // Linear predictor (log-odds of pMHC+ among CD8)
    //
    // Components:
    // - global intercept + sample covariates (sex, age)
    // - sample random intercept
    // - global severity effect
    // - peptide and HLA baseline effects (mild)
    // - peptide- and HLA-specific severity shifts
    // - residual peptide-HLA pair effect (not severity-dependent)
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
