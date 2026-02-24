data {
  // -----------------------------
  // Dimensions
  // -----------------------------
  int<lower=1> N;                 // number of observations (testable sample × peptide × HLA combinations)
  int<lower=1> S;                 // number of samples/patients (here: TP1 samples)
  int<lower=1> J_pair;            // number of unique peptide-HLA pairs
  int<lower=1> J_pep;             // number of unique peptides
  int<lower=1> J_allele;          // number of unique HLA alleles
  int<lower=1> K_age;             // dimension of age spline basis

  // -----------------------------
  // Binomial outcomes
  // -----------------------------
  array[N] int<lower=0> y;        // successes: number of CD8 T cells that are multimer+ for the given peptide-HLA pair
  array[N] int<lower=0> n;        // trials: total CD8 count for each observation

  // -----------------------------
  // Observation-level indices linking each observation to sample, peptide, HLA, and pair
  // (used for indirect indexing of sample-level covariates, e.g. sex_M[sample_of_obs[n_i]])
  // -----------------------------
  array[N] int<lower=1, upper=S> sample_of_obs;
  array[N] int<lower=1, upper=J_pep> pep_of_obs;
  array[N] int<lower=1, upper=J_allele> allele_of_obs;
  array[N] int<lower=1, upper=J_pair> pair_of_obs;

  // -----------------------------
  // Sample-level covariates (potential confounders)
  // -----------------------------
  array[S] int<lower=0, upper=1> severity;   // 0 = mild, 1 = severe
  array[S] int<lower=0, upper=1> sex_M;      // 0 = female, 1 = male
  matrix[S, K_age] B_age;                    // age spline basis (recommended centered/scaled in R)
}

parameters {
  // Hierarchical structure:
  // Sample-, peptide-, HLA-, and pair-level effects are modeled hierarchically
  // with estimated scale parameters (SDs). This induces partial pooling ("shrinkage"):
  // levels with limited data are pulled more toward the overall mean (here zero on the
  // log-odds scale), while well-supported levels are allowed larger deviations.
  
  // =============================
  // Global regression coefficients (shared across observations; "fixed effects")
  // =============================
  real b0;                    // global intercept (baseline log-odds)
  real b_sexM;                // male effect
  vector[K_age] b_age;        // spline coefficients for age
  real b_sev;                 // global severity shift (shared across all observations)

  // =============================
  // Patient/sample-level deviations (hierarchical; non-centered parameterization; "random effects")
  // Captures overall between-patient differences in baseline multimer+ level
  // (including biological and technical variation not otherwise modeled)
  // =============================
  real<lower=0> sd_samp;      // SD of sample random intercepts
  vector[S] z_samp;           // standard-normal latent variables

  // =============================
  // Hierarchical SDs for peptide / HLA / pair structure
  // =============================
  real<lower=0> sd_pep0;      // SD of peptide baseline effects (mild)
  real<lower=0> sd_pepSev;    // SD of peptide-specific severity shifts
  real<lower=0> sd_hla0;      // SD of HLA baseline effects (mild)
  real<lower=0> sd_hlaSev;    // SD of HLA-specific severity shifts
  real<lower=0> sd_pair0;     // SD of residual peptide-HLA pair effects (baseline only)

  // =============================
  // Standard-normal latent deviations for hierarchical effects 
  // Non-centered; used to construct group-level effects (“random effects”)
  // =============================
  vector[J_pep] z_pep0;       // peptide baseline effects
  vector[J_pep] z_pepSev;     // peptide severity-shift effects
  vector[J_allele] z_hla0;    // HLA baseline effects
  vector[J_allele] z_hlaSev;  // HLA severity-shift effects
  vector[J_pair] z_pair0;     // residual pair effects (no severity interaction)
}

transformed parameters {
  // Convert latent deviations + SDs (non-centered parameterization) into group-level effects.
  vector[S] samp_re = sd_samp * z_samp;
  
  vector[J_pep]    re_pep0   = sd_pep0   * z_pep0;      // peptide effect in mild
  vector[J_pep]    re_pepSev = sd_pepSev * z_pepSev;    // peptide-specific shift in severe
  vector[J_allele] re_hla0   = sd_hla0   * z_hla0;      // HLA effect in mild
  vector[J_allele] re_hlaSev = sd_hlaSev * z_hlaSev;    // HLA-specific shift in severe
  vector[J_pair]   re_pair0  = sd_pair0  * z_pair0;     // residual pair effect (baseline only)
}

model {
  // =============================
  // Priors: global regression coefficients ("fixed effects")
  // Weakly informative priors on log-odds scale
  // =============================
  b0     ~ normal(0, 1.5);
  b_sexM ~ normal(0, 0.7);
  b_sev  ~ normal(0, 0.7);
  b_age  ~ normal(0, 0.7);

  // =============================
  // Priors: hierarchical SDs
  // Half-normal due to lower=0 constraints
  // =============================
  sd_samp   ~ normal(0, 0.7);

  sd_pep0   ~ normal(0, 0.5);
  sd_pepSev ~ normal(0, 0.35);
  sd_hla0   ~ normal(0, 0.6);
  sd_hlaSev ~ normal(0, 0.4);
  sd_pair0  ~ normal(0, 0.35);

  // =============================
  // Priors: non-centered latent variables
  // =============================
  z_samp   ~ std_normal();
  z_pep0   ~ std_normal();
  z_pepSev ~ std_normal();
  z_hla0   ~ std_normal();
  z_hlaSev ~ std_normal();
  z_pair0  ~ std_normal();

  // =============================
  // Likelihood (binomial logistic regression)
  // =============================
  // Precompute sample-level contributions once per leapfrog evaluation
  // instead of recomputing inside the observation loop.
  vector[S] age_term = B_age * b_age;
  vector[S] sample_base;      // terms shared across all observations from sample s
  vector[S] sample_sev_term;  // sample-specific contribution of the global severity term (b_sev * severity[s])
  for (s in 1:S) {
    sample_base[s] =
      b0
      + b_sexM * sex_M[s]
      + age_term[s]
      + samp_re[s];
    sample_sev_term[s] = b_sev * severity[s];
  }

  for (n_i in 1:N) {
    int s = sample_of_obs[n_i];   // sample/patient index
    int p = pep_of_obs[n_i];      // peptide index
    int a = allele_of_obs[n_i];   // HLA allele index
    int j = pair_of_obs[n_i];     // peptide-HLA pair index
    int sev = severity[s];        // sample-level severity (0 mild / 1 severe)

    // Linear predictor (log-odds that a CD8 cell in this observation is multimer+ for this peptide-HLA pair)
    //
    // Components:
    // - sample-level baseline (intercept + sex + age + sample random intercept)
    // - global severity shift
    // - peptide and HLA baseline effects (mild)
    // - peptide- and HLA-specific severity shifts
    // - residual peptide-HLA pair effect (not severity-dependent)
    real eta =
      sample_base[s]
      + sample_sev_term[s]
      + re_pep0[p] + re_pepSev[p] * sev
      + re_hla0[a] + re_hlaSev[a] * sev
      + re_pair0[j];

    target += binomial_logit_lpmf(y[n_i] | n[n_i], eta);
  }
}
