# ================================================================================
# 03_compute_delta_p_from_posterior_draws.R
# Compute posterior average predictive comparisons (APCs) from Stan posterior draws
#
# Definition:
#   - delta-p = p_severe - p_mild
#   - APCs are obtained by averaging delta-p across observations within each group
#     (peptide, HLA allele, or peptide-HLA pair)
#
# Input files:
#   - data/processed/standata.rds
#   - data/processed/index_map.rds
#   - results/stanfit.rds
#
# Main tasks:
#   - Extract posterior draws from fitted Stan model
#   - Reconstruct linear predictors for severity = 0 and severity = 1
#   - Compute observation-level delta-p for each posterior draw
#   - Aggregate to posterior APCs for peptide, HLA, and peptide-HLA pair levels
#   - Convert outputs to tidy long-format tables
#
# Output files (written to results/):
#   - pep_post.rds
#   - hla_post.rds
#   - pair_post.rds
# ================================================================================

library(tidyverse)

standata  = readRDS("data/processed/standata.rds")
indexmap  = readRDS("data/processed/index_map.rds")
fit       = readRDS("results/stanfit.rds")

# Logistic / inverse-logit function
# (Equivalent to plogis(x); written explicitly for clarity)
logistic = function(x) 1 / (1 + exp(-x))

# =============================================================================
# 1) Extract posterior draws from fitted Stan model
#    Convert to plain data frame so parameters can be indexed easily
# =============================================================================

dr = as.data.frame(posterior::as_draws_df(fit$draws()))
D  = nrow(dr)  # number of posterior draws

# =============================================================================
# 2) Pull observation-level indices and covariates from standata
#    Each observation n is one patient/sample × peptide-HLA measurement
# =============================================================================

N = standata$N

# Indices (length N): which sample / peptide / allele / pair each observation belongs to
s    = standata$sample_of_obs
pep  = standata$pep_of_obs
hla  = standata$allele_of_obs
pair = standata$pair_of_obs

# Number of levels for each grouping
J_pep  = standata$J_pep
J_hla  = standata$J_allele
J_pair = standata$J_pair

# Observation-level covariates (looked up from sample-level values)
sex_s = standata$sex_M[s]
B     = standata$B_age[s, , drop = FALSE]   # N x K_age age spline basis at obs level

# Number of observations contributing to each APC (observation-weighted averaging)
n_pep_obs  = tabulate(pep,  nbins = J_pep)
n_hla_obs  = tabulate(hla,  nbins = J_hla)
n_pair_obs = tabulate(pair, nbins = J_pair)

# =============================================================================
# 3) Extract posterior parameter blocks (names match Stan model)
# =============================================================================

# Fixed effects
b0     = dr$b0
b_sev  = dr$b_sev
b_sexM = dr$b_sexM
b_age  = as.matrix(dr[, grep("^b_age\\[", names(dr)), drop = FALSE])

# Sample random intercept (non-centered)
sd_samp = dr$sd_samp
z_samp  = as.matrix(dr[, grep("^z_samp\\[", names(dr)), drop = FALSE])

# Hierarchical SDs
sd_pep0   = dr$sd_pep0
sd_pepSev = dr$sd_pepSev
sd_hla0   = dr$sd_hla0
sd_hlaSev = dr$sd_hlaSev
sd_pair0  = dr$sd_pair0

# Non-centered random effects
z_pep0   = as.matrix(dr[, grep("^z_pep0\\[", names(dr)), drop = FALSE])
z_pepSev = as.matrix(dr[, grep("^z_pepSev\\[", names(dr)), drop = FALSE])
z_hla0   = as.matrix(dr[, grep("^z_hla0\\[", names(dr)), drop = FALSE])
z_hlaSev = as.matrix(dr[, grep("^z_hlaSev\\[", names(dr)), drop = FALSE])
z_pair0  = as.matrix(dr[, grep("^z_pair0\\[", names(dr)), drop = FALSE])

# =============================================================================
# 4) Helper: rowsum with full index alignment (1..J)
#    rowsum() omits groups with zero entries; this restores full length.
# =============================================================================

rowsum_full = function(x, g, J) {
  rs = rowsum(x, g, reorder = FALSE)
  out = numeric(J)
  idx = as.integer(rownames(rs))
  out[idx] = as.numeric(rs)
  out
}

# =============================================================================
# 5) Compute APC posteriors in one pass over draws
#
#    For each posterior draw:
#      - reconstruct random effects
#      - compute eta0 and eta1 for each observation (severity set to 0/1)
#      - convert to probabilities p0 and p1
#      - compute dp_obs = p1 - p0
#      - aggregate dp_obs to peptide / HLA / pair APCs
# =============================================================================

pep_mat  = matrix(NA_real_, nrow = D, ncol = J_pep)
hla_mat  = matrix(NA_real_, nrow = D, ncol = J_hla)
pair_mat = matrix(NA_real_, nrow = D, ncol = J_pair)

for (d in seq_len(D)) {
  
  # --- Reconstruct random effects for this draw (non-centered parameterization)
  samp_re_obs = (sd_samp[d] * z_samp[d, ])[s]   # sample RE expanded to obs level
  
  re_pep0   = sd_pep0[d]   * z_pep0[d, ]
  re_pepSev = sd_pepSev[d] * z_pepSev[d, ]
  re_hla0   = sd_hla0[d]   * z_hla0[d, ]
  re_hlaSev = sd_hlaSev[d] * z_hlaSev[d, ]
  re_pair0  = sd_pair0[d]  * z_pair0[d, ]
  
  # --- Age contribution at observation level
  age_eta = as.numeric(B %*% b_age[d, ])
  
  # --- Linear predictor with severity forced to mild (sev = 0)
  # Matches the Stan model after replacing severity[s] by 0
  eta0 = b0[d] +
    b_sexM[d] * sex_s +
    age_eta +
    samp_re_obs +
    re_pep0[pep] +
    re_hla0[hla] +
    re_pair0[pair]
  
  # --- Linear predictor with severity forced to severe (sev = 1)
  # Adds global + peptide-specific + allele-specific severity shifts
  eta1 = eta0 +
    b_sev[d] +
    re_pepSev[pep] +
    re_hlaSev[hla]
  
  # --- Observation-level APC on probability scale
  # Δp = expected severe fraction - expected mild fraction
  dp_obs = logistic(eta1) - logistic(eta0)
  
  # --- Aggregate to group-specific APCs (observation-weighted means)
  pep_mat[d, ]  = rowsum_full(dp_obs, pep,  J_pep)  / n_pep_obs
  hla_mat[d, ]  = rowsum_full(dp_obs, hla,  J_hla)  / n_hla_obs
  pair_mat[d, ] = rowsum_full(dp_obs, pair, J_pair) / n_pair_obs
  
}

# =============================================================================
# 6) Metadata tables (labels for peptide / HLA / pair)
# =============================================================================

# Peptide and HLA labels from index map
pep_labels  = indexmap$peptide$label
hla_labels  = indexmap$allele$label
pair_labels = indexmap$pair$label  # expected format: "PEPTIDE | A01:01" (or similar)

# Pair metadata: keep pair label, and split into peptide + hla columns
pair_meta = tibble(
  pair_idx = seq_along(pair_labels),
  pair = pair_labels
) %>%
  tidyr::separate(
    pair,
    into   = c("peptide", "hla"),
    sep    = " \\| ",
    remove = FALSE
  )


# =============================================================================
# 7) Convert matrices to long tidy tables
#    (draw x group matrices -> long tibbles for plotting/summaries)
# =============================================================================

# Peptide posterior draws
pep_post = tibble(
  draw    = rep(seq_len(D), times = J_pep),
  pep_idx = rep(seq_len(J_pep), each = D),
  dp      = as.vector(pep_mat)
) %>%
  mutate(peptide = pep_labels[pep_idx]) %>%
  select(draw, peptide, dp)

# HLA posterior draws
hla_post = tibble(
  draw    = rep(seq_len(D), times = J_hla),
  hla_idx = rep(seq_len(J_hla), each = D),
  dp      = as.vector(hla_mat)
) %>%
  mutate(hla = hla_labels[hla_idx]) %>%
  select(draw, hla, dp)

# Pair posterior draws (with peptide + HLA columns attached)
pair_post = tibble(
  draw     = rep(seq_len(D), times = J_pair),
  pair_idx = rep(seq_len(J_pair), each = D),
  dp       = as.vector(pair_mat)
) %>%
  left_join(pair_meta, by = "pair_idx") %>%
  select(draw, pair, peptide, hla, dp)

# =============================================================================
# 8) Inspect and save
# =============================================================================

# Quick checks
print(pep_post)
print(hla_post)
print(pair_post)

# Save outputs
saveRDS(pep_post,  "results/pep_post.rds")
saveRDS(hla_post,  "results/hla_post.rds")
saveRDS(pair_post, "results/pair_post.rds")
