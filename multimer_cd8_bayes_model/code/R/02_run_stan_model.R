# ================================================================================
# 02_run_stan_model.R
# Compile and fit the Stan mixed-effects logistic regression model
#
# Input file:
#   - data/processed/standata.rds
#
# Main tasks:
#   - Load Stan data list (standata)
#   - Compile Stan model
#   - Run MCMC sampling with CmdStanR
#   - Save fitted model object and CmdStan CSV chain files
#   - Run CmdStan diagnostics
#
# Output files (written to results/):
#   - stanfit.rds
#   - CmdStan CSV chain files (one per chain)
#
# Runtime:
#   - Approximately 11-12 hours on a MacBook Pro M3 Max using 4 cores
# ================================================================================

library(cmdstanr)

# ================================================================================
# 1) Load Stan data list
# ================================================================================
sd = readRDS("data/processed/standata.rds")

# ================================================================================
# 2) Compile Stan model
# ================================================================================
mod = cmdstan_model("code/Stan/m3_binom_pair_samplere.stan")

# ================================================================================
# 3) Run MCMC sampling with CmdStanR
# ================================================================================
stanfit = mod$sample(
  data = sd,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.9,
  max_treedepth = 14,
  output_dir = "results/"
)

# ================================================================================
# 4) Save fitted model object
# ================================================================================
saveRDS(stanfit, "results/stanfit.rds")

# ================================================================================
# 5) Run CmdStan diagnostics
# ================================================================================
stanfit$cmdstan_diagnose()
