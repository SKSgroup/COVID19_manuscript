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

# Load Stan data list
sd = readRDS("data/processed/standata.rds")

# Compile Stan model
mod = cmdstan_model("code/Stan/m3_binom_pair_samplere.stan")

# Fit Stan model, save fit object and CSV files with draws to disk
# Note: this takes about 2 hours on a MacBook Pro M3 Max using 4 cores
stanfit = mod$sample(
  data = sd,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 15,
  output_dir = "results/"
)
saveRDS(stanfit, "results/stanfit.rds")

# Run diagnostics to check fit converged without issues
stanfit$cmdstan_diagnose()
