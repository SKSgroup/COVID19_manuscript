# ================================================================================
# 02_run_stan_model.R
# Compile and fit the Stan mixed-effects logistic regression model.
#
# Input:
#   - data/processed/standata.rds
#
# Outputs (results/):
#   - stanfit.rds
#   - CmdStan CSV chain files
#
# Runtime:
#   - ~8-9 hours on MacBook Pro M3 Max (4 cores) in the author’s setup
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
