library(tidyverse)
library(cmdstanr)

# Load Stan data list
sd = readRDS("data/processed/standata.rds")

# Compile Stan model
mod = cmdstan_model("code/Stan/m3_binom_pair_samplere.stan")

# Fit Stan model, save fit object and CSV files with draws to disk
# Note: this takes about 8 hours on a MacBook Pro M3 Max using 4 cores
stanfit = mod$sample(
  data = sd,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.8,
  max_treedepth = 10,
  output_dir = "results/"
)
saveRDS(stanfit, "results/stanfit.rds")

# Run diagnostics to check fit converged without issues
stanfit$cmdstan_diagnose()