# ================================================================================
# 05_create_hla_figure.R
# Create posterior interval figure for HLA-level APCs (delta-p)
#
# Input files:
#   - results/hla_post.rds
#
# Main tasks:
#   - Summarize HLA-level posterior delta-p values
#   - Order HLA alleles by posterior median delta-p
#   - Build draw matrix for bayesplot::mcmc_intervals()
#   - Create interval plot (including pseudo-log x-axis version)
#
# Output files (written to figures/):
#   - hla_pseudolog.pdf
# ================================================================================

library(tidyverse)
library(bayesplot)

# ================================================================================
# 0) Load posterior APC draws (delta-p)
#    - hla_post : draw-wise HLA APCs (draw, hla, dp)
# ================================================================================
hla_post = readRDS("results/hla_post.rds")

# Region of practical equivalence (ROPE) around zero for delta-p
rope = 1e-4

# ================================================================================
# Helper: summarize posterior draws for delta-p within groups
# Returns means, quantiles, and posterior probabilities relative to 0 / ROPE
# ================================================================================
summarise_dp = function(data, group_var, rope = 1e-4) {
  data %>%
    group_by({{ group_var }}) %>%
    summarise(
      dp_mean       = mean(dp),
      dp_q10        = quantile(dp, probs = 0.10),
      dp_q50        = quantile(dp, probs = 0.50),
      dp_q90        = quantile(dp, probs = 0.90),
      p_gt_0        = mean(dp > 0),
      p_gt_rope     = mean(dp > rope),
      p_lt_neg_rope = mean(dp < -rope),
      p_in_rope     = mean(dplyr::between(dp, -rope, rope)),
      .groups = "drop"
    )
}

# ================================================================================
# 1) Summarize HLA-level posterior APCs
# ================================================================================
hla_summary = summarise_dp(hla_post, hla, rope = rope) %>% 
  arrange(desc(dp_q50), desc(p_gt_rope))

print(hla_summary)

# Order HLA alleles by posterior median delta-p (largest positive at top of plot)
hla_order = hla_summary$hla

# ================================================================================
# 2) Build draw matrix for bayesplot::mcmc_intervals()
#    bayesplot expects one column per parameter and one row per draw
# ================================================================================
df_plot = hla_post %>%
  pivot_wider(names_from = hla, values_from = dp) %>%
  select(draw, all_of(hla_order))

# ================================================================================
# 3) Create posterior interval plot for HLA alleles
# ================================================================================
# Linear-scale version (useful for checking absolute magnitudes)
p_hla = mcmc_intervals(df_plot %>% select(-draw), prob_outer = 0.9) +
  labs(x = expression(Delta * "p (severe - mild): change in pMHC+ fraction of CD8")) +
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2)

p_hla

# Pseudo-log x-axis to show both very small and large effects more clearly
p_hla_pseudolog = p_hla +
  scale_x_continuous(
    transform = scales::pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(-1, -0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.1, 1)
  )

p_hla_pseudolog

# ================================================================================
# 4) Save figure (PDF)
# ================================================================================
pdf("figures/hla_pseudolog.pdf", width = 7, height = 3, useDingbats = FALSE)
print(p_hla_pseudolog)
dev.off()
