# ================================================================================
# 04_create_peptide_hla_pair_figure.R
# Create posterior interval figure for selected peptide-HLA pair APCs (delta-p).
#
# Inputs:
#   - results/pep_post.rds
#   - results/pair_post.rds
#
# Output (figures/):
#   - peptide_hla_25pairs_13pos_7neg_pseudolog.pdf
# ================================================================================

library(tidyverse)
library(bayesplot)

# ================================================================================
# 0) Load posterior APC draws (delta-p)
#    - pep_post : draw-wise peptide APCs (draw, peptide, dp)
#    - pair_post: draw-wise peptide-HLA pair APCs (draw, pair, peptide, hla, dp)
# ================================================================================
pep_post  = readRDS("results/pep_post.rds")
pair_post = readRDS("results/pair_post.rds")

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
# 1) Summarize peptide-level posterior APCs and select peptides for pair plot
# ================================================================================
pep_summary = summarise_dp(pep_post, peptide, rope = rope)

# Top 13 peptides with strongest positive delta-p (and high posterior support)
pep_sel_pos = pep_summary %>%
  filter(p_gt_rope > 0.95) %>%
  arrange(desc(dp_q50), desc(p_gt_rope)) %>%
  slice_head(n = 13)

print(pep_sel_pos)

# Top 7 peptides with strongest negative delta-p (higher in mild than severe)
pep_sel_neg = pep_summary %>%
  filter(p_lt_neg_rope > 0.95) %>%
  arrange(dp_q50, desc(p_lt_neg_rope)) %>%   # most negative first for selection
  slice_head(n = 7) %>%
  arrange(desc(dp_q50), desc(p_lt_neg_rope)) # reverse so final combined plot order is easier

print(pep_sel_neg)

# Combined set of 20 peptides used to define the pair plot
pep_sel = bind_rows(pep_sel_pos, pep_sel_neg)

# ================================================================================
# 2) Summarize peptide-HLA pair posterior APCs for selected peptides
# ================================================================================
pair_summary = pair_post %>%
  semi_join(pep_sel %>% select(peptide), by = "peptide") %>%
  summarise_dp(pair, rope = rope) %>%
  arrange(desc(dp_q50))

print(pair_summary)

# Order pairs by posterior median delta-p (largest positive at top of plot)
pair_order = pair_summary$pair

# Sanity check: should be 25 pairs (20 peptides, 5 of them with two HLAs)
message("Selected peptide-HLA pairs: ", length(pair_order))

# ================================================================================
# 3) Build draw matrix for bayesplot::mcmc_intervals()
#    bayesplot expects one column per parameter/pair and one row per draw
# ================================================================================
df_plot = pair_post %>%
  filter(pair %in% pair_order) %>%
  select(draw, pair, dp) %>%
  mutate(pair = factor(pair, levels = pair_order)) %>%
  pivot_wider(names_from = pair, values_from = dp) %>%
  select(draw, all_of(pair_order))

# ================================================================================
# 4) Create posterior interval plot for selected peptide-HLA pairs
# ================================================================================
# Linear-scale version (useful for checking absolute magnitudes)
p_pair = mcmc_intervals(df_plot %>% select(-draw), prob_outer = 0.9) +
  labs(x = expression(Delta * "p (severe - mild): change in pMHC+ fraction of CD8")) +
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2)

p_pair

# Pseudo-log x-axis to show both very small and large effects more clearly
p_pair_pseudolog = p_pair +
  ggplot2::scale_x_continuous(
    transform = scales::pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(-1, -0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.1, 1)
  )

p_pair_pseudolog

# ================================================================================
# 5) Save figure (PDF)
# ================================================================================
pdf("figures/peptide_hla_25pairs_13pos_7neg_pseudolog.pdf", width = 7, height = 6, useDingbats = FALSE)
print(p_pair_pseudolog)
dev.off()
