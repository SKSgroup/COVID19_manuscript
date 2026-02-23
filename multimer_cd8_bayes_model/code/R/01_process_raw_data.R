# =========================================================================
# Read raw data from original excel files
# Extract and format relevant variables
# Encode age as cubic spline with 3 df
# Encode sex as binary (0=F, 1=M)
# Encode severity as binary (0=mild, 1=severe)
# Add indices for sample, peptide, hla, hla-peptide pair 
# (to index arrays of parameters)
# Create Stan data list and save to disk
# =========================================================================

library(tidyverse)
library(readxl)
library(janitor)
library(splines)

# =========================================================================
# Create data frame with selected patient covariates: severity, sex, age
# Encode:
#   sex as binary: sex_M (0=F, 1=M)
#   severity as binary severity_bin: 0=mild, 1=severe
#   age as cubic spline with 3 df
# =========================================================================
df_cov = read_excel("data/raw/Supplementary Tables_eBioMedicine.xlsx",
                        sheet=1, range = "B3:L76") %>% 
  clean_names() %>%  
  select(patient_id, severity, age, sex) %>% 
  rename(patient = patient_id) %>% 
  fill(severity, .direction = "down") %>%
  mutate(
    severity = case_when(
      str_detect(severity, regex("^\\s*severe\\b", ignore_case = TRUE)) ~ "severe",
      str_detect(severity, regex("^\\s*mild\\b",   ignore_case = TRUE)) ~ "mild"
    ),
    severity_bin = as.integer(severity=="severe"),
    age = age %>%
      str_trim() %>%            
      parse_integer(),            # -> integer (NA where not parseable)
    sex = str_trim(sex),
    sex_M = as.integer(sex == "M"),
    patient = as.integer(str_remove(patient, "^AP-")),  # "AP-01" -> 1
  ) %>%
  relocate(patient, severity, sex, sex_M, age)
B_age = ns(df_cov$age, df = 3)
B_age = scale(B_age, center = TRUE, scale = FALSE)
colnames(B_age) = c("age_s1", "age_s2", "age_s3")
df_cov = bind_cols(df_cov, as_tibble(B_age)) %>% 
  relocate(patient, severity, severity_bin, sex, sex_M, age)

df_cov
saveRDS(df_cov, "data/processed/df_cov.rds")

# =============================================================================================
# Create data frame with observed T-cell counts, one row per patient x peptide x HLA observation
#    cd8_count: total number CD8+ T-cells in sample (same for all observations from patient)
#    p_mhc_count: number of T-cells responding to this specific peptide-HLA pair
#    multimer_count: number of T-cells recognising any peptide-HLA pair in this patient sample
# =============================================================================================

# For groups of peptides recognised by same T-cell: 
# create mapping so we can replace group members by representative peptide
# (first member of group)
pep_groups = list(
  c("TTDPSFLGRY", "HTTDPSFLGRY", "TTDPSFLGRYM"),
  c("FTSDYYQLY",  "YFTSDYYQLY"),
  c("CTDDNALAYY", "CTDDNALAYYN", "TDDNALAYY"),
  c("VATSRTLSYY", "ATSRTLSYY")
)
pep_map = c()
for (g in pep_groups) {
  pep_map[g] = g[1]
}

df_counts = read_excel("data/raw/COVID_counts.xlsx") %>%
  clean_names() %>% 
  filter(timepoint == "01") %>% 
  select(patient, hla, peptide, cd8_count, est_freq_adjusted) %>% 
  mutate(patient = as.integer(patient)) %>% 
  mutate(cd8_count = as.integer(cd8_count),
         p_mhc_count = as.integer(round(cd8_count*est_freq_adjusted/100))
  ) %>%
  select(-est_freq_adjusted) %>% 
  mutate(
    peptide_raw = peptide,
    peptide = coalesce(unname(pep_map[peptide_raw]), peptide_raw)
  ) %>% 
  group_by(patient, hla, peptide) %>%
  summarise(
    cd8_count = first(cd8_count),                 # constant within patient+timepoint
    p_mhc_count = sum(p_mhc_count, na.rm = TRUE), # sum counts across merged peptides
    .groups = "drop"
  ) %>% group_by(patient) %>%
  mutate(multimer_count = sum(p_mhc_count, na.rm = TRUE)) %>%
  ungroup()  %>%
  mutate(
    hla = str_replace(as.character(hla), "^HLA-", ""),
    pair = paste(peptide, hla, sep = " | ")
  ) %>% 
  relocate(patient, hla, peptide, pair, cd8_count, multimer_count, p_mhc_count)

df_counts
saveRDS(df_counts, "data/processed/df_counts.rds")

# =============================================================================================
# Merge df_cov + df_counts to create df with all model-relevant data, one row per observation
# Also add IDs (consecutive indices) for sample, pep, hla, pair.
# IDs are assigned by factor level order (alphabetical); mappings are saved in index_map
# Used to create data structures for Stan, while ensuring alignment of indices etc.
# =============================================================================================
df_model = df_counts %>% 
  left_join(df_cov, by="patient") %>% 
  mutate(
    sample_id = as.integer(factor(patient)),
    pep_id    = as.integer(factor(peptide)),
    hla_id = as.integer(factor(hla)),
    pair_id   = as.integer(factor(pair))
  ) 

df_model
saveRDS(df_model, "data/processed/df_model.rds")

# =============================================================================================
# Create sample-level table (one row per sample/patient)
# =============================================================================================
df_samp = df_model %>%
  distinct(sample_id, patient, severity_bin, sex_M, age_s1, age_s2, age_s3) %>%
  arrange(sample_id)

df_samp
saveRDS(df_samp, "data/processed/df_samp.rds")

# ============================================================
# Optional sanity checks 
# ============================================================
RUN_CHECKS = TRUE

if (RUN_CHECKS) {
  stopifnot(!anyNA(df_model$sample_id))
  stopifnot(!anyNA(df_model$pep_id))
  stopifnot(!anyNA(df_model$hla_id))
  stopifnot(!anyNA(df_model$pair_id))
  
  stopifnot(all(df_model$p_mhc_count >= 0))
  stopifnot(all(df_model$cd8_count >= 0))
  stopifnot(all(df_model$p_mhc_count <= df_model$cd8_count))
  
  stopifnot(nrow(df_samp) == dplyr::n_distinct(df_model$sample_id))
}

# ============================================================
# Create Stan data list
# ============================================================

standata = list(
  N = nrow(df_model),
  S = nrow(df_samp),
  J_pair = n_distinct(df_model$pair_id),
  J_pep = n_distinct(df_model$pep_id),
  J_allele = n_distinct(df_model$hla_id),
  K_age = df_samp %>%
    select(starts_with("age_s")) %>%
    as.matrix() %>% 
    ncol(),
  
  y = as.integer(df_model$p_mhc_count),
  n = as.integer(df_model$cd8_count),
  
  sample_of_obs = df_model$sample_id,
  pep_of_obs = df_model$pep_id,
  allele_of_obs = df_model$hla_id,
  pair_of_obs = df_model$pair_id,
  
  severity = df_samp$severity_bin,
  sex_M = df_samp$sex_M,
  B_age = df_samp %>% 
    select(age_s1, age_s2, age_s3) %>% 
    as.matrix()
)

standata
saveRDS(standata, "data/processed/standata.rds")

# ============================================================
# Index maps for interpreting Stan indices in post processing
# ============================================================

index_map = list(
  peptide = df_model %>%
    distinct(pep_id, peptide) %>%
    arrange(pep_id) %>%
    transmute(idx = pep_id, label = peptide),
  
  allele = df_model %>%
    distinct(hla_id, hla) %>%
    arrange(hla_id) %>%
    transmute(idx = hla_id, label = hla),
  
  pair = df_model %>%
    distinct(pair_id, pair, pep_id, hla_id) %>%
    arrange(pair_id) %>%
    transmute(idx = pair_id, label = pair, pep_id = pep_id, hla_id = hla_id)
)

index_map
saveRDS(index_map, "data/processed/index_map.rds")
