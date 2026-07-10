#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 14: CollecTRI pseudobulk MPR vs pCR, CD8
# Purpose: test whether the 20 CollecTRI TFs already established
# as discriminant between MPR and non-MPR (Script 12b) continue
# their trend, reverse, or show no relationship when comparing
# MPR to pCR specifically, at the patient level.
# Input:  Objects/Bloc3_12_tf_acts_raw.rds (raw per-cell TF scores)
#         Bloc3_06_seu_final_annotated.rds (metadata: sampleID, pathological_response)
# Output: Results/Tables/BLOC3/Bloc3_14_CollecTRI_pseudobulk_MPR_vs_pCR.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load raw TF activity scores (already computed, Script 12b)
tf_acts <- readRDS(file.path(DATA_DIR, "Objects/Bloc3_12_tf_acts_raw.rds"))

# 2) Same Top 20 TFs already established as discriminant (Script 12b)
top20_tfs <- c("MTF1","ELK4","HOPX","HSF4","HSF2","FOXA1","HSF1","NFAT5",
               "NFYB","NFYC","FOXF2","KLF4","MYC","MLX","RFX5","RFXANK",
               "RFXAP","MLXIP","DAXX","FOXO3")

tf_scores <- tf_acts %>%
  filter(statistic == "ulm", source %in% top20_tfs) %>%
  select(source, condition, score)

# 3) Load metadata (sampleID, pathological_response, cell_type) and restrict to MPR + pCR
seu <- readRDS(file.path(DATA_DIR, "Objects/Bloc3_12b_seu_CD8_TF.rds"))
meta <- seu@meta.data %>%
  tibble::rownames_to_column("condition") %>%
  select(condition, functional.cluster, sampleID, pathological_response) %>%
  filter(functional.cluster %in% c("CD8.TEX", "CD8.TPEX"),
         pathological_response %in% c("MPR", "pCR"))

# 4) Merge, compute patient-level pseudobulk mean per TF x state x patient
merged <- tf_scores %>%
  inner_join(meta, by = "condition")

pseudobulk <- merged %>%
  group_by(source, functional.cluster, sampleID, pathological_response) %>%
  summarise(mean_score = mean(score), .groups = "drop")

# 5) Wilcoxon MPR vs pCR, within each TF x CD8 state
results <- pseudobulk %>%
  group_by(source, functional.cluster) %>%
  summarise(
    n_MPR = sum(pathological_response == "MPR"),
    n_pCR = sum(pathological_response == "pCR"),
    median_MPR = median(mean_score[pathological_response == "MPR"]),
    median_pCR = median(mean_score[pathological_response == "pCR"]),
    p_value = if (n_MPR >= 3 & n_pCR >= 3) {
      wilcox.test(mean_score[pathological_response == "MPR"], mean_score[pathological_response == "pCR"])$p.value
    } else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(functional.cluster) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(functional.cluster, p_adj)

fwrite(results, file.path(OUT_TAB, "Bloc3_14_CollecTRI_pseudobulk_MPR_vs_pCR.csv"))
message("Saved: Bloc3_14_CollecTRI_pseudobulk_MPR_vs_pCR.csv (",
        sum(results$p_adj < 0.05, na.rm = TRUE), " significant of ", nrow(results), " tested)")
print(results)