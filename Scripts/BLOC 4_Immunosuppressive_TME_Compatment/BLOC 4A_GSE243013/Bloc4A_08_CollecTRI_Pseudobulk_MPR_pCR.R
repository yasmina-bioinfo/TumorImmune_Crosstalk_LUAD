#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 08: CollecTRI pseudobulk MPR vs pCR, TAMs
# Purpose: characterize pCR on its own terms relative to MPR,
# using the same Top 20 discriminant TFs already established
# (non-MPR vs MPR, Script 07), as a screening step to flag where
# pCR diverges from MPR (or does not) at the TF activity level.
# Input:  Objects/Bloc4A_07_tf_acts_raw.rds
#         Objects/Bloc4A_07_seu_TAMs_TF.rds
# Output: Results/Tables/Bloc4A_08_CollecTRI_pseudobulk_MPR_vs_pCR.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load raw TF activity scores (already computed, Script 07)
tf_acts <- readRDS(file.path(DATA_DIR, "Objects/Bloc4A_07_tf_acts_raw.rds"))

# 2) Same Top 20 TFs already established as discriminant (Script 07)
top20_tfs <- c("RFXAP","HSF1","RELA","RFXANK","NFKB","STAT1","NFKB1","REL",
               "RFX5","JUN","IRF1","CIITA","HIF1A","HSF2","STAT3","CREB1",
               "AP1","DDIT3","SPI1","EGR1")

tf_scores <- tf_acts %>%
  filter(statistic == "ulm", source %in% top20_tfs) %>%
  select(source, condition, score)

# 3) Load metadata, restrict to MPR + pCR
seu <- readRDS(file.path(DATA_DIR, "Objects/Bloc4A_07_seu_TAMs_TF.rds"))
meta <- seu@meta.data %>%
  tibble::rownames_to_column("condition") %>%
  select(condition, tam_short, sampleID, `pathological_response`) %>%
  rename(PathResponse = `pathological_response`) %>%
  filter(PathResponse %in% c("MPR", "pCR"))

# 4) Merge, compute patient-level pseudobulk mean per TF x subtype x patient
merged <- tf_scores %>% inner_join(meta, by = "condition")

pseudobulk <- merged %>%
  group_by(source, tam_short, sampleID, PathResponse) %>%
  summarise(mean_score = mean(score), .groups = "drop")

# 5) Wilcoxon MPR vs pCR, within each TF x TAM subtype
results <- pseudobulk %>%
  group_by(source, tam_short) %>%
  summarise(
    n_MPR = sum(PathResponse == "MPR"),
    n_pCR = sum(PathResponse == "pCR"),
    median_MPR = median(mean_score[PathResponse == "MPR"]),
    median_pCR = median(mean_score[PathResponse == "pCR"]),
    p_value = if (n_MPR >= 3 & n_pCR >= 3) {
      wilcox.test(mean_score[PathResponse == "MPR"], mean_score[PathResponse == "pCR"])$p.value
    } else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(tam_short) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(tam_short, p_adj)

fwrite(results, file.path(OUT_TAB, "Bloc4A_08_CollecTRI_pseudobulk_MPR_vs_pCR.csv"))
message("Saved: Bloc4A_08_CollecTRI_pseudobulk_MPR_vs_pCR.csv (",
        sum(results$p_adj < 0.05, na.rm = TRUE), " significant of ", nrow(results), " tested)")
print(results)