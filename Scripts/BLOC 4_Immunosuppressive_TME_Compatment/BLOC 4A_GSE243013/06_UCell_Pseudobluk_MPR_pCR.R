#!/usr/bin/env Rscript
# ============================================================
# GSE243013  Bloc4A Script 06: Pseudobulk UCell scores MPR vs pCR
# Patient-level comparison to test MPR as intermediate state hypothesis
# NOTE: non-MPR excluded, focus on MPR (n=10) vs pCR (n=11)
# to formally test whether partial responders are biologically distinct
# from complete responders at individual patient level
# Pseudobulk: mean UCell score per patient per signature
# Wilcoxon rank-sum test MPR vs pCR (n=10 vs n=11)
# References: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)
# Input:  Objects/Bloc4A_05_seu_TAMs_UCell.rds
# Output: Results/Figures/BLOC4A_TAMs/Bloc4A_pseudobulk_MPR_pCR_*.png
#         Results/Tables/BLOC4A/Bloc4A_pseudobulk_MPR_pCR_scores.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggpubr)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4A_05_seu_TAMs_UCell.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC4A")

# 1) Load UCell object
message("Loading UCell TAMs object...")
seu_TAM <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_TAM))

# 2) Define score columns
score_cols <- c("M2_immunosuppressive_UCell", "M1_inflammatory_UCell",
                "SPP1_signature_UCell", "IFN_response_UCell")

# 3) Subset MPR and pCR only
# NOTE: non-MPR excluded, focus on MPR vs pCR to test intermediate state hypothesis
# n=10 MPR patients, n=11 pCR patients
message("Subsetting MPR and pCR...")
seu_sub <- subset(seu_TAM,
                  subset = pathological_response %in% c("MPR", "pCR"))
message("Cells MPR + pCR: ", ncol(seu_sub))
print(table(seu_sub$pathological_response))

# 4) Pseudobulk : mean UCell score per patient per signature
# NOTE: pseudobulk aggregates cell-level scores to patient-level means
# This avoids pseudo-replication and enables proper patient-level statistics
message("Computing pseudobulk scores...")
pseudobulk <- seu_sub@meta.data %>%
  filter(!is.na(sampleID), !is.na(pathological_response)) %>%
  select(sampleID, pathological_response, all_of(score_cols)) %>%
  group_by(sampleID, pathological_response) %>%
  summarise(across(all_of(score_cols), mean, .names = "{.col}"),
            n_cells = n(),
            .groups = "drop")

message("Patients per group:")
print(table(pseudobulk$pathological_response))
print(pseudobulk)

# Save pseudobulk table
fwrite(as.data.frame(pseudobulk),
       file.path(OUT_TAB, "Bloc4A_pseudobulk_MPR_pCR_scores.csv"))
message("Saved: Bloc4A_pseudobulk_MPR_pCR_scores.csv")

# 5) Wilcoxon test MPR vs pCR, patient-level
# NOTE: n=10 MPR, n=11 pCR, sufficient for Wilcoxon
message("Wilcoxon test MPR vs pCR (patient-level)...")
pvals_patient <- sapply(score_cols, function(s) {
  wilcox.test(pseudobulk[[s]][pseudobulk$pathological_response == "MPR"],
              pseudobulk[[s]][pseudobulk$pathological_response == "pCR"])$p.value
})
print(round(pvals_patient, 4))

# 6) Dotplot, patient-level scores with mean + p-value
# NOTE: each dot = one patient, individual variability visible
message("Generating dotplot...")

format_p <- function(p) {
  ifelse(round(p, 4) == 0, "p < 0.0001", paste0("p = ", round(p, 4)))
}

pval_df <- data.frame(
  signature = score_cols,
  label     = sapply(score_cols, function(s) format_p(pvals_patient[s]))
)

pseudobulk_long <- pseudobulk %>%
  pivot_longer(cols = all_of(score_cols),
               names_to  = "signature",
               values_to = "score")

p_dot <- ggplot(pseudobulk_long,
                aes(x = pathological_response, y = score,
                    color = pathological_response)) +
  geom_boxplot(width = 0.3, alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  scale_color_manual(values = c("MPR" = "#4393C3", "pCR" = "#1A7A1A")) +
  facet_wrap(~signature, scales = "free_y", ncol = 2) +
  geom_text(data = pval_df,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 2, size = 3.5, color = "black") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x    = element_blank(),
        strip.text      = element_text(size = 11, face = "bold")) +
  ylab("Mean UCell score (patient-level)")

ggsave(file.path(OUT_FIG, "Bloc4A_pseudobulk_MPR_pCR_dotplot.png"),
       p_dot, width = 10, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4A_pseudobulk_MPR_pCR_dotplot.png")

message("DONE Bloc4A Script 06")