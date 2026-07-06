#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 13: Pseudobulk Hallmark UCell scores MPR vs pCR — CD8 T cells
# Patient-level comparison to test MPR as intermediate state hypothesis
# NOTE: non-MPR excluded, focus on MPR (n=10) vs pCR (n=11)
# Uses the same 21 Hallmark signatures established for the MPR vs non-MPR
# analysis (sig_filtered), restricted to CD8.TEX/CD8.TPEX for consistency
# with the main narrative's focus on these two states.
# Pseudobulk: mean UCell score per patient per signature
# Wilcoxon rank-sum test MPR vs pCR (n=10 vs n=11), BH-corrected
# References: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)
# Input:  Objects/Bloc3_GSE243013_12c_seu_CD8_Hallmark.rds
# Output: Results/Figures/CD8/Bloc3_pseudobulk_CD8_MPR_pCR_dotplot.png
#         Results/Tables/Bloc3_pseudobulk_CD8_MPR_pCR_Hallmark_scores.csv
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
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_GSE243013_12c_seu_CD8_Hallmark.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load Hallmark CD8 object
message("Loading Hallmark CD8 object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))

# 2) Define score columns — same 21 signatures used for MPR vs non-MPR
sig_filtered <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_COMPLEMENT"
)

# 3) Subset MPR and pCR, CD8.TEX/CD8.TPEX only
# NOTE: restricted to TEX/TPEX to match the scope of the main CD8 narrative
message("Subsetting MPR and pCR, CD8.TEX/CD8.TPEX only...")
seu_sub <- subset(seu_CD8,
                  subset = pathological_response %in% c("MPR", "pCR") &
                    functional.cluster %in% c("CD8.TEX", "CD8.TPEX"))
message("Cells MPR + pCR: ", ncol(seu_sub))
print(table(seu_sub$pathological_response))

# 4) Pseudobulk: mean UCell score per patient per signature
# NOTE: pseudobulk aggregates cell-level scores to patient-level means
# This avoids pseudo-replication and enables proper patient-level statistics
message("Computing pseudobulk scores...")
pseudobulk <- seu_sub@meta.data %>%
  filter(!is.na(sampleID), !is.na(pathological_response)) %>%
  select(sampleID, pathological_response, all_of(sig_filtered)) %>%
  group_by(sampleID, pathological_response) %>%
  summarise(across(all_of(sig_filtered), mean, .names = "{.col}"),
            n_cells = n(),
            .groups = "drop")

message("Patients per group:")
print(table(pseudobulk$pathological_response))
print(pseudobulk)

# Save pseudobulk table
fwrite(as.data.frame(pseudobulk),
       file.path(OUT_TAB, "Bloc3_pseudobulk_CD8_MPR_pCR_Hallmark_scores.csv"))
message("Saved: Bloc3_pseudobulk_CD8_MPR_pCR_Hallmark_scores.csv")

# 5) Wilcoxon test MPR vs pCR, patient-level, BH-corrected
message("Wilcoxon test MPR vs pCR (patient-level)...")

pval_summary <- lapply(sig_filtered, function(s) {
  mpr <- pseudobulk[[s]][pseudobulk$pathological_response == "MPR"]
  pcr <- pseudobulk[[s]][pseudobulk$pathological_response == "pCR"]
  test <- wilcox.test(mpr, pcr)
  data.frame(signature   = s,
             p_value     = test$p.value,
             median_MPR  = median(mpr),
             median_pCR  = median(pcr))
}) %>% bind_rows()

pval_summary <- pval_summary %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

fwrite(pval_summary,
       file.path(OUT_TAB, "Bloc3_pseudobulk_CD8_MPR_pCR_Hallmark_results.csv"))
message("Saved: Bloc3_pseudobulk_CD8_MPR_pCR_Hallmark_results.csv")
print(pval_summary)

# 6) Dotplot, patient-level scores with mean + adjusted p-value
# NOTE: each dot = one patient, individual variability visible
message("Generating dotplot...")

format_p <- function(p) {
  ifelse(round(p, 4) == 0, "p_adj < 0.0001", paste0("p_adj = ", round(p, 4)))
}

pval_df <- pval_summary %>%
  mutate(label = format_p(p_adj),
         signature = gsub("HALLMARK_", "", signature)) %>%
  select(signature, label)

pseudobulk_long <- pseudobulk %>%
  pivot_longer(cols = all_of(sig_filtered),
               names_to  = "signature",
               values_to = "score") %>%
  mutate(signature = gsub("HALLMARK_", "", signature))

pval_df <- pval_df %>%
  mutate(signature = gsub("HALLMARK_", "", signature))

p_dot <- ggplot(pseudobulk_long,
                aes(x = pathological_response, y = score,
                    color = pathological_response)) +
  geom_boxplot(width = 0.3, alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = c("MPR" = "#4393C3", "pCR" = "#1A7A1A")) +
  facet_wrap(~signature, scales = "free_y", ncol = 4) +
  geom_text(data = pval_df,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 2, size = 2.8, color = "black") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x    = element_blank(),
        strip.text      = element_text(size = 9, face = "bold")) +
  ylab("Mean UCell Hallmark score (patient-level)")

ggsave(file.path(OUT_FIG, "Bloc3_pseudobulk_CD8_MPR_pCR_dotplot.png"),
       p_dot, width = 16, height = 20, dpi = 300, bg = "white")
message("Saved: Bloc3_pseudobulk_CD8_MPR_pCR_dotplot.png")

message("DONE Bloc3 Script 13")