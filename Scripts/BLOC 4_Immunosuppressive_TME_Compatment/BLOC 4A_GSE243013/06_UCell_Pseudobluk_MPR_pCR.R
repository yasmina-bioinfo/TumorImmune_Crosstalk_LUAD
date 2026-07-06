#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 06: Pseudobulk Hallmark UCell scores MPR vs pCR — TAMs
# Patient-level comparison to test MPR as intermediate state hypothesis
# NOTE: recomputes Hallmark UCell scores from the original (pre-filter) object,
# since GSE243013_seu_TAMs_Hallmark.rds excludes pCR at creation (Script 05c)
# Pseudobulk: mean UCell score per patient per signature
# Wilcoxon rank-sum test MPR vs pCR (n=10 vs n=11), BH-corrected
# References: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)
# Input:  Objects/Bloc4A_05_seu_TAMs_UCell.rds
# Output: Results/Figures/BLOC4A_TAMs/Bloc4A_pseudobulk_MPR_pCR_dotplot.png
#         Results/Tables/BLOC4A/Bloc4A_pseudobulk_MPR_pCR_Hallmark_results.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggpubr)
  library(BiocParallel)
})
register(SnowParam(workers = 1))

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4A_05_seu_TAMs_UCell.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC4A")

# 1) Load original TAMs object (includes MPR, non-MPR, pCR)
message("Loading original TAMs object (pre-filter, includes pCR)...")
seu_TAM <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_TAM))
print(table(seu_TAM$pathological_response))
DefaultAssay(seu_TAM) <- "RNA"

# 2) Define the 25 Hallmark signatures used for MPR vs non-MPR
sig_filtered <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_NOTCH_SIGNALING"
)

# 3) Load Hallmark gene sets and restrict to the 25 already-established signatures
message("Loading Hallmark signatures...")
hallmark_sets_all <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
hallmark_sets <- hallmark_sets_all[sig_filtered]

# 4) Subset MPR and pCR only (before scoring, to save compute time)
message("Subsetting MPR and pCR...")
seu_sub <- subset(seu_TAM,
                  subset = pathological_response %in% c("MPR", "pCR"))
message("Cells MPR + pCR: ", ncol(seu_sub))
print(table(seu_sub$pathological_response))

# 5) Compute UCell Hallmark scores (25 signatures) on this subset
message("Computing UCell Hallmark scores (25 signatures)...")
seu_sub <- AddModuleScore_UCell(seu_sub,
                                features = hallmark_sets,
                                name     = "")
message("Hallmark UCell scores added.")

# 6) Pseudobulk: mean UCell score per patient per signature
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

fwrite(as.data.frame(pseudobulk),
       file.path(OUT_TAB, "Bloc4A_pseudobulk_MPR_pCR_Hallmark_scores.csv"))
message("Saved: Bloc4A_pseudobulk_MPR_pCR_Hallmark_scores.csv")

# 7) Wilcoxon test MPR vs pCR, patient-level, BH-corrected
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
       file.path(OUT_TAB, "Bloc4A_pseudobulk_MPR_pCR_Hallmark_results.csv"))
message("Saved: Bloc4A_pseudobulk_MPR_pCR_Hallmark_results.csv")
print(pval_summary)

# 8) Dotplot, patient-level scores with mean + adjusted p-value
message("Generating dotplot...")

format_p <- function(p) {
  ifelse(round(p, 4) == 0, "p_adj < 0.0001", paste0("p_adj = ", round(p, 4)))
}

pval_df <- pval_summary %>%
  mutate(label     = format_p(p_adj),
         signature = gsub("HALLMARK_", "", signature)) %>%
  select(signature, label)

pseudobulk_long <- pseudobulk %>%
  pivot_longer(cols = all_of(sig_filtered),
               names_to  = "signature",
               values_to = "score") %>%
  mutate(signature = gsub("HALLMARK_", "", signature))

p_dot <- ggplot(pseudobulk_long,
                aes(x = pathological_response, y = score,
                    color = pathological_response)) +
  geom_boxplot(width = 0.3, alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  scale_color_manual(values = c("MPR" = "#4393C3", "pCR" = "#1A7A1A")) +
  facet_wrap(~signature, scales = "free_y", ncol = 5) +
  geom_text(data = pval_df,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 2, size = 2.6, color = "black") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x    = element_blank(),
        strip.text      = element_text(size = 8, face = "bold")) +
  ylab("Mean UCell Hallmark score (patient-level)")

ggsave(file.path(OUT_FIG, "Bloc4A_pseudobulk_MPR_pCR_dotplot.png"),
       p_dot, width = 20, height = 20, dpi = 300, bg = "white")
message("Saved: Bloc4A_pseudobulk_MPR_pCR_dotplot.png")

message("DONE Bloc4A Script 06")