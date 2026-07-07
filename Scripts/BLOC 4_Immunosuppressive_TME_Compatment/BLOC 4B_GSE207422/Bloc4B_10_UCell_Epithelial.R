#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 10: UCell Hallmark scoring on epithelial cells
# Unbiased approach, all 50 Hallmark signatures tested (custom signatures removed)
# 5 groups: tumor_MPR, tumor_NMPR, normal_MPR, normal_NMPR, Ciliated
# Reference: Liberzon et al., Cell Systems 2015 (MSigDB Hallmark)
# Input:  Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds
# Output: Results/Tables/Bloc4B_10_UCell_Epithelial_scores.csv
#         Results/Tables/Bloc4B_10_UCell_Epithelial_wilcox_results.csv
#         Results/Figures/BLOC4B_Epithelial_TAMs/UCell_Epithelial/
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(BiocParallel)
})
register(SnowParam(workers = 1))

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/UCell_Epithelial")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load epithelial object with final groups
message("Loading epithelial object...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds"))
message("Total cells: ", ncol(seu_Epi))
print(table(seu_Epi$final_group, useNA = "always"))

# 2) Load ALL 50 Hallmark signatures (unbiased)
message("Loading all 50 Hallmark signatures...")
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
hallmark_cols <- names(hallmark_sets)
message("Total Hallmark signatures: ", length(hallmark_sets))

# 3) Compute UCell scores (all 50)
message("Computing UCell Hallmark scores...")
seu_Epi <- AddModuleScore_UCell(seu_Epi,
                                features = hallmark_sets,
                                name     = "")
message("Hallmark UCell scores added.")

# 4) Save full per-cell scores table
scores_df <- seu_Epi@meta.data %>%
  dplyr::select(final_group, PathResponse, TME_cell_type,
                all_of(hallmark_cols))
scores_df$cell <- colnames(seu_Epi)

fwrite(as.data.frame(scores_df),
       file.path(OUT_TAB, "Bloc4B_10_UCell_Epithelial_scores.csv"))
message("Saved: Bloc4B_10_UCell_Epithelial_scores.csv")

# 5) Wilcoxon tests — 3 comparisons, BH-corrected within each comparison
message("Running Wilcoxon tests for the 3 comparisons...")

run_comparison <- function(meta, filter_expr, group_col, level_a, level_b, label) {
  sub_meta <- meta %>% filter(!!rlang::parse_expr(filter_expr))
  results <- lapply(hallmark_cols, function(sig) {
    a <- sub_meta[sub_meta[[group_col]] == level_a, sig]
    b <- sub_meta[sub_meta[[group_col]] == level_b, sig]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) < 3 | length(b) < 3) return(NULL)
    test <- wilcox.test(a, b)
    data.frame(comparison = label, signature = sig, p_value = test$p.value,
               median_A = median(a), median_B = median(b))
  }) %>% bind_rows()
  names(results)[names(results) == "median_A"] <- paste0("median_", level_a)
  names(results)[names(results) == "median_B"] <- paste0("median_", level_b)
  results
}

meta <- seu_Epi@meta.data

res_tumor    <- run_comparison(meta, "final_group %in% c('tumor_MPR','tumor_NMPR')",
                               "final_group", "tumor_MPR", "tumor_NMPR", "tumor_MPR_vs_NMPR")
res_normal   <- run_comparison(meta, "final_group %in% c('normal_MPR','normal_NMPR')",
                               "final_group", "normal_MPR", "normal_NMPR", "normal_MPR_vs_NMPR")
res_ciliated <- run_comparison(meta, "final_group == 'Ciliated'",
                               "PathResponse", "MPR", "NMPR", "Ciliated_MPR_vs_NMPR")

wilcox_results <- bind_rows(res_tumor, res_normal, res_ciliated) %>%
  group_by(comparison) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(comparison, p_adj)

fwrite(wilcox_results,
       file.path(OUT_TAB, "Bloc4B_10_UCell_Epithelial_wilcox_results.csv"))
message("Saved: Bloc4B_10_UCell_Epithelial_wilcox_results.csv")

for (comp in unique(wilcox_results$comparison)) {
  n_sig <- sum(wilcox_results$comparison == comp & wilcox_results$p_adj < 0.05, na.rm = TRUE)
  message(comp, ": ", n_sig, "/", length(hallmark_cols), " signatures significant (p_adj < 0.05)")
}

# 6) Discriminant signatures (post-hoc biological relevance filter)
irrelevant_sets <- c(
  "HALLMARK_ADIPOGENESIS", "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_MYOGENESIS",
  "HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_SPERMATOGENESIS",
  "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_COAGULATION", "HALLMARK_PEROXISOME"
)

wilcox_discriminant <- wilcox_results %>%
  filter(!signature %in% irrelevant_sets, p_adj < 0.05)

message("Discriminant signatures retained (biologically relevant, p_adj<0.05): ",
        length(unique(wilcox_discriminant$signature)))

fwrite(wilcox_discriminant,
       file.path(OUT_TAB, "Bloc4B_10_UCell_Epithelial_wilcox_discriminant.csv"))
message("Saved: Bloc4B_10_UCell_Epithelial_wilcox_discriminant.csv")

# 7) Dotplot — descriptive visualization, all 50 signatures x all groups
message("Producing dotplot...")

seu_sub <- subset(seu_Epi, cells = colnames(seu_Epi)[!is.na(seu_Epi$final_group)])
Idents(seu_sub) <- "final_group"

p1 <- DotPlot(seu_sub,
              features = hallmark_cols,
              group.by = "final_group",
              cols = c("lightblue", "red"),
              dot.scale = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 18, color = "black", face = "bold")) +
  ggtitle("UCell Hallmark scores — Epithelial compartment GSE207422 (all 50 signatures)") +
  xlab("") + ylab("")

ggsave(file.path(OUT_FIG, "Bloc4B_10_Dotplot_UCell_Epithelial_all50.png"),
       p1, width = 30, height = 8, dpi = 300)
message("Saved: Dotplot (all 50, descriptive)")

# 8) Violin plots- discriminant signatures only (supplementary)
message("Producing violin plots for discriminant signatures...")

group_colors <- c(
  "tumor_MPR"   = "#E63946",
  "normal_MPR"  = "#457B9D",
  "tumor_NMPR"  = "#F4A261",
  "normal_NMPR" = "#2A9D8F",
  "Ciliated"    = "#ADB5BD"
)

for (sig in unique(wilcox_discriminant$signature)) {
  tryCatch({
    p <- VlnPlot(seu_sub, features = sig,
                 group.by = "final_group",
                 cols = group_colors,
                 pt.size = 0) +
      ggtitle(gsub("HALLMARK_", "", sig)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    fname <- paste0("Bloc4B_10_Violin_", gsub("HALLMARK_", "", sig), ".png")
    ggsave(file.path(OUT_FIG, fname), p, width = 8, height = 5, dpi = 150)
  }, error = function(e) {
    message("Violin failed for ", sig, ": ", e$message)
  })
}
message("Saved: Violin plots (discriminant signatures)")

# 9) Save updated Seurat object
saveRDS(seu_Epi, file.path(DATA_DIR, "Objects/Bloc4B_10_seu_Epithelial_UCell.rds"))
message("Saved: Bloc4B_10_seu_Epithelial_UCell.rds")

message("DONE UCell Epithelial")