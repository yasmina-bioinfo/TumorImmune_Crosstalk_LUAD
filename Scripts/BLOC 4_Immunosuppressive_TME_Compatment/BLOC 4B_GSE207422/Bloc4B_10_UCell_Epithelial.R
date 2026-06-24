#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 10: UCell scoring on epithelial cells
# 5 groups: tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, Ciliated
# Signatures: 7 MSigDB Hallmark + 2 custom (HSF1_targets, Antigen_presentation)
# References: Liberzon et al. Cell Systems 2015 (MSigDB Hallmark)
#             Mendillo et al. Cell 2012 (HSF1 targets)
#             Hu et al. Genome Medicine 2023 (Antigen presentation)
# Input:  Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/UCell_Epithelial/
#         Results/Tables/Bloc4B_10_UCell_Epithelial_scores.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs/UCell_Epithelial")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# 1) Load epithelial object with final groups
message("Loading epithelial object...")
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds"))
message("Total cells: ", ncol(seu_Epi))
message("Group distribution:")
print(table(seu_Epi$final_group, useNA = "always"))

# 2) Load MSigDB Hallmark signatures
message("Loading MSigDB Hallmark signatures...")
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")

get_hallmark <- function(name) {
  hallmarks %>% filter(gs_name == name) %>% pull(gene_symbol) %>% unique()
}

# 3) Define all signatures
signatures <- list(
  # MSigDB Hallmark
  Proliferation_E2F     = get_hallmark("HALLMARK_E2F_TARGETS"),
  Proliferation_G2M     = get_hallmark("HALLMARK_G2M_CHECKPOINT"),
  Apoptosis             = get_hallmark("HALLMARK_APOPTOSIS"),
  EMT                   = get_hallmark("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  IFN_gamma_response    = get_hallmark("HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  IL6_JAK_STAT3         = get_hallmark("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  TNFA_NFkB             = get_hallmark("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  WNT_beta_catenin      = get_hallmark("HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
  Notch_signaling       = get_hallmark("HALLMARK_NOTCH_SIGNALING"),
  Unfolded_protein      = get_hallmark("HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  Hypoxia               = get_hallmark("HALLMARK_HYPOXIA"),
  
  # Custom signatures
  HSF1_targets          = c("HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", 
                            "HSPB1", "HSPA5", "HSPH1", "DNAJB1"),
  Antigen_presentation  = c("HLA-DRA", "HLA-DRB1", "CD74", "CIITA")
)

message("Signatures defined: ", length(signatures))

# 4) Run UCell
message("Running UCell...")
seu_Epi <- AddModuleScore_UCell(seu_Epi, features = signatures, name = "_UCell")
message("UCell done.")

# 5) Save scores
score_cols <- paste0(names(signatures), "_UCell")

# Check which columns were created
message("Score columns created:")
print(intersect(score_cols, colnames(seu_Epi@meta.data)))

scores_df <- seu_Epi@meta.data %>%
  dplyr::select(final_group, PathResponse, TME_cell_type, 
                any_of(score_cols))
scores_df$cell <- colnames(seu_Epi)

fwrite(as.data.frame(scores_df), 
       file.path(OUT_TAB, "Bloc4B_10_UCell_Epithelial_scores.csv"))
message("Saved: Bloc4B_10_UCell_Epithelial_scores.csv")

# 6) Dotplot , all signatures x all groups
message("Producing dotplot...")

seu_sub <- subset(seu_Epi, cells = colnames(seu_Epi)[!is.na(seu_Epi$final_group)])
Idents(seu_sub) <- "final_group"

# Get actual score column names
score_cols_present <- intersect(score_cols, colnames(seu_sub@meta.data))

p1 <- DotPlot(seu_sub, 
              features = score_cols_present,
              group.by = "final_group",
              cols = c("lightblue", "red"),
              dot.scale = 8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 18, color = "black", face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)) +
  ggtitle("UCell scores — Epithelial compartment GSE207422") +
  xlab("") + ylab("")

ggsave(file.path(OUT_FIG, "Bloc4B_10_Dotplot_UCell_Epithelial.png"),
       p1, width = 20, height = 8, dpi = 300)
message("Saved: Dotplot")

# 7) Violin plots per signature (supplementary)
message("Producing violin plots...")

group_colors <- c(
  "tumor_MPR"   = "#E63946",
  "normal_MPR"  = "#457B9D",
  "tumor_NMPR"  = "#F4A261",
  "normal_NMPR" = "#2A9D8F",
  "Ciliated"    = "#ADB5BD"
)

for (sig in score_cols_present) {
  tryCatch({
    p <- VlnPlot(seu_sub, features = sig,
                 group.by = "final_group",
                 cols = group_colors,
                 pt.size = 0) +
      ggtitle(gsub("_UCell", "", sig)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    fname <- paste0("Bloc4B_10_Violin_", gsub("_UCell", "", sig), ".png")
    ggsave(file.path(OUT_FIG, fname), p, width = 8, height = 5, dpi = 150)
  }, error = function(e) {
    message("Violin failed for ", sig, ": ", e$message)
  })
}
message("Saved: Violin plots")

# Save updated Seurat object
saveRDS(seu_Epi, file.path(DATA_DIR, "Objects/Bloc4B_10_seu_Epithelial_UCell.rds"))
message("Saved: Bloc4B_10_seu_Epithelial_UCell.rds")

message("DONE UCell Epithelial")