#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 07: Epithelial cells extraction
# Extracts epithelial populations from TME object for downstream analysis
# Includes Ciliated_epithelial as normal reference for CopyKAT
# Input:  C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/objects/04_TME_MPR_NMPR.rds
# Output: Objects/Bloc4B_07_seu_Epithelial.rds
#         Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_07_UMAP_Epithelial.png
#         Results/Tables/Bloc4B_07_Epithelial_cellcounts.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(scales)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy/objects/04_TME_MPR_NMPR.rds"
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load TME object
message("Loading TME object...")
seu_TME <- readRDS(IN_OBJ)
message("Total cells: ", ncol(seu_TME))

# 2) Extract epithelial populations
# NOTE: Ciliated_epithelial included as normal reference for CopyKAT
# Tumor_epithelial subtypes expected to be aneuploid (malignant)
epithelial_types <- c("Tumor_epithelial",
                      "Tumor_epithelial_AT2",
                      "Tumor_epithelial_basal",
                      "Tumor_epithelial_EMT",
                      "Ciliated_epithelial")

seu_Epi <- subset(seu_TME,
                  subset = TME_cell_type %in% epithelial_types)
message("Epithelial cells extracted: ", ncol(seu_Epi))
print(table(seu_Epi$TME_cell_type))
print(table(seu_Epi$PathResponse))

# 3) Cell counts table
cell_counts <- as.data.frame(table(seu_Epi$TME_cell_type, seu_Epi$PathResponse))
colnames(cell_counts) <- c("cell_type", "response", "n_cells")
fwrite(cell_counts, file.path(OUT_TAB, "Bloc4B_07_Epithelial_cellcounts.csv"))
message("Saved: Bloc4B_07_Epithelial_cellcounts.csv")

# 4) UMAP by cell type
png(file.path(OUT_FIG, "Bloc4B_07_UMAP_Epithelial_celltype.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_Epi,
              reduction = "umap",
              group.by  = "TME_cell_type",
              label     = TRUE,
              repel     = TRUE,
              raster    = FALSE) +
        theme_bw() +
        theme(legend.text = element_text(size = 12)) +
        labs(title = "Epithelial cells — GSE207422"))

      dev.off()
message("Saved: Bloc4B_07_UMAP_Epithelial_celltype.png")

# 5) UMAP by response
png(file.path(OUT_FIG, "Bloc4B_07_UMAP_Epithelial_response.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_Epi,
              reduction = "umap",
              group.by  = "PathResponse",
              raster    = FALSE) +
        scale_color_manual(values = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
        theme_bw() +
        theme(legend.text = element_text(size = 12)) +
        labs(title = "Epithelial cells by response — GSE207422"))
dev.off()
message("Saved: Bloc4B_07_UMAP_Epithelial_response.png")

# 6) Barplot proportions by response

df_prop <- as.data.frame(table(seu_Epi$TME_cell_type, seu_Epi$PathResponse)) %>%
  rename(cell_type = Var1, response = Var2, n = Freq) %>%
  group_by(response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

epi_colors <- c(
  "Ciliated_epithelial"    = "#74ADD1",
  "Tumor_epithelial"       = "#D73027",
  "Tumor_epithelial_AT2"   = "#1A7A1A",
  "Tumor_epithelial_basal" = "#FC8D59",
  "Tumor_epithelial_EMT"   = "#762A83"
)

p_bar <- ggplot(df_prop, aes(x = response, y = prop, fill = cell_type)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = epi_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_classic() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(size = 13, face = "bold"),
        legend.title  = element_blank()) +
  ylab("Epithelial subtype proportion") +
  labs(title = "Epithelial composition by response — GSE207422")

ggsave(file.path(OUT_FIG, "Bloc4B_07_Barplot_Epithelial_proportions.png"),
       p_bar, width = 8, height = 6, dpi = 300, bg = "white")
message("Saved: Bloc4B_07_Barplot_Epithelial_proportions.png")

# 7) Save epithelial object
message("Saving epithelial object...")
saveRDS(seu_Epi, file.path(OUT_OBJ, "Bloc4B_07_seu_Epithelial.rds"))
message("Saved: Objects/Bloc4B_07_seu_Epithelial.rds")
message("DONE Bloc4B Script 07")

