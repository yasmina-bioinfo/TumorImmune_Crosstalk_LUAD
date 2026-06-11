#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 11: scRepertoire TCR analysis
# TCR clonotype diversity and expansion analysis on CD8 T cells
# NOTE: Run in WSL (R 4.6.0) — scRepertoire requires gsl >= R 4.5.0
#       not available on Windows R 4.4
# Input:  Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
#         Data/GSE243013_T_with_TCR_annotation.csv.gz
# Output: Results/Figures/CD8/Bloc3_scRepertoire_*.png
#         Results/Tables/Bloc3_scRepertoire_summary.csv
# ============================================================

library(Seurat)
library(scRepertoire)
library(data.table)
library(ggplot2)
library(dplyr)

# Paths
DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
TCR_FILE <- file.path(DATA_DIR, "Data/GSE243013_T_with_TCR_annotation.csv.gz")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load CD8 ProjecTILs object
message("Loading CD8 ProjecTILs object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))

# 2) Load TCR data
message("Loading TCR annotation file...")
tcr <- fread(TCR_FILE)
message("TCR file: ", nrow(tcr), " cells")

# 3) Filter TCR to LUAD CD8 cells only
# Keep only cells present in our CD8 Seurat object
cd8_barcodes <- colnames(seu_CD8)
tcr_cd8 <- tcr[cellID %in% cd8_barcodes]
message("TCR matched to CD8 object: ", nrow(tcr_cd8), " cells")

# 4) Add TCR metadata to Seurat object
# NOTE: same barcode mapping approach as previous scripts
tcr_map <- setNames(tcr_cd8$clonotype, tcr_cd8$cellID)
test <- tcr_map[colnames(seu_CD8)]
names(test) <- colnames(seu_CD8)
seu_CD8$clonotype <- test

expansion_map <- setNames(tcr_cd8$expansion, tcr_cd8$cellID)
test2 <- expansion_map[colnames(seu_CD8)]
names(test2) <- colnames(seu_CD8)
seu_CD8$expansion <- test2

clono_num_map <- setNames(tcr_cd8$clonotype_number, tcr_cd8$cellID)
test3 <- clono_num_map[colnames(seu_CD8)]
names(test3) <- colnames(seu_CD8)
seu_CD8$clonotype_number <- test3

message("TCR coverage in CD8 object:")
print(table(!is.na(seu_CD8$clonotype)))

# 5) Clonotype expansion by ProjecTILs state
message("Clonotype expansion by CD8 state...")
expansion_summary <- seu_CD8@meta.data %>%
  filter(!is.na(functional.cluster), !is.na(expansion)) %>%
  group_by(functional.cluster, expansion) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(functional.cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

print(expansion_summary)

# 6) Barplot expansion by CD8 state
p_exp_state <- ggplot(expansion_summary,
                      aes(x = functional.cluster, y = prop, fill = expansion)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = c("expanded" = "#D73027", "non-expanded" = "#4393C3")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Proportion") +
  theme_classic() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y   = element_text(size = 11),
        legend.title  = element_blank(),
        legend.text   = element_text(size = 10)) +
  labs(title = "Clonotype expansion by CD8 state — GSE243013 LUAD")

ggsave(file.path(OUT_FIG, "Bloc3_scRepertoire_expansion_by_state.png"),
       p_exp_state, width = 8, height = 5, dpi = 300, bg = "white")
message("Saved: Bloc3_scRepertoire_expansion_by_state.png")

# 7) Clonotype expansion by pathological response
expansion_response <- seu_CD8@meta.data %>%
  filter(!is.na(pathological_response), !is.na(expansion)) %>%
  group_by(pathological_response, expansion) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pathological_response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_exp_response <- ggplot(expansion_response,
                         aes(x = pathological_response, y = prop, fill = expansion)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = c("expanded" = "#D73027", "non-expanded" = "#4393C3")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Proportion") +
  theme_classic() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(size = 13, face = "bold"),
        axis.text.y   = element_text(size = 11),
        legend.title  = element_blank(),
        legend.text   = element_text(size = 10)) +
  labs(title = "Clonotype expansion by pathological response — GSE243013 LUAD")

ggsave(file.path(OUT_FIG, "Bloc3_scRepertoire_expansion_by_response.png"),
       p_exp_response, width = 7, height = 5, dpi = 300, bg = "white")
message("Saved: Bloc3_scRepertoire_expansion_by_response.png")

# 8) Clonal diversity by response group
# NOTE: diversity computed per patient to avoid N bias
diversity_summary <- seu_CD8@meta.data %>%
  filter(!is.na(sampleID), !is.na(clonotype)) %>%
  group_by(sampleID, pathological_response) %>%
  summarise(
    n_cells     = n(),
    n_clonotypes = n_distinct(clonotype),
    diversity    = n_distinct(clonotype) / n(),
    .groups = "drop"
  )

message("Diversity summary by response:")
print(diversity_summary %>%
        group_by(pathological_response) %>%
        summarise(mean_diversity = mean(diversity),
                  sd_diversity   = sd(diversity),
                  .groups = "drop"))

# 9) Save summary tables
fwrite(expansion_summary, file.path(OUT_TAB, "Bloc3_scRepertoire_expansion_state.csv"))
fwrite(expansion_response, file.path(OUT_TAB, "Bloc3_scRepertoire_expansion_response.csv"))
fwrite(diversity_summary, file.path(OUT_TAB, "Bloc3_scRepertoire_diversity.csv"))
message("Saved: scRepertoire summary tables")

# 10) Save updated Seurat object with TCR metadata
message("Saving updated CD8 object with TCR metadata...")
saveRDS(seu_CD8, file.path(DATA_DIR, "Objects/Bloc3_11_seu_CD8_TCR.rds"))
message("Saved: Objects/Bloc3_11_seu_CD8_TCR.rds")
message("DONE Bloc3 Script 11")