#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc2 Script 04: Azimuth reference-based annotation
# Environment: Ubuntu 24.04 LTS (WSL2) + R 4.6.0 + VS Code
# NOTE: This script must be run in WSL ; NOT in RStudio/Windows
#
# INSTALLATION CHALLENGES (documented for full reproducibility) :
#
# 1. Azimuth on Windows (RStudio) — FAILED
#    Dependencies BSgenome.Hsapiens.UCSC.hg38 and EnsDb.Hsapiens.v86
#    not installable on Windows. Documented in Script 03 (sctype).
#
# 2. WSL setup : Ubuntu 24.04 + R upgrade
#    Default Ubuntu R = 4.3.3. TFMPvalue requires R >= 4.5.0.
#    Solution: added noble-cran40 CRAN repository → upgraded to R 4.6.0
#    Commands used:
#      wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc |
#        sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
#      sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/"
#      sudo apt-get update && sudo apt-get install -y r-base r-base-dev
#
# 3. System libraries required (Ubuntu):
#    sudo apt-get install -y libgsl-dev libgsl27 libfontconfig1-dev
#      libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev
#      libtiff5-dev libjpeg-dev libnode-dev build-essential
#      libcurl4-openssl-dev libssl-dev libxml2-dev
#
# 4. R dependencies installed in WSL:
#    BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86",
#      "glmGamPoi", "JASPAR2020", "TFBSTools"))
#    install.packages(c("DT", "ggplot2", "googlesheets4", "patchwork", "shiny",
#      "shinyBS", "shinydashboard", "shinyjs", "plotly", "remotes"))
#    BiocManager::install(c("Seurat", "SeuratObject", "Signac"))
#    remotes::install_github("mojaveazure/seurat-disk")
#    remotes::install_github("satijalab/seurat-data")
#    remotes::install_github("bnprks/BPCells/r")
#    remotes::install_github("satijalab/azimuth", ref = "master")
#    Total installation time: ~3 hours
#
# 5. lungref.SeuratData (833 Mb) , download challenges
#    InstallData("lungref") failed twice due to connection timeout.
#    Solution: wget with resume option (-c):
#      wget -c http://seurat.nygenome.org/src/contrib/lungref.SeuratData_2.0.0.tar.gz
#    Then installed locally:
#      install.packages("lungref.SeuratData_2.0.0.tar.gz", repos=NULL, type="source")
#    Total download time: ~2 hours with 3 connection interruptions (wget resumed each time)
#
# 6. SeuratObject version conflict
#    lungref reference built with Seurat 4.x ; Key<- function removed in v5
#    Solution: downgrade SeuratObject to 4.1.3
#      remotes::install_version("SeuratObject", version="4.1.3",
#        repos="https://cran.r-project.org/")
#
# 7. BPCells path conflict : BLOCKING on 16GB RAM laptop
#    Seurat object stores Windows absolute path to BPCells matrix.
#    In WSL, this path does not resolve.
#    Workaround attempted: open_matrix_dir() with WSL path + as(mat, "dgCMatrix")
#    Result: process killed : insufficient RAM to convert 31831 x 298567 matrix
#    CONCLUSION: Azimuth with this dataset requires a compute server (>32GB RAM)
#    This script is retained for documentation and future server-based execution.
#
# WHAT WAS SUCCESSFULLY VALIDATED IN WSL (despite RAM constraints) :
#    - Azimuth installed and functional in WSL
#    - lungref reference (584,884 cells, human lung) loaded successfully
#    - reference <- LoadData("lungref", type="azimuth") executed without error
#    - RunAzimuth() blocked by RAM constraint on local machine
#
# Reference: Hao et al., bioRxiv 2022
#            https://www.biorxiv.org/content/10.1101/2022.02.24.481684v1
# ============================================================

library(Azimuth)
library(SeuratData)
library(Seurat)
library(BPCells)

#  TO RUN ON SERVER :
# Uncomment and execute the following when running on a server with >32GB RAM

# # 1) Load Seurat object
# seu <- readRDS("/path/to/Objects/Bloc2_03_seu_sctype.rds")
#
# # 2) Rebuild BPCells matrix with server path
# # NOTE: BPCells stores absolute Windows path in the object
# # On server, must reopen matrix from correct path and reassign
# mat <- open_matrix_dir("/path/to/Objects/01_LUAD_matrix_BPCells")
# seu[["RNA"]] <- CreateAssayObject(counts = as(mat, "dgCMatrix"))
#
# # 3) Run Azimuth with lung reference
# seu <- RunAzimuth(seu, reference = "lungref")
#
# # 4) Check annotations
# table(seu$predicted.annotation.l1)
# table(seu$predicted.annotation.l2)
#
# # 5) UMAP colored by Azimuth annotation
# DimPlot(seu, group.by = "predicted.annotation.l1", label = TRUE) + theme_bw()
#
# # 6) Save
# saveRDS(seu, "Objects/Bloc2_04_seu_azimuth.rds")

message("Azimuth script : server execution required. See header for full documentation.")
message("Validated: Azimuth installed, lungref loaded. RunAzimuth blocked by RAM on local machine.")