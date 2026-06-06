# Canonical Cell Type Markers — TME Reference

## CD8 T cells
### General
- CD8A, CD8B, CD3D, CD3E, CD3G

### Exhausted / Cytotoxic TRM (Cluster 2)
- Exhaustion: PDCD1, LAG3, HAVCR2, CXCL13
- Cytotoxic: GZMB, GZMA, GZMH, PRF1, IFNG, NKG7
- Tissue-resident (TRM): ITGA1, ITGAE, CXCR6, ZNF683, HOPX

### Effector Memory (Cluster 3)
- GZMK, EOMES, KLRG1, CCL4, CCL4L2
- Note: GZMK+ distinct from GZMB+ exhausted
- Note: absence of PDCD1/LAG3/HAVCR2 distinguishes from exhausted

### CD4 T cells
- CD4, CD3D, CD3E, IL7R, TCF7

### Naive / Memory activated CD4 T cells (Cluster 1)
- TCF7, LEF1 — naive/stem-like signature
- SATB1 — transcriptional regulator of T cells
- CD40LG — CD4 T helper marker
- KLRB1 — NKT or activated CD4 T cells

### Naive / Resting CD4 T cells (Cluster 10)
- IL7R, CD28, ITK, CAMK4, RORA, STAT4
- Note: absence of activation/exhaustion markers
- Note: distinct from Cluster 1 by absence of TCF7/LEF1/CD40LG

## CD4 T follicular helper (Tfh)
- CXCL13, BTLA, MAF, ICOS, CD40LG, TCF7
- Note: CXCL13+TOX+PDCD1 = exhausted Tfh-like in tumors

## Regulatory T cells (Tregs)
- FOXP3, IL2RA (CD25), CTLA4, IKZF2, TIGIT

### Tumor-infiltrating activated Tregs (Cluster 4)
- CCR8, RTKN2, LAYN, TNFRSF4, TNFRSF18, TNFRSF9 : highly suppressive tumor Tregs
- Note: CCR8+ marks highly suppressive tumor Tregs

## NK cells
- GNLY, NKG7, KLRD1, KLRF1, NCR1

### Cytotoxic NK cells
- GNLY, NKG7, KLRD1, KLRF1, PRF1, GZMB
- S1PR5, CX3CR1, FGFBP2, SPON2 : circulating cytotoxic NK
- FCGR3A (CD16) : cytotoxic NK marker
- TBX21 : NK/effector transcription factor
- Note: HOPX+ = tissue-resident NK subset

## NKT / γδ T cells
- TRDC : gamma-delta T cell marker (TCR delta chain)
- XCL1, XCL2, KIR2DL4, KLRC1, KLRC3
- Note: coexpression TRDC + TRBC1 = mixed NKT/γδ population
- Note: to confirm with SingleR and Azimuth

## B cells
- CD19, MS4A1 (CD20), CD79A, CD79B, IGHM
- CXCR5 : follicular B cells marker

### Naive / Transitional B cells (Cluster 12)
- IGHD, TCL1A, CR2, FCER2, BACH2 : naive B cell markers
- BTK, PAX5 : B cell signaling and lineage
- Note: distinct from Cluster 0 (memory/activated) by IGHD+ TCL1A+

### Memory / Activated B cells (Cluster 0)
- TNFRSF13B, CD40, CXCR5 : memory/activated markers

### Memory B cells (Cluster 16)
- TNFRSF13C, FCRL1, BACH2, PLCG2 : memory B cell markers
- Note: absence of IGHD (not naive) and secreted Ig (not plasma cell)
- Note: distinct from Cluster 0; to confirm with SingleR/Azimuth

## Macrophages / Monocytes
- CD68, CD14, LYZ, CSF1R, MRC1 (CD206)

## Dendritic cells
- ITGAX (CD11c), CLEC9A, FCER1A, CD1C, LILRA4

## Malignant / Epithelial cells
- EPCAM, KRT7, KRT19, KRT18, MUC1

## Fibroblasts
- COL1A1, COL1A2, FAP, ACTA2, DCN

## Endothelial cells
- PECAM1 (CD31), VWF, CDH5, CLDN5, RAMP2

## Mast cells
- CPA3, TPSB2, TPSAB1, TPSD1 : tryptases (highly specific)
- KIT (CD117), GATA2 : mast cell identity
- HDC, MS4A2, HPGDS : histamine/IgE/prostaglandin pathway
- IL1RL1 (ST2) : activated mast cells

## Monocytes / Macrophages
### Inflammatory monocytes
- S100A8, S100A9, S100A12, VCAN, FCN1, TREM1
- CXCL8, IL1B, OSM : pro-inflammatory cytokines
- Note: S100A8/A9 high = classical inflammatory monocyte signature

### Tumor-Associated Macrophages (TAMs) — M2-like
- CD68, CD14, MRC1, CSF1R, CD163
- TREM2, APOE, APOC1 : immunosuppressive/lipid-associated TAMs
- SPP1, MARCO, VSIG4 : pro-tumoral markers
- Note: TREM2+ TAMs associated with immunotherapy resistance in NSCLC

### Classical Monocytes CD14+ (Cluster 19)
- VCAN, FCN1, LYZ, CD36, ANPEP (CD13) : classical monocyte markers
- LILRA1/2/5, LILRB2 : LILR receptors
- Note: distinct from Cluster 6 (inflammatory) by absence of TREM1/CXCL8/IL1B

## Proliferating cells (Cycling)
- MKI67, TOP2A, BIRC5, CDK1, CCNA2, CCNB2
- FOXM1, E2F1 — proliferation transcription factors
- TYMS, RRM2, MCM2, MCM4 — DNA replication (S phase)
- Note: transversal cluster — multiple cell types in active proliferation
- Note: not a cell type but a proliferative state

## Dendritic cells (DCs)
- FCER1A, CD1C, CST3, CSF1R, CD86, CD80 : general DC markers
- ITGAX (CD11c), CLEC9A, LILRA4 : previously listed

### Plasmacytoid Dendritic cells (pDCs) (Cluster 18)
- LILRA4, CLEC4C, IL3RA (CD123), LAMP5, TCL1A
- DNASE1L3, PLD4 : pDC function markers
- Note: major producers of type I interferon
- Note: rare but highly distinctive population

### Mature / Migratory DCs (mregDC) (Cluster 14)
- LAMP3, FSCN1 — migratory DC signature
- CD207, CD1A, CD1E : Langerhans/cDC1-like
- CLEC10A, CLEC4A : cDC2 markers
- AXL : transitional DC
- Note: LAMP3+ FSCN1+ = mregDC, immunosuppressive in tumors

## Plasma cells
- IGHG1-4, IGHA1-2, IGKC, IGLC2-3 : secreted immunoglobulins
- JCHAIN, SDC1 (CD138), CD38 : plasma cell identity
- XBP1, MZB1 : plasmacytic differentiation
- TNFRSF17 (BCMA) : plasma cell survival receptor
- Note: distinct from B cells by absence of CD19/MS4A1 and high Ig secretion

## IFN-stimulated cells (ISG-high) (Cluster 17)
- IFIT1/2/3, ISG15, MX1/2, OAS1/2/3 : interferon-stimulated genes
- STAT1, STAT2, IRF7 : IFN signaling transcription factors
- Note: transversal cluster : multiple cell types in IFN-high state
- Note: STAT2 high : connects to portfolio narrative (STAT2-high non-MPR program)