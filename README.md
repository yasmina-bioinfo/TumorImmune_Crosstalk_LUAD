# TumorImmune_Crosstalk_LUAD

> **How does tumor cell plasticity shape CD8 T cell functional states and their response to anti-PD-1 in lung adenocarcinoma?**

This project investigates the crosstalk between malignant epithelial cells and CD8 T cells in the tumor microenvironment (TME) of lung adenocarcinoma (LUAD), using neoadjuvant chemo-immunotherapy data with matched pathological response (MPR vs. non-MPR).

This work is an extension to validate findings from [CD8_NSCLC_scRNAseq](https://github.com/yasmina-bioinfo/CD8_NSCLC_scRNAseq) in an independent, larger cohort with paired TCR sequencing.

---

## Dataset

| Field | Details |
|---|---|
| **Accession** | [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013) |
| **Publication** | Liu et al., *Cell* 2025 |
| **Histology** | LUAD only (n = 63, filtered from 234 NSCLC) |
| **Modalities** | scRNA-seq + scTCR-seq |
| **Treatment** | Neoadjuvant chemo + anti-PD-1 |
| **Response variable** | Pathological response: MPR vs. non-MPR |
| **Timing** | Post-treatment surgical resection |

---

## Biological question

Building on a CD8_Exhausted_Terminal enrichment signal in MPR patients (OR = 3.36, GSE207422), this project asks:

1. What transcriptional programs define malignant epithelial heterogeneity in LUAD under treatment pressure?
2. Do specific tumor cell states associate with impaired CD8 T cell function in non-MPR patients?
3. What ligand-receptor axes mediate epithelial–CD8 crosstalk differentially in MPR vs. non-MPR?
4. Which transcription factors drive these programs in each compartment?

---

## Analytical pipeline

### Block 0 — Data acquisition
- Download GSE243013 (expression matrix + TCR metadata)
- Filter to LUAD samples (n = 63)
- Harmonize clinical metadata (MPR / non-MPR / pCR classification)

### Block 1 — QC and preprocessing
- Quality control, normalization, HVG selection (`Seurat`)
- Batch correction across patients (`Harmony`)
- Global TME clustering

### Block 2 — Global TME annotation
- Cluster marker identification (top markers per cluster)
- Manual annotation + automated validation (`SingleR`)

### Block 3 — CD8 T cell focus
- Robust state annotation on reference atlas (`ProjecTILs`)
- Functional module scoring: exhaustion, cytotoxicity, memory (`UCell`)
- TCR integration: clonotype expansion, repertoire diversity (`scRepertoire`)
- MPR vs. non-MPR comparison of CD8 state composition

### Block 4 — Malignant epithelial focus
- Isolation and sub-clustering of epithelial population
- Malignancy confirmation (`CopyKAT`)
- Functional scoring: EMT, IFN stress, proliferation (`UCell`)
- MPR vs. non-MPR comparison of tumor cell programs

### Block 5 — TME intercellular communication
- CellChat object stratified by MPR vs. non-MPR (`CellChat`)
- Global signaling network comparison
- Focus on epithelial ↔ CD8 ligand-receptor axes

### Block 6 — Transcription factor activity
- TF activity inference on CD8 and epithelial compartments separately
- Network: CollecTRI (successor to DoRothEA, 12-source curated GRN)
- Method: `run_ulm` via `decoupleR`
- Differential TF activity: MPR vs. non-MPR

---

## Tools and references

| Tool | Purpose | Reference |
|---|---|---|
| Seurat | Preprocessing, clustering | Hao et al., *Cell* 2021 |
| Harmony | Batch correction | Korsunsky et al., *Nat Methods* 2019 |
| ProjecTILs | CD8 state annotation | Andreatta et al., *Nat Comm* 2021 |
| UCell | Gene module scoring | Andreatta & Carmona, *Comp Struct Biotech J* 2021 |
| scRepertoire | TCR clonotype analysis | Borcherding et al., *F1000Research* 2020 |
| CopyKAT | Malignant cell inference | Gao et al., *Nat Biotech* 2021 |
| CellChat | Cell-cell communication | Jin et al., *Nat Comm* 2021 |
| decoupleR | TF activity inference | Badia-i-Mompel et al., *Bioinformatics Advances* 2022 |
| CollecTRI | TF regulon network | Müller-Dott et al., *Nucleic Acids Research* 2023 |
| SingleR | Automated annotation | Aran et al., *Nat Immunol* 2019 |

---

## Progress

| Block | Status |
|---|---|
| Block 0 — Data acquisition |  Pending |
| Block 1 — QC and preprocessing |  Pending |
| Block 2 — Global TME annotation |  Pending |
| Block 3 — CD8 T cell focus |  Pending |
| Block 4 — Malignant epithelial focus |  Pending |
| Block 5 — TME communication (CellChat) |  Pending |
| Block 6 — TF activity (CollecTRI) |  Pending |

---

## Repository structure

```
TumorImmune_Crosstalk_LUAD/
├── README.md
├── Scripts/
│   ├── 00_data_download.R
│   ├── 01_QC_preprocessing.R
│   ├── 02_TME_annotation.R
│   ├── 03_CD8_ProjecTILs_UCell_TCR.R
│   ├── 04_Epithelial_CopyKAT_UCell.R
│   ├── 05_CellChat_MPR_NMPR.R
│   └── 06_TF_activity_CollecTRI.R
├── Data/
│   └── metadata_LUAD.csv
├── Figures/
└── Results/
```

---

## Relationship to prior work

This project is a direct continuation of [CD8_NSCLC_scRNAseq](https://github.com/yasmina-bioinfo/CD8_NSCLC_scRNAseq), which characterized CD8 T cell heterogeneity across two NSCLC datasets (GSE131907, GSE207422) and identified a CD8_Exhausted_Terminal enrichment in MPR patients (OR = 3.36, p_adj < 0.001) and a STAT2-high exhaustion program in non-MPR patients.

The present project extends this work by:
- Validating CD8 exhaustion findings in a larger independent cohort (n = 61 LUAD vs. n = 8)
- Adding the malignant epithelial compartment to investigate tumor-immune crosstalk
- Incorporating paired TCR sequencing to validate T cell state annotations
- Upgrading TF inference from DoRothEA to CollecTRI for improved regulon coverage

---

## Author

**Myriam Yasmina Soumahoro**   
[GitHub](https://github.com/yasmina-bioinfo)
