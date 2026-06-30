# TumorImmune_Crosstalk_LUAD

## Publication status 
> **Preprint in preparation** : Multi-cohort single-cell analysis of tumor immune crosstalk and anti-PD-1 response in NSCLC.
> Targer: bioRxiv : manuscript in preparation.

> **How does tumor cell plasticity and heterogeneity shape CD8 T cell functional states and their response to anti-PD-1 in lung adenocarcinoma?**

This project investigates the crosstalk between malignant epithelial cells and CD8 T cells in the tumor microenvironment (TME) of lung 
adenocarcinoma (LUAD), using neoadjuvant chemo-immunotherapy data with matched pathological response (MPR vs. non-MPR; MPR vs pCR and non-MPR 
vs pCR).
Although histological subtype was not explicitly reported by Hu et al. 2023, the presence of AT2 epithelial cells (Tumor_epithelial_AT2) in the 
TME annotation is consistent with a predominantly LUAD composition, as AT2 cells represent the recognized cell of origin of lung adenocarcinoma 
(Sainz de Aja et al. 2021; Xing et al. 2021). KRT5 and TP63 expression, canonical squamous markers, was restricted to the basal epithelial 
cluster and absent in AT2 cells, further supporting an adenocarcinoma histology. Cross-cohort comparison was therefore restricted to LUAD 
patients in GSE243013 (n=63), excluding LUSC and other NSCLC subtypes, to minimize histological confounding. LUAD and LUSC exhibit distinct 
tumor immune microenvironments, with LUAD showing higher CD8 T cell infiltration and exhaustion signatures, ensuring that the observed CD8 
exhaustion and TAM immunosuppressive signals reflect LUAD-specific biology.

This work is an extension to validate findings from [CD8_NSCLC_scRNAseq](https://github.com/yasmina-bioinfo/CD8_NSCLC_scRNAseq) in an 
independent, larger cohort with paired TCR sequencing.

---

## Dataset

## Dataset

| Field | Details |
|---|---|
| **Accession** | [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013) |
| **Publication** | Liu et al., *Cell* 2025 |
| **Histology** | LUAD only (n = 63, filtered from 234 NSCLC) |
| **Modalities** | scRNA-seq + scTCR-seq |
| **Treatment** | Neoadjuvant anti-PD-1 + chemotherapy |
| **Response variable** | Pathological response: MPR vs. non-MPR vs. pCR |
| **Timing** | Post-treatment surgical resection |

| Field | Details |
|---|---|
| **Accession** | [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) |
| **Publication** | Hu et al., *Cancer Cell* 2023 |
| **Histology** | NSCLC (predominantly LUAD, n = 13) |
| **Modalities** | scRNA-seq |
| **Treatment** | Neoadjuvant anti-PD-1 + chemotherapy |
| **Response variable** | Pathological response: MPR vs. NMPR |
| **Timing** | Post-treatment surgical resection |
---

## Biological question

Building on a CD8_Exhausted_Terminal enrichment signal in MPR patients (OR = 3.36, GSE207422), this project asks:

1. What transcriptional programs define malignant epithelial heterogeneity in LUAD under treatment pressure?
2. Do specific tumor cell states associate with impaired CD8 T cell function in non-MPR patients?
3. What ligand-receptor axes mediate epithelial–CD8 crosstalk differentially in MPR vs. non-MPR, MPR vs pCR and non-MPR and pCR?
4. Which transcription factors drive these programs in each compartment?

---

## Analytical pipeline

### Block 0 : Data acquisition
- Download GSE243013 (expression matrix + TCR metadata)
- Filter to LUAD samples (n = 63)
- Harmonize clinical metadata (MPR / non-MPR / pCR classification)

### Block 1 : QC and preprocessing
- Quality control, normalization, HVG selection (`Seurat`)
- Batch correction across patients (`Harmony`)
- Global TME clustering

### Block 2 : Global TME annotation
- Cluster marker identification (top markers per cluster)
- Manual annotation + automated validation (`SingleR`, `sctype`)
- Azimuth annotation pending (RAM constraint, server execution required)

### Block 3 : CD8 T cell focus
- Robust state annotation on reference atlas (`ProjecTILs`)
- Functional module scoring: exhaustion, cytotoxicity, memory (`UCell`)
- TCR integration: clonotype expansion, repertoire diversity (`scRepertoire`)
- TF activity inference: top 20 TFs by variance across CD8 states (`CollecTRI/decoupleR`)
- MPR vs. non-MPR, pCR vs MPR, and non_MPR vs pCR comparisons of CD8 state composition
- Cross-cohort validation: ProjecTILs and CollecTRI applied to GSE207422 CD8 T cells (planned)

### Block 4 : Immunosuppressive TME Compartment 
- TAM subtype annotation and reclustering (GSE243013 + GSE207422)
- Functional scoring: M2, M1, SPP1, IFN signatures (`UCell`)
- Patient-level pseudobulk analysis MPR vs pCR (`UCell`)
- TF activity inference on TAMs (`CollecTRI`) in PROGRESS 
- Malignancy epithelial analysis : CopyKAT, SCEVAN, CytoTRACE 2, UCell, CollecTRI (GSE207422 only)

### Block 5 : TME intercellular communication
- Intercellular communication inference across full TME using LIANA+ (consensus across CellChatDB, OmniPath and other LR databases)
- Differential cell-cell communication analysis across MPR vs. non-MPR vs. pCR using MultiNicheNet (multi-sample multi-condition design)
- Focus on TAM ↔ CD8, TAM ↔ Epithelial, CD8 ↔ Epithelial ligand-receptor axes
- Full TME opened to detect unexpected interactions, interpretation prioritized on axes established in Blocs 3 and 4

---

## Tools and references

| Tool | Purpose | Reference |
|---|---|---|
| Seurat | Preprocessing, clustering | Hao et al., *Cell* 2021 |
| Harmony | Batch correction | Korsunsky et al., *bioRxiv* 2019 |
| ProjecTILs | CD8 state annotation | Andreatta et al., *Nat Comm* 2021 |
| UCell | Gene module scoring | Andreatta & Carmona, *Comp Struct Biotech J* 2021 |
| scRepertoire | TCR clonotype analysis | Borcherding et al., *F1000Research* 2020 |
| CopyKAT | Malignant cell inference | Gao et al., *Nat Biotech* 2021 (via ResearchGate) |
| decoupleR | TF activity inference | Badia-i-Mompel et al., *Bioinformatics Advances* 2022 |
| CollecTRI | TF regulon network | Müller-Dott et al., *Nucleic Acids Research* 2023 |
| SingleR | Automated annotation | Aran et al., *Nat Immunol* 2019 |
| sctype | Automated cell type annotation | Ianevski et al., *Nature Communications* 2022 |
| Azimuth | Reference-based annotation (not installable on Windows, replaced by sctype) |  Hao et al., *bioRxiv* 2022 |
| SCEVAN | CNV inference (malignant vs normal epithelial) | De Falco et al., Nature Communications 2023 |
| CytoTRACE 2 | Differentiation potency inference (epithelial cells) | Kang et al., Nature 2024 |
| LIANA+ | Intercellular communication inference (consensus) | Dimitrov et al., *Nat Comm* 2022 |
| CellChat v2 | Differential cell-cell communication analysis MPR vs NMPR | Jin et al., *Nature Communications* 2021 |
---

## Progress

## Progress log

### Bloc 0 : Data acquisition and preparation
- Dataset: GSE243013, filtered to LUAD (n = 63 patients)
- Cell-level metadata: 336,685 LUAD cells
- Response groups: non-MPR (n=42), MPR (n=10), pCR (n=11)
- Clinical metadata harmonized and saved

### Bloc 1 : QC and preprocessing
- Cells before QC filtering: 336,685
- Cells after QC filtering: 298,567 (~11% removed)
- Thresholds: nFeature > 200, nCount 500–12,000, percent.mt < 7%
- Note: nFeature max = 4,000 for non-MPR, 3,000 for MPR/pCR
- Normalization: LogNormalize (scale factor = 10,000)
- HVG selection: 3,000 genes (vst method)
- PCA: 50 PCs computed, 40 retained (ElbowPlot inspection)
- Batch correction: Harmony by sampleID (inter-patient)
- Clustering: resolution = 0.5 and 20 clusters identified
- UMAP: generated and validated (Harmony integration confirmed)

### Bloc 2 : TME annotation

#### Script 01: FindAllMarkers complete with 20 clusters annotated manually (top50 markers)
    - 18/20 clusters annotated with high confidence
    - 2/20 clusters with medium confidence (Cluster 5: CD4 Tfh/exhausted, Cluster 9: NKT/γδ)
#### Script 02: SingleR complete (HumanPrimaryCellAtlas, label.main + label.fine)

#### Script 03: sctype complete (Lung tissue reference; limited resolution for immune subtypes)

#### Script 04: Azimuth, installed in WSL (Ubuntu 24.04, R 4.6.0) after 3h installation
  - lungref reference (584,884 cells) loaded successfully
  - RunAzimuth blocked by RAM constraint (16GB), server execution pending
#### Script 05: Consensus annotation complete with 20 clusters annotated, 3 pending Azimuth confirmation
  (Clusters 5, 9, 13, 18)

#### Script 06: TME visualization complete
  - UMAP final annotation + split by response
  - Barplot proportions (Chi-2 p < 2.2e-16)
  - Fisher post-hoc: MPR vs non-MPR, pCR vs non-MPR, MPR vs pCR; post-hoc pairwise comparisons saved in Results/Tables/Bloc2_fisher_posthoc.csv

### Bloc 3 : T cells analysis

**Methodological update (Bloc 3):** ProjecTILs and CollecTRI applied to GSE207422 for cross-dataset validation of CD8 exhaustion states and transcription factor activity in anti-PD-1 response [IN PROGESS]

> Note: Bloc 3 scripts without a dataset suffix correspond to GSE243013 cohort ( n=63 patients). 
> GSE207422 cohort (n= 13 patients) cross-validation scripts are explicitly labelled `Bloc3_GSE207422_`.

#### GSE243013
#### Script 01: T cells subsetting (clusters 1,2,3,4,5,9,10,11,17) : 172,110 cells
  ElbowPlot inspection = 30 PCs retained

#### Script 02: Harmony + clustering (resolution=0.4) + UMAP = 16 T cell clusters

#### Script 03: ProjecTILs annotation (human CD8 reference)
  - 108,456/172,110 cells (63%) filtered by scGate (non-CD8 pure)
  - 63,654 cells projected = 7 CD8 states identified:
    CD8.EM (25,575), CD8.CM (18,181), CD8.TEX (15,115), CD8.TPEX (2,011),
    CD8.NaiveLike (1,342), CD8.MAIT (1,273), CD8.TEMRA (157)
  - STACAS alignment failed (RAM) = direct projection used
  - CD8.TEX enriched in non-MPR, consistent with portfolio narrative

  - **METHODOLOGICAL NOTE (Scripts 4, 5, 6)**: ProjecTILs was run before canonical marker-based annotation.
  The recommended workflow is to first identify CD8 clusters by canonical markers, then project only confirmed CD8 cells onto the reference atlas.
  Running ProjecTILs on the full T cell subset (including CD4, Tregs, NKT, MAIT) caused scGate to incompletely filter non-CD8 populations, resulting in artifactual CD8 labels for non-CD8 clusters.
  **For future analyses: subset confirmed CD8 clusters first, then run ProjecTILs.**
  Final annotation relies on canonical markers for non-CD8 clusters.

#### Script 07 : CD8 T cells subsetting (corrected workflow)
  - **METHODOLOGICAL CORRECTION**: Script 03 ran ProjecTILs on full T cell subset (172,110 cells including CD4, Tregs, NKT) before canonical marker annotation.
  This caused scGate to incompletely filter non-CD8 populations (63% removed) and generated artifactual CD8 labels for non-CD8 clusters.
  - Corrected workflow: canonical marker annotation (Script 04) performed first to identify confirmed CD8 clusters (1, 2, 5, 6, 10), then ProjecTILs applied to CD8 pure subset only.
  - CD8 confirmed clusters: 75,622 cells
  - ElbowPlot inspection = 25 PCs retained (vs 30 for full T cell subset)

#### Script 08 : ProjecTILs on confirmed CD8 cells (corrected workflow)
  - CD8 input: 75,622 cells (clusters 1, 2, 5, 6, 10)
  - scGate filtering: 18,035/75,622 cells removed (24%) vs 63% in Script 03
  - Cells projected: 57,587
  - Runtime: ~42 minutes (15h33 - 16h15)
  - Note: Script 03 total time ~3h including reference download (~2h) + projection

| CD8 state | Script 03 (T cells full) | Script 08 (CD8 pure) |
|---|---|---|
| CD8.EM | 25,575 | 22,459 |
| CD8.CM | 18,181 | 16,944 |
| CD8.TEX | 15,115 | 15,675 |
| CD8.TPEX | 2,011 | 1,821 |
| CD8.NaiveLike | 1,342 | 446 |
| CD8.MAIT | 1,273 | 47 |
| CD8.TEMRA | 157 | 195 |

Key observation: CD8.TEX stable across both attempts (~15,000 cells), robust signal.
CD8.MAIT dropped from 1,273 to 47,  confirms cleaner CD8 population in Script 08.

#### Script 09 : CD8 ProjecTILs barplot (proportions by response)
  - Chi-2 p < 2.2e-16 : CD8 state distribution highly significantly different across groups
  - CD8.TEX: non-MPR (~30%) > MPR (~25%) > pCR (~15%): gradient confirmed statistically
  - CD8.CM: pCR (~40%) = MPR (~35%) > non-MPR (~30%) : memory enriched in responders
  - CD8.EM: similar across groups (~35-40%) : functional status to be assessed by UCell
  - CD8.TPEX: visible in pCR only, consistent with reactivation hypothesis
  - CD8.MAIT: 47 cells (~0.08%), negligible, not visible on barplot. Reduction from 1,273 (Script 03) to 47 (Script 08) confirms cleaner CD8 population
  - Fisher post-hoc results saved in Results/Tables/Bloc3_CD8_fisher_posthoc.csv

#### Script 10 : UCell scoring on CD8 T cells
  **Exhaustion scores (mean), CD8.TEX:**
    - non-MPR = 0.347 > pCR = 0.249 > MPR = 0.220

  **Exhaustion scores (mean), CD8.TPEX:**
    - non-MPR = 0.268 > pCR = 0.225 > MPR = 0.197

#### Script 11: scRepertoire TCR analysis (WSL, R 4.6.0)
  - NOTE: scRepertoire requires gsl >= R 4.5.0, not available on Windows R 4.4
  Installed and run in WSL environment
  - TCR file: GSE243013_T_with_TCR_annotation.csv.gz (434,458 cells)
  - Matched to CD8 object: cells with TCR data in our CD8 subset (57,587 cells)

  **Clonotype expansion by CD8 state:**
  - CD8.TEX: 85.4% expanded, highest, confirms tumor-reactive identity
  - CD8.TPEX: 79.4% expanded, also tumor-reactive
  - CD8.EM: 70.8% expanded
  - CD8.NaiveLike: 16.5% expanded, coherent with naive/quiescent state

  **Clonotype expansion by pathological response:**
  - MPR: ~70% expanded
  - non-MPR: ~73% expanded
  - pCR: ~65% expanded, more non-expanded clones, higher polyclonality

  **Clonal diversity (mean per patient):**
  - MPR: 0.490 ± 0.184
  - pCR: 0.458 ± 0.118
  - non-MPR: 0.438 ± 0.168
  - Responders show higher clonal diversity than non-responders

#### Script 12 : CollecTRI TF activity on CD8 T cells (GSE243013)
  - Tool: decoupleR run_ulm + CollecTRI network (43,159 interactions, 1,186 TFs)
  - RAM constraint: restricted to CD8.TEX, CD8.TPEX, CD8.EM, max 10,000 cells/state
  - Total analyzed: ~21,821 cells, set.seed(42) for reproducibility
  - Top 20 TFs by variance: HSF1, HSF2, RFXAP, RFXANK, ELK4, MYC, RFX5, NFYC, NFKB, RLF, DMTF1, RELA, CIITA, NFYB, MLXIP, JUN, HOPX, NFKB1, DAXX, TBX21
  - Violin plots: Top 6 by variance (objective, no confirmation bias) Top 6: HSF1, HSF2, RFXAP, RFXANK, ELK4, MYC
    - ELK4 rank 5/772, enriched MPR CD8.TEX 
    - TBX21 rank 20/772 , present MPR TEX WITH ELK4 (functional) and non-MPR TEX WITHOUT ELK4 (abortive)
    - STAT2 rank 59/772 , not in Top 20, IFN program confirmed by UCell only

  **Heatmap state-specific observations:**
  - CD8.TEX non-MPR : TBX21 dominant WITHOUT ELK4 → abortive cytotoxic program
  - CD8.TPEX non-MPR : MYC, NFKB1, RELA, HSF1/2 → chronic proliferation
  - CD8.TEX MPR : TBX21 + ELK4 + MHC II program → functional cytotoxic program
  - CD8.TPEX MPR : ELK4 dominant, TBX21 quasi-absent → pure reactivation program
  - CD8.TEX pCR : HOPX dominant, ELK4 absent → post-response quiescence

#### Script 13 : Pseudobulk UCell scores MPR vs pCR / CD8 T cells
  - Input: Objects/Bloc3_11_seu_CD8_UCell.rds
  - Patient-level comparison: MPR (n=10) vs pCR (n=11); non-MPR excluded
  - Pseudobulk: mean UCell score per patient per signature
  - Wilcoxon rank-sum test (patient-level)
  - All signatures ns (p>0.05) : Cytotoxicity p=0.5573, Exhaustion p=0.5116, Memory p=0.1734, TPEX p=0.6047
  - Inter-patient heterogeneity similar between MPR and pCR, contrasts with TAMs where pCR showed greater heterogeneity than MPR (Bloc4A Script 06)

#### GSE207422 

#### Script 01 : ProjecTILs on annotated CD8 object (08_CD8_MPR_NMPR.rds)
  - Input: GSE207422 portfolio object (pre-annotated, Harmony already applied)
  - scGate filtering: 24% removed (consistent with GSE243013 corrected workflow)
  - 7 CD8 states identified: CD8.TEX, CD8.TPEX, CD8.EM, CD8.CM, CD8.NaiveLike, CD8.TEMRA, CD8.MAIT
  - CD8.TEX dominant in both MPR and NMPR (~55%), proportions more similar than GSE243013
  - CD8.TPEX slightly enriched in MPR , consistent with reactivation hypothesis

#### Script 02 : CD8 ProjecTILs barplot (proportions by response)
  - Chi-2 p = 4.37e-13 : CD8 state distribution significantly different between MPR and NMPR
  - CD8.TEX : MPR (~55%) ≈ NMPR (~54%) , no major proportional difference
  - CD8.TEMRA : enriched in NMPR, terminal differentiation, non-reactivable
  - CD8.TPEX : slightly enriched in MPR , reactivation capacity
  - Note: higher p-value vs GSE243013 (p < 2.2e-16) reflects smaller cohort (n=13) and subtler proportional differences

#### Script 03 : UCell scoring on CD8 T cells
  **Exhaustion scores by response (all CD8):**
  - NMPR > MPR : deeper exhaustion in non-responders, consistent with GSE243013

  **Cytotoxicity scores:**
  - MPR slightly higher , residual effector capacity preserved

  **TPEX_UCell scores:**
  - MPR > NMPR : precursor program more active in responders

  **Memory_UCell scores:**
  - NMPR > MPR : quiescence/anergy rather than functional memory

  **Focus CD8.TEX and CD8.TPEX:**
  - CD8.TEX MPR : high exhaustion AND high cytotoxicity scores confirms intra-TEX heterogeneity (co-expression exhaustion/effector) 
  - CD8.TPEX MPR : higher TPEX_UCell than NMPR; preserved plasticity in responders 
  - Simultaneous elevation of TEX and TPEX scores in MPR confirms reactivation hypothesis across both datasets

#### Script 04 : CollecTRI TF activity on CD8 T cells (GSE207422)
  - Tool: decoupleR run_ulm + CollecTRI network (local CSV)
  - Subset: CD8.TEX, CD8.TPEX, CD8.EM, max 10,000 cells/state, set.seed(42)
  - Top 20 TFs by variance: MYC, NFKB, JUN, SP1, HSF1, E2F4, E2F1, STAT1, RFXAP, RFXANK, ELK4, RFX5, HSF2, ABL1, SRSF2, TFDP1, NFYC, ZBTB4, CIITA, DOT1L
  - Violin plots: Top 6 by variance (objective, no confirmation bias)
    Top 6: MYC, NFKB, JUN, SP1, HSF1, E2F4
    - ELK4 rank 11/719 — present Top 20, signal weaker than GSE243013
    - STAT2 absent Top 20 both datasets : IFN program confirmed by UCell only

  **Heatmap state-specific observations:**
  - CD8.TEX NMPR : MHC II program (RFXAP/RFXANK/RFX5/CIITA) + STAT1 dominant
  - CD8.TPEX NMPR : MYC, E2F1, E2F4, SRSF2, TFDP1 → chronic proliferation 
  - CD8.TEX MPR : ELK4 dominant + NFYC + ABL1 + DOT1L
  - CD8.TPEX MPR : ELK4 + ABL1 + NFYC → ELK4-driven reactivation 

### Bloc 4 : Immunosuppressive TME compartment / TAMs and malignant epithelial cells

**Methodological note (Bloc 4 CollecTRI):** Following correction applied to 
Bloc 3 CollecTRI CD8 analysis, violin plots for TAM CollecTRI will display 
Top 6 TFs by cross-cell variance (objective selection, no confirmationbias). 
Scripts Bloc4A_07 and Bloc4B_06 to be updated accordingly. 

  ### Bloc 4A : GSE243013
#### Script 01 : TAMs subsetting and ElbowPlot
  - Subset from cluster 7 (M2-like/immunosuppressive) of global TME annotation
  - 14,901 cells extracted
  - Normalization + HVG (2000) + Scaling + PCA (30 PCs)
  - ElbowPlot inspection, 20 PCs retained

#### Script 02 : Harmony + Clustering + UMAP (TAMs)
  - dims = 1:20 based on ElbowPlot inspection (Script 01)
  - Resolution = 0.1 , 8 clusters identified
  NOTE: resolution 0.4 gave 19 clusters (over-fragmentation)
  resolution 0.1 gives 8 biologically meaningful TAM subtypes
  - Harmony batch correction by sampleID

#### Script 03 : FindAllMarkers TAMs
  - Top 20 markers per cluster (vs top 50 for TME, TAMs more homogeneous)
  - Wilcoxon test, min.pct = 0.25, logFC ≥ 0.25, only.pos = TRUE
  - Output: Results/Markers/ (à préciser selon ton dossier)
  - See TAM_markers_biological_roles.md for full annotation

#### Script 04 : TAMs final annotation + UMAP + Barplot
  - 7 TAM subtypes annotated (cluster 6 excluded : lymphocyte contamination)
  - Chi-2 p < 2.2e-16, TAM subtype distribution highly significantly different across groups

   **Key observations:**
     - Stress-response TAMs (MARCO+/PPARG+/HSP-high): enriched in non-MPR, consistent with hypoxic/metabolically stressed immunosuppressive TME
     - Tissue-resident M2-like TAMs: enriched in non-MPR, chronic anti-inflammatory program blocking CD8 T cell function
     - IFN-stimulated TAMs (PD-L1+/IDO1+): enriched in non-MPR, direct CD8 suppression via PD-L1/PD-1 axis and tryptophan depletion (IDO1) = triple TAM-mediated immunosuppression in non-MPR supports H2

#### Script 05 : UCell scoring on TAMs (GSE243013)
- Input: Bloc4A_04_seu_TAMs_annotated.rds (7 TAM subtypes)
- Signatures: M2_immunosuppressive (MRC1, CD163, TGFB1, IL10, VEGFA, CD274, IDO1, CSF1R), M1_inflammatory (TNF, IL1B, IL6, CXCL10, NOS2), SPP1_signature (SPP1, GPNMB, APOE, TREM2),IFN_response (ISG15, IFIT1, IFIT3, CXCL9, CXCL10)
- Pairwise Wilcoxon tests with Bonferroni correction (3 comparisons per signature)
  Bonferroni threshold: p < 0.0167
  References: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)
- Short labels used for violin subtype readability

#### Script 06 : Pseudobulk UCell scores MPR vs pCR (patient-level)
- Input: Bloc4A_05_seu_TAMs_UCell.rds
- Subset: MPR (n=10 patients, 867 cells) and pCR (n=11 patients, 1019 cells)
  NOTE: non-MPR excluded, focus on MPR vs pCR to test intermediate state hypothesis
- Pseudobulk: mean UCell score per patient per signature
  NOTE: aggregation to patient-level avoids pseudo-replication
- Wilcoxon rank-sum test MPR vs pCR (patient-level)
  References: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)

#### Script 07 : CollecTRI TF activity on TAMs (GSE243013)
  - Tool: decoupleR run_ulm + CollecTRI network (local CSV)
  - All 7 TAM subtypes included in full (range: 248–4,316 cells/subtype, set.seed(42))
  - Top 20 TFs by variance across cells
  - Short labels applied for readability (consistent with Script 05)
  - Violin plots: Top 6 by variance (objective, no confirmation bias)
    Top 6: RFXAP, HSF1, RELA, RFXANK, NFKB, STAT1
  - Output: 3 heatmaps (MPR/non-MPR/pCR) + key TFs violin
  - Objects saved: Bloc4A_07_seu_TAMs_TF.rds

  **Heatmap state-specific observations:**
  - IFN-stimulated : STAT1 constitutive all conditions : not discriminating
  - Resident M2 MPR : ELK4 dominant : pro-immunogenic signal
  - Resident M2 non-MPR : HIF1A, NFKB, stress program dominant
  - Resident M2 pCR : ELK4 + DMTF1
  - Monocyte FCN1+ : HIF1A enriched non-MPR : hypoxic adaptation
  - Stress-response non-MPR : MYC + HSF1/HSF2 dominant 
  - Proliferating MPR/non-MPR : KAT6B dominant : chromatin remodeling
  - Proliferating pCR : MYC very strong
  - Classical-Mono MPR : ELK4 + MHC II dominant 
  - Classical-Mono pCR : ELK4 + MHC II present but weaker
  - LAMs MPR : KAT6B + AEBP1
  - LAMs non-MPR : MHC II + KAT6B + AEBP1 : heterogeneous

### Bloc 4B : GSE207422 TAMs and Epithelial cells

#### Script 01 : TAMs extraction + UMAP + Barplot (initial annotation)
- TAM populations extracted from TME object (04_TME_MPR_NMPR.rds)
  NOTE: annotation already performed during global TME annotation (TME_cell_type column)
- Three TAM subtypes identified: TAM_like (n=7,969),TAM_like_MRC1 (n=5,229), TAM_like_SPP1 (n=362)
- NOTE: TAM_like_SPP1 small population (n=362), interpret with caution
- Chi-2 p < 2.2e-16

#### Script 02 : UCell scoring on TAMs (initial : TAM_like not yet reclustered)
- Signatures: M2_immunosuppressive (MRC1, CD163, TGFB1, IL10, VEGFA, CD274, IDO1, CSF1R), M1_inflammatory (TNF, IL1B, IL6, CXCL10, NOS2), SPP1_signature (SPP1, GPNMB, APOE, TREM2),IFN_response (ISG15, IFIT1, IFIT3, CXCL9, CXCL10)
- References: Chen et al. 2021 (PMC8053174), Italiani & Boraschi 2019 (PMC6543837)
- NOTE: UCell will be rerun after TAM_like reclustering (Script 05) for complete interpretation

#### Script 03 : TAM_like reclustering and annotation
- TAM_like (n=7,969) reclustered, initial annotation insufficiently granular for biological interpretation
- Preprocessing: 2000 HVGs, 15 PCs retained (stdev inspection), resolution = 0.3
- 8 clusters identified, cluster 7 excluded (T cell contamination: CD3G, CD3D, TRAC, GZMA)
- Top 20 markers per cluster (FindAllMarkers, Wilcoxon)
- 7 TAM_like subtypes annotated:
  | Cluster | Annotation |
  |---|---|
  | 0 | TAM_like_resident_M2 (iron metabolism/anti-inflammatory) |
  | 1 | TAM_like_IFN (PD-L1+/IDO1+/CXCL9+) |
  | 2 | TAM_like_monocyte (classical inflammatory) |
  | 3 | TAM_like_lipid (CCL18+/AKR+) |
  | 4 | TAM_like_stress (HSP-high/M1-like) |
  | 5 | TAM_like_regulatory (glucocorticoid-responsive) |
  | 6 | TAM_like_M2 (SIGLEC8+/CCL18+) |
  | 7 | EXCLUDE : T cell contamination |

#### Script 04 : Combined TAM annotation + UMAP + Barplot
- TAM_like subclusters (Script 03) merged with TAM_like_MRC1 and TAM_like_SPP1
- TAM_like_MRC1 and TAM_like_SPP1 kept as independent populations, not merged with TAM_like subclusters
- Final: 9 TAM subtypes analyzed (T cell contamination excluded)
- Chi-2 p < 2.2e-16

#### Script 05 : UCell scoring on TAMs combined (post-reclustering)
- Input: Bloc4B_04_seu_TAMs_combined.rds (9 TAM subtypes)
- Signatures: M2_immunosuppressive (MRC1, CD163, TGFB1, IL10, VEGFA, CD274, IDO1, CSF1R), M1_inflammatory (TNF, IL1B, IL6, CXCL10, NOS2),SPP1_signature (SPP1, GPNMB, APOE, TREM2), IFN_response (ISG15, IFIT1, IFIT3, CXCL9, CXCL10)
- Wilcoxon test computed manually on cell-level scores (sapply)
  NOTE: stat_compare_means (ggpubr) tested but returns p=1 on aggregated data
  Solution: manual Wilcoxon + geom_text annotation on barplot
- References: Chen et al. 2021 (PMC8053174), Italiani & Boraschi 2019 (PMC6543837)

#### Script 06 : CollecTRI TF activity on TAMs (GSE207422)
  - Tool: decoupleR run_ulm + CollecTRI network (local CSV)
  - All 9 TAM subtypes included in full (range: 171–5,229 cells/subtype, set.seed(42))
  - Top 20 TFs by variance across cells
  - Short labels applied for readability (consistent with Bloc4B UCell script)
  - Violin plots: Top 6 by variance (objective, no confirmation bias)
    Top 6: RFXAP, RFXANK, HSF1, RFX5, CIITA, RELA
  - Output: 2 heatmaps (MPR/NMPR) + key TFs violin
  - Objects saved: Bloc4B_06_seu_TAMs_TF.rds

  **Heatmap state-specific observations:**
  - IFN-stimulated MPR : IRF1, STAT1, IRF5, REL dominant : strong IFN program
  - MRC1+ M2-like MPR : MHC II dominant (CIITA, RFX5, RFXANK, RFXAP) 
  - SPP1+ immunosuppressive NMPR : highly active, heterogeneous : complex immunosuppressive program 
  - Monocyte-derived NMPR : HIF1A dominant : hypoxic adaptation 
  - Stress-response : active both conditions : constitutive program

#### Script 08 : CNV inference on epithelial cells

**CopyKAT, first attempt (Bloc4B_08_CopyKAT.R):**
- Tool: CopyKAT v1.2.5 (Gao et al., Nature Genetics 2021)
- Normal reference: Ciliated_epithelial (565 cells)
- Full epithelial object (8,944 cells), win.size=15, LOW.DR=0.05, UP.DR=0.1, n.cores=1
- RAM constraint: process killed at step 4/7 despite WSL 13GB and win.size reduction.
Kaggle: R package installation blocked (no internet access).
Google Colab (12GB): session crashed. Full dataset analysis not feasible.
- Script retained for methodological transparency

**CopyKAT by condition : revised approach (Scripts 08_MPR, 08_NMPR, 08h):**
- Methodological revision: initial per-subtype runs (Scripts 08d-08g) were replaced by
  per-condition runs for cross-tool consistency with SCEVAN
- MPR run: all epithelial MPR cells (1,612 cells) + Ciliated reference
- NMPR run: per epithelial subtype x NMPR condition due to RAM constraints (7,332 cells total):
  AT2_NMPR, Basal_NMPR, EMT_NMPR, Tumor_epithelial_NMPR; each run 500-2,000 cells
- Same parameters: win.size=15, LOW.DR=0.05, UP.DR=0.1, n.cores=1
- Original per-subtype scripts (08d-08g) retained for methodological transparency

**SCEVAN (Scripts 08a-08b), cross-validation:**
- Tool: SCEVAN v1.0.3 (De Falco et al., Nature Communications 2023)
- Normal reference: Ciliated_epithelial (565 cells)
- MPR and NMPR analyzed separately, par_cores=1, SUBCLONES=FALSE, ClonalCN=FALSE
- After internal QC: 8,405 cells and 8,162 genes retained

**SCEVAN vs CopyKAT : cross-validation by condition (Ciliated cells excluded):**

| Condition | SCEVAN tumor | CopyKAT aneuploid | Concordance |
|-----------|-------------|-------------------|-------------|
| MPR | ~32% | ~38% | Consistent |
| NMPR | ~47% | ~75% | Consistent |
| MPR vs NMPR | NMPR > MPR | NMPR > MPR | Directional concordance |

**Key observations:**
- Both tools agree on directional difference: higher tumor/aneuploid fraction in NMPR than MPR, consistent with clinical non-response
- CopyKAT systematically detects higher aneuploid fractions than SCEVAN, reflecting its higher sensitivity to low-level copy number fluctuations, a known property of the tool
- SCEVAN proportions (~32% MPR, ~47% NMPR) are more conservative and better aligned with expected residual tumor burden after neoadjuvant treatment
- MPR CopyKAT run shows higher not.defined fraction (~25%) compared to NMPR, reflecting the global run design vs per-subtype design, a methodological asymmetry documented here
- UMAP spatial concordance: both tools identify the same tumor cell cluster in the lower UMAP region, confirming localization of malignant epithelial cells
- CNV clonal profiles (SCEVAN ClonalCNProfile) are identical between MPR and NMPR, suggesting treatment response is not driven by chromosomal architecture differences but by transcriptional and microenvironmental factors
- EMT and Tumor_epithelial subtypes exclusive to NMPR, transcriptional rather than genomic drivers of non-response

**Previous per-subtype cross-validation (Scripts 08d-08g, retained for transparency):**

| Cell type | SCEVAN tumor | SCEVAN normal | CopyKAT aneuploid | CopyKAT diploid |
|-----------|-------------|---------------|-------------------|-----------------|
| Tumor_epithelial | 3,359 (99%) | 10 (0.3%) | 3,369 (99%) | 0 |
| Tumor_epithelial_AT2 | 0 (0%) | 1,661 (95%) | 1,642 (94%) | 4 (0.2%) |
| Tumor_epithelial_basal | 533 (35%) | 689 (45%) | 744 (49%) | 460 (30%) |
| Tumor_epithelial_EMT | 0 (0%) | 1,600 (94%) | 643 (38%) | 940 (55%) |

#### Script 09a : CytoTRACE2 on epithelial cells (SCEVAN labels)
- Tool: CytoTRACE2 v1.1.0 (Kang et al., Nature Methods 2025)
- Input: Bloc4B_07_seu_Epithelial.rds (n=8,944 cells)
- SCEVAN labels used as annotation (primary CNV tool)
- 6 groups: tumor_MPR (n=506), normal_MPR (n=493), tumor_NMPR (n=3,386), normal_NMPR (n=1,867), EMT_NMPR (n=1,707), Ciliated (n=565)
- Parameters: species="human", ncores=1, seed=42, batch_size=3000, smooth_batch_size=1000
- Batching used to manage RAM constraints (8,944 cells, 24,292 genes)
- Score range: 0 (differentiated) to 1 (totipotent)
- Output: Bloc4B_09a_CytoTRACE2_SCEVAN_scores.csv
          Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds
          Bloc4B_09a_UMAP_CytoTRACE2_Score.png
          Bloc4B_09a_Violin_CytoTRACE2_SCEVAN.png

#### Script 09b : CytoTRACE2 on epithelial cells (CopyKAT labels)
- Tool: CytoTRACE2 v1.1.0 (Kang et al., Nature Methods 2025)
- Input: Bloc4B_09a_seu_Epithelial_CytoTRACE2.rds
- CopyKAT labels used for cross-validation of EMT classification
- 6 groups: aneuploid_MPR, diploid_MPR, aneuploid_NMPR, diploid_NMPR, EMT_NMPR, Ciliated
- CytoTRACE2 scores reused from Script 09a, no recomputation
- Output: Bloc4B_09b_CytoTRACE2_CopyKAT_scores.csv
          Bloc4B_09b_seu_Epithelial_CytoTRACE2_CopyKAT.rds
          Bloc4B_09b_UMAP_CytoTRACE2_Score.png
          Bloc4B_09b_Violin_CytoTRACE2_CopyKAT.png

#### Script 09c : Final epithelial group classification
- EMT_NMPR reclassified as normal_NMPR based on:
  1. SCEVAN highest specificity (0.75) among CNV tools (Lanucara et al., Biomedicines 2024)
  2. CytoTRACE2 low stemness scores under both SCEVAN and CopyKAT labels
  3. CopyKAT tendency to overestimate tumor fractions in benchmarks
- Final 5 groups: tumor_MPR (n=506), normal_MPR (n=493), tumor_NMPR (n=3,386),normal_NMPR (n=3,574), Ciliated (n=565)
- Output: Bloc4B_09c_Final_groups.csv
          Bloc4B_09c_seu_Epithelial_FinalGroups.rds
          Bloc4B_09c_Violin_CytoTRACE2_FinalGroups.png

#### Script 10 : UCell scoring on epithelial cells

- Tool: UCell v2.x (Andreatta & Carmona, iScience 2021)
- Input: Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds (n=8,944 cells)
- 5 final groups: tumor_MPR (n=506), normal_MPR (n=493), tumor_NMPR (n=3,386), normal_NMPR (n=3,574 including reclassified EMT), Ciliated (n=565)
- 13 signatures: 11 MSigDB Hallmark + 2 custom (HSF1_targets, Antigen_presentation)
- Signatures imported via msigdbr package (Liberzon et al., Cell Systems 2015)
- Output: Results/Figures/BLOC4B_Epithelial_TAMs/UCell_Epithelial/
          Results/Tables/Bloc4B_10_UCell_Epithelial_scores.csv
          Objects/Bloc4B_10_seu_Epithelial_UCell.rds
- Figures: 1 dotplot (main figure) + 13 violin plots (supplementary)

#### Script 11 : CollecTRI TF activity on epithelial cells

- Tool: decoupleR run_ulm (Badia-i-Mompel et al., Bioinformatics Advances 2022)
- Network: CollecTRI (Müller-Dott et al., Nucleic Acids Research 2023), 43,159 interactions, 1,186 TFs
- Input: Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds (n=8,944 cells)
- 5 final groups: tumor_MPR (n=506), normal_MPR (n=493), tumor_NMPR (n=3,386), normal_NMPR (n=3,574), Ciliated (n=565)
- No downsampling applied (manageable dataset size)
- Top 20 TFs by cross-cell variance (objective selection, no confirmation bias)
- Top 6 TFs: MYC, RFXAP, RFXANK, HSF2, RFX5, CIITA
- Output: Results/Figures/BLOC4B_Epithelial_TAMs/CollecTRI_Epithelial/
          Results/Tables/Bloc4B_11_CollecTRI_Epithelial_TF_activity.csv
          Objects/Bloc4B_11_seu_Epithelial_TF.rds
- Figures: 1 heatmap (main figure) + violin plots split tumor/normal (supplementary)

## Bloc 5 : TME Intercellular Communication

### Biological question
Does the transcriptional divergence between MPR and NMPR/non-MPR compartments (CD8, TAMs, epithelial) reflect distinct intercellular communication programs that reinforce or undermine anti-PD-1 response? Specifically, do responders and non-responders conditions differ in the directionality, specificity, and coordination of bidirectional communication loops between CD8 T cells, TAMs, and epithelial cells, and does this bidirectional crosstalk constitute a self-reinforcing immunosuppressive network in non-responders?

### Cells of interest — guided by Blocs 3 and 4 findings

Intercellular communication analysis is guided by discriminant cell states identified in Blocs 3 and 4 (CollecTRI/decoupleR TF activity and UCell functional scoring). Cell subtypes of interest were selected based on differential transcriptional programs across response conditions. Full TME was opened for inference; interpretation is prioritized on these axes without excluding unexpected interactions.

**CD8 T cells — both datasets:**
- CD8.TEX and CD8.TPEX : most discriminant between MPR and NMPR/non-MPR (ELK4 enriched MPR cross-dataset, exhaustion program enriched non-MPR)
- CD8.EM, CD8.CM, CD8.TEMRA, CD8.NaiveLike : secondary context

**TAMs — GSE207422:**
- MPR priority : IFN-stimulated (IRF1/STAT1), Lipid-associated (HIF1A)
- NMPR priority : Stress-response (HIF1A/NF-KB), SPP1+ immunosuppressive (NF-KB/STAT3), Monocyte-derived (NF-KB/RELA)

**TAMs — GSE243013:**
- MPR priority : IFN-stimulated (IRF1/STAT1/REL), Resident M2 (ELK4)
- non-MPR priority : Resident M2 (HIF1A/NF-KB), Monocyte FCN1+ (HIF1A/NF-KB),Stress-response (MYC/NF-KB)

**Epithelial — GSE207422 only:**
- Excluded from main preprint narrative (single dataset)
- Retained in perspectives : tumor_NMPR shows enriched oncogenic/stress program (E2F1/E2F4/MYC/HIF1A) distinct from immune compartment signals

**Cross-dataset observation:**
HSF1 transcriptional activity is detected in TAMs and CD8 T cells across both conditions and both datasets, indicating a constitutive stress response program not exclusive to non-responders in immune compartments. HSF1 discriminant signal is restricted to tumor epithelial cells in GSE207422 (single dataset, not cross-validated). The intercellular link between HSF1 epithelial activity and immune compartment programs remains to be resolved, a question amenable to spatial transcriptomics addressed in perspectives.

### Dataset coverage
GSE207422 : CD8 + TAMs + Epithelial (tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, Ciliated) , full bidirectional analysis across all three compartments
- GSE243013 : CD8 + TAMs (no epithelial compartment available) , bidirectional CD8 ↔ TAMs analysis; pCR retained in analysis, excluded from main narrative

### Tools
- LIANA+ v0.1.14 (Dimitrov et al., Nature Communications 2022), consensus LR inference (sca + natmi + connectome)
- CellChat v2 (Jin et al., Nature Communications 2021), differential communication analysis

### Design rationale
Full TME opened for inference: interpretation prioritized on bidirectional CD8 ↔ TAMs, CD8 ↔ Epithelial, TAMs↔Epithelial axes established in Blocs 3 and 4. Both directions of each axis are analyzed to characterize communication loops, not unidirectional signals. Unexpected interactions retained only if they corroborate findings from prior blocs.

---

### Script 01 : LIANA+ GSE207422, intercellular communication inference

- Input: Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds + Bloc4B_04_seu_TAMs_combined.rds + 
  Bloc4B_09c_seu_Epithelial_FinalGroups.rds
- Merged TME object: 32,813 cells (CD8 n=7,423 | TAMs n=18,073 | Epithelial n=7,317)
- Cell types removed: CD8.MAIT (< 5 cells)
- Methods: sca + natmi + connectome aggregated via liana_aggregate()
- Resource: Consensus (CellChatDB + OmniPath + others)
- Conversion: SingleCellExperiment (SCE) for Seurat v5 compatibility
- Per condition runs: MPR (n=5,892 cells) + NMPR (n=26,921 cells)
- Output: Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_aggregated.csv
          Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv
          Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_NMPR_aggregated.csv
          Objects/Bloc5_01_seu_TME_GSE207422.rds
          Results/Figures/BLOC5_Communication/GSE207422/Bloc5_01_LIANA_dotplot_CD8_to_TME.png
- Supplementary figure: top 20 interactions from CD8.TEX and CD8.TPEX toward TAM_like_IFN and TAM_like_stress (selected based on Bloc 4A/4B findings as most discriminant TAM subtypes between MPR and NMPR) and all epithelial groups.Target selection was biologically motivated by prior results, not exhaustive.

### Script 02 : LIANA+ GSE207422, main axis figures (supplementary)

- Input: Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv + NMPR_aggregated.csv
- 6 bidirectional figures (MPR vs NMPR comparative dotplots), full TME:
  1. CD8 → TAMs (Bloc5_02_LIANA_CD8_TAMs.png)
  2. TAMs → CD8 (Bloc5_02_LIANA_TAMs_CD8.png)
  3. CD8 → Epithelial (Bloc5_02_LIANA_CD8_Epithelial.png)
  4. Epithelial → CD8 (Bloc5_02_LIANA_Epithelial_CD8.png)
  5. TAMs → Epithelial (Bloc5_02_LIANA_TAMs_Epithelial.png)
  6. Epithelial → TAMs (Bloc5_02_LIANA_Epithelial_TAMs.png)
- 1 supplementary figure: all compartments (Bloc5_01_LIANA_dotplot_CD8_to_TME.png)
- TAM short labels applied for readability (consistent with Blocs 4A/4B)
- Epithelial labels simplified: tumor_MPR + tumor_NMPR → Tumor epithelial; normal_MPR + normal_NMPR → Normal epithelial
- Top 5 interactions per target group (aggregate_rank <= 0.05)
- All 9 TAM subtypes + all CD8 subtypes + all epithelial groups included

### Script 02b : LIANA+ GSE207422, preprint figures - cells of interest

- Input: Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv + NMPR_aggregated.csv
- Same 6 bidirectional axes as Script 02 but filtered on cells of interest guided by Blocs 3 and 4 CollecTRI findings
- CD8 priority: CD8.TEX, CD8.TPEX
- TAMs priority: IFN-stimulated, Stress-response, SPP1+ immunosuppressive, Monocyte-derived, Lipid-associated
- Epithelial: Tumor epithelial + Normal epithelial + Ciliated (all kept)
- Output: Results/Figures/BLOC5_Communication/GSE207422/Preprint/
  1. Bloc5_02b_LIANA_CD8_TAMs_preprint.png
  2. Bloc5_02b_LIANA_TAMs_CD8_preprint.png
  3. Bloc5_02b_LIANA_CD8_Epithelial_preprint.png
  4. Bloc5_02b_LIANA_Epithelial_CD8_preprint.png
  5. Bloc5_02b_LIANA_TAMs_Epithelial_preprint.png
  6. Bloc5_02b_LIANA_Epithelial_TAMs_preprint.png
- Top 5 interactions per target group (aggregate_rank <= 0.05)

### Script 03 : CellChat v2 GSE207422, differential communication analysis

- Input: Bloc5_01_seu_TME_GSE207422.rds
- CellChatDB.human full database
- Per condition objects: cellchat_MPR.rds + cellchat_NMPR.rds
- Parameters: nboot = 100, min.cells = 5
- Merged object: Bloc5_03_cellchat_merged_GSE207422.rds
- Differential analysis: netVisual_diffInteraction() + rankNet() per axis
- Output: Bloc5_03_cellchat_MPR.rds
          Bloc5_03_cellchat_NMPR.rds
          Bloc5_03_cellchat_merged_GSE207422.rds
          Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_MPR_interactions.csv
          Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_NMPR_interactions.csv
- Note: Stress-response and Monocyte-derived TAM subtypes absent from CellChat 
output, likely excluded by min.cells = 5 threshold or no significant interactions detected. Present in LIANA+ only.
- Note: MultiNicheNet abandoned because n=3 MPR patients insufficient for pseudobulk differential analysis. Script conserved with explanatory note.

### Script 04 : CellChat v2 GSE207422, custom figures (supplementary)

- Input: Bloc5_03_CellChat_GSE207422_MPR_interactions.csv + NMPR_interactions.csv
- Custom ggplot2 figures from exported CellChat CSV (not netVisual_bubble)
- 6 bidirectional figures — full TME (supplementary):
  1. CD8 → TAMs (Bloc5_04_CellChat_CD8_TAMs.png)
  2. TAMs → CD8 (Bloc5_04_CellChat_TAMs_CD8.png)
  3. CD8 → Epithelial (Bloc5_04_CellChat_CD8_Epithelial.png)
  4. Epithelial → CD8 (Bloc5_04_CellChat_Epithelial_CD8.png)
  5. TAMs → Epithelial (Bloc5_04_CellChat_TAMs_Epithelial.png)
  6. Epithelial → TAMs (Bloc5_04_CellChat_Epithelial_TAMs.png)
- Heatmaps: number of interactions MPR vs NMPR (Bloc5_04_CellChat_heatmap_MPR.png + NMPR.png)
- Barplot: interaction counts per condition (Bloc5_04_CellChat_barplot.png)
- Opacity encoding: p < 0.001 = opaque | p < 0.01 = semi | p < 0.05 = light

### Script 04b : CellChat v2 GSE207422, preprint + discussion figures

- Input: Bloc5_03_CellChat_GSE207422_MPR_interactions.csv + NMPR_interactions.csv
- Cells of interest guided by Blocs 3 and 4 CollecTRI findings
- CD8 priority: CD8.TEX, CD8.TPEX
- TAMs priority: IFN-stimulated, SPP1+ immunosuppressive, Lipid-associated, 
  Resident M2 (Stress-response and Monocyte-derived absent from CellChat output)
- Epithelial: Tumor epithelial + Normal epithelial + Ciliated (all kept)
- Preprint figures (CD8 ↔ TAMs - main narrative):
  1. CD8 → TAMs (Bloc5_04b_CellChat_CD8_TAMs_preprint.png)
  2. TAMs → CD8 (Bloc5_04b_CellChat_TAMs_CD8_preprint.png)
- Discussion figures (epithelial axes — GSE207422 only):
  3. CD8 → Epithelial (Bloc5_04b_CellChat_CD8_Epithelial_discussion.png)
  4. Epithelial → CD8 (Bloc5_04b_CellChat_Epithelial_CD8_discussion.png)
  5. TAMs → Epithelial (Bloc5_04b_CellChat_TAMs_Epithelial_discussion.png)
  6. Epithelial → TAM(Bloc5_04b_CellChat_Epithelial_TAMs_discussion.png)
- Output: Results/Figures/BLOC5_Communication/GSE207422/Preprint/

### Script 05 : LIANA+ GSE243013, intercellular communication inference

- Input: Objects/Bloc3_03_seu_Tcells_ProjecTILs.rds + 
         Objects/Bloc4A_04_seu_TAMs_annotated.rds
- Merged TME object: CD8 + TAMs (n=TBD cells)
- Condition column: pathological_response (MPR | non-MPR | pCR)
- Methods: sca + natmi + connectome aggregated via liana_aggregate()
- Resource: Consensus (CellChatDB + OmniPath + others)
- Conversion: SingleCellExperiment (SCE) for Seurat v5 compatibility
- Full object run: MPR + non-MPR + pCR (39,730 interactions)
- Per condition runs: MPR (n=8,936 cells) | non-MPR (n=56,977 cells) | pCR (n=12,118 cells)
- Output: Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_aggregated.csv
          Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_MPR_aggregated.csv
          Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_non_MPR_aggregated.csv
          Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_pCR_aggregated.csv
          Objects/Bloc5_03_seu_TME_GSE243013.rds
          Results/Figures/BLOC5_Communication/GSE243013/Bloc5_03_LIANA_dotplot_CD8_to_TAMs.png
- TAM short labels applied for readability (consistent with Bloc 4A)
- 7 TAM subtypes: Resident M2, Lipid-associated, Inflammatory Mono-derived,
  Stress-response, Proliferating, IFN-stimulated, Classical Mono-derived

### Script 06 : LIANA+ GSE243013, main axis figures

- Input: Bloc5_03_LIANA_GSE243013_MPR_aggregated.csv + non_MPR_aggregated.csv + pCR_aggregated.csv
- Two series of figures:
  Preprint (MPR vs non-MPR):
  1. CD8 → TAMs (Bloc5_04_LIANA_CD8_TAMs_MPR_nonMPR.png)
  2. TAMs → CD8 (Bloc5_04_LIANA_TAMs_CD8_MPR_nonMPR.png)
  Discussion (MPR vs pCR):
  3. CD8 → TAMs (Bloc5_04_LIANA_CD8_TAMs_MPR_pCR.png)
  4. TAMs → CD8 (Bloc5_04_LIANA_TAMs_CD8_MPR_pCR.png)
- Top 5 interactions per target group (aggregate_rank <= 0.05)
- TAM short labels applied for readability

### Script 07 : CellChat v2 GSE243013, differential communication analysis

- Input: Objects/Bloc5_03_seu_TME_GSE243013.rds
- CellChatDB.human full database
- Per condition objects: MPR, non-MPR, pCR
- Parameters: nboot = 100, min.cells = 5
- Short TAM labels applied: "Monocyte FCN1+" harmonized across blocs
- Merged objects: MPR vs non-MPR + MPR vs pCR (separate)
- Output: Objects/Bloc5_07_cellchat_MPR_GSE243013.rds
          Objects/Bloc5_07_cellchat_nonMPR_GSE243013.rds
          Objects/Bloc5_07_cellchat_pCR_GSE243013.rds
          Objects/Bloc5_07_cellchat_merged_MPR_nonMPR_GSE243013.rds
          Objects/Bloc5_07_cellchat_merged_MPR_pCR_GSE243013.rds
          Results/Tables/BLOC5/GSE243013/Bloc5_07_CellChat_GSE243013_MPR_interactions.csv
          Results/Tables/BLOC5/GSE243013/Bloc5_07_CellChat_GSE243013_nonMPR_interactions.csv
          Results/Tables/BLOC5/GSE243013/Bloc5_07_CellChat_GSE243013_pCR_interactions.csv
- Note: MultiNicheNet abandoned — n=3 MPR insufficient for pseudobulk.
- Note: Figures generated from CSV (not netVisual_bubble) to avoid blank page rendering issue with CellChat native figures.

### Script 08b : CellChat v2 GSE243013, preprint + discussion figures

- Input: Results/Tables/BLOC5/GSE243013/Bloc5_07_CellChat_GSE243013_*_interactions.csv
- Cells of interest guided by Blocs 3 and 4 CollecTRI findings
- CD8 priority: CD8.TEX, CD8.TPEX
- TAMs priority: IFN-stimulated, Stress-response, Resident M2, Lipid-associated, Monocyte FCN1+
- Preprint figures (CD8 <-> TAMs MPR vs non-MPR):
  1. CD8 → TAMs (Bloc5_08b_CellChat_CD8_TAMs_MPR_nonMPR_preprint.png)
  2. TAMs → CD8 (Bloc5_08b_CellChat_TAMs_CD8_MPR_nonMPR_preprint.png)
- Discussion figures (MPR vs pCR):
  3. CD8 → TAMs (Bloc5_08b_CellChat_CD8_TAMs_MPR_pCR_preprint.png)
  4. TAMs → CD8 (Bloc5_08b_CellChat_TAMs_CD8_MPR_pCR_preprint.png)
- Output: Results/Figures/BLOC5_Communication/GSE243013/Preprint/
- Note: Script 08 (full figures) skipped — went directly to 08b preprint figures on cells of interest for analytical efficiency.

---

## Preliminary observations : TME composition and CD8 states (Bloc 2-3)

- CD8.TEX gradient: non-MPR > MPR > pCR (visual UMAP), confirmed statistically by ProjecTILs barplot (Chi-2 p<2.2e-16)
- non-MPR dominance (n=42) requires proportional analysis, absolute counts not comparable
- GSE207422 preliminary analysis identified CD8_Exhausted_Terminal enriched in MPR (OR=3.36, different annotation tool). Current harmonized analysis applies identical tools (ProjecTILs, UCell, CollecTRI) to both datasets for direct cross-dataset comparison.

## Preliminary observations : CD8 T cell states (Bloc 3)

- CD8.TEX gradient confirmed: non-MPR (~30%) > MPR (~25%) > pCR (~15%)
- CD8.TPEX enriched in pCR, consistent with reactivation of precursor-exhausted cells
- CD8.EM similar across groups, functional status assessed by UCell (Script 10)
- Proposed gradient: non-MPR (TEX dominant) → MPR (TEX + TPEX) → pCR (TPEX dominant)
Hypothesis: anti-PD-1 reactivated a fraction of TEX toward TPEX in responders

## Preliminary observations : UCell CD8 scoring (bloc 3, GSE243013)

- Exhaustion gradient confirmed: non-MPR (0.347) > pCR (0.249) > MPR (0.220) in CD8.TEX
- CD8.EM exhaustion: non-MPR (0.128) > MPR (0.107) ≈ pCR (0.106), EM dysfunction in non-responders
- CD8.TPEX retains higher cytotoxicity than CD8.TEX across all groups,preserved precursor function
- Memory score highest in CD8.NaiveLike = expected

**Hypotheses:**
- H1 (CollecTRI): superior intrinsic CD8 reactivation capacity in responders
- H2 (CellChat/NicheNet): immunosuppressive TME determines reactivation failure
- H3: interaction between CD8 intrinsic state and TME context

## Preliminary observations : TCR analysis (Bloc 3, GSE243013)

- CD8.TEX 85.4% expanded, confirms tumor-reactive identity
- Clonal diversity: MPR (0.490) > pCR (0.458) > non-MPR (0.438), responders maintain broader TCR repertoire
- pCR more polyclonal, broader anti-tumor response

## Preliminary observations : GSE207422 CD8 cross-validation (Bloc 3)

**ProjecTILs:**
- 7 CD8 states identified, consistent with GSE243013
- CD8.TEX dominant MPR ~55% ≈ NMPR ~54% , no gradient (absent pCR group, n=13)
- CD8.TPEX slightly enriched MPR, CD8.TEMRA enriched NMPR

**UCell scoring:**
- Exhaustion NMPR > MPR, Cytotoxicity MPR > NMPR 
- TPEX_UCell MPR > NMPR, Memory_UCell NMPR > MPR (anergy interpretation)
- CD8.TEX MPR : simultaneous high exhaustion + cytotoxicity, intra-TEX heterogeneity confirmed 

**Cross-dataset consensus:**
- TPEX more active in responders, TEMRA enriched non-responders
- Intra-TEX heterogeneity confirmed both datasets
- Memory_UCell higher non-responders, consistent with anergy/quiescence

## Preliminary observations : CollecTRI CD8 (Bloc 3)

**Methodological correction:** Initial violin plots used manually selected TFs based on prior DoRothEA hypotheses-confirmation bias risk. Corrected to Top 6 by cross-cell variance (objective). Top 10 and Top 15 tested, conclusions unchanged.

**Visualization approach:** Heatmaps row-scaled (z-score) show state-specific TF patterns within each CD8 state per condition. Violin plots (Top 6 by variance, absolute scores) show global MPR vs non-MPR differences. The two are complementary.

**Two complementary visualizations:**
- **Violin plots (Top 6 by variance, absolute scores)** , global MPR vs non-MPR comparison across CD8.TEX and CD8.TPEX
- **Heatmaps (scale=row, Top 20 by variance)** , state-specific TF patterns within each CD8 state per condition

**Violin plots, cross-dataset consensus (objective Top 6):**
- ELK4 enriched MPR CD8.TEX = cross-dataset
- MYC/HSF1/HSF2 enriched non-responders = cross-dataset
- STAT2/TBX21 : not in objective Top 6 : DoRothEA preliminary only

**Heatmap observations , state-specific patterns:**

GSE243013 :
- CD8.EM non-MPR : MHC II program dominant (CIITA, RFXAP, RFX5, RFXANK)
- CD8.TEX non-MPR : TBX21 dominant WITHOUT ELK4 → abortive cytotoxic program
- CD8.TPEX non-MPR : MYC, NFKB1, RELA, HSF1/2 → chronic proliferation without cytotoxic engagement
- CD8.TEX MPR : MHC II (RFXAP/RFXANK/RFX5) + TBX21 + ELK4 → functional cytotoxic program
- CD8.TPEX MPR : ELK4 dominant, TBX21 quasi-absent → pure reactivation program
- CD8.TEX pCR : HOPX dominant, ELK4 absent → post-response quiescence
- CD8.EM pCR : most transcriptionally diverse state, nearly all 20 TFs active → polyfonctional memory surveillance

GSE207422 :
- CD8.TEX NMPR : MHC II (RFXAP/RFXANK/RFX5/CIITA) + STAT1 + ZBTB4
- CD8.TPEX NMPR : MYC, E2F1, E2F4, SRSF2, TFDP1 → chronic proliferation = consistent with GSE243013
- CD8.TEX MPR : ELK4 dominant + NFYC + ABL1 + DOT1L
- CD8.TPEX MPR : ELK4 + HSF1 + ABL1 + NFYC → ELK4-driven reactivation = consistent

**Cross-dataset conclusion (violin + heatmap combined):**
- ELK4 most robust discriminating TF in MPR , dominant in TPEX MPR both datasets, present in TEX MPR GSE207422
- TBX21 in non-MPR CD8.TEX WITHOUT ELK4 = abortive cytotoxic program; TBX21 in MPR CD8.TEX WITH ELK4 = functional cytotoxic program
- ABL1 and NFYC as consistent ELK4 co-activators in MPR TEX/TPEX (GSE207422) and MPR TPEX (GSE243013)
- non-MPR TPEX : chronic proliferation program (MYC/E2F/HSF) without cytotoxic coordination : confirmed cross-dataset
- pCR CD8.TEX : HOPX-dominant quiescence : post-complete response state

## Preliminary observations : Pseudobulk CD8 MPR vs pCR (Bloc3 Script 13)

All UCell signatures showed no significant difference between MPR and pCR at patient level (Cytotoxicity p=0.5573, Exhaustion p=0.5116, Memory p=0.1734, TPEX p=0.6047). Inter-patient heterogeneity appeared similar between the two groups, contrasting with TAMs where pCR showed greater heterogeneity than MPR (Bloc4A Script 06). This pattern may suggest biological similarity between MPR 
and pCR at the CD8 level, though this interpretation requires caution given the small sample size (n=10 vs n=11).

**Perspectives:**
- TAM-CD8 interaction analysis (CellChat/NicheNet) may help resolve whether the immunosuppressive TAM context differentially constrains CD8 function in MPR vs pCR; planned in subsequent blocs
- Larger balanced cohort with sufficient statistical power would be required to formally test CD8-level differences between MPR and pCR at patient level
- Spatial transcriptomics could resolve TAM/CD8 co-localization patterns in the TME to complement these transcriptomic observations

## Preliminary observations : Pseudobulk CD8 MPR vs pCR (Bloc 3, GSE243013)

- All UCell signatures ns (Cytotoxicity p=0.558, Exhaustion p=0.512, Memory p=0.173, TPEX p=0.605)
- Inter-patient heterogeneity similar between MPR and pCR  may reflect biological similarity at CD8 level (n=10 vs n=11, interpret with caution)
- MPR vs pCR distinction may be driven by TAM compartment rather than CD8 intrinsic state, to be explored via CellChat/NicheNet (Bloc 5)

## Preliminary observations : TAMs UCell (Bloc 4, GSE207422)

**UCell scores (9 subtypes combined):**
- M2_immunosuppressive: NMPR > MPR, p<0.0001 
- SPP1_signature: NMPR > MPR, p<0.0001 
- IFN_response: MPR > NMPR, p=3e-04 
- M1_inflammatory: ns (p=0.2848) limited power (MPR n=3)

**Key subtype enrichments:**
- NMPR: TAM_like_lipid, monocyte, regulatory, SPP1, stress-immunosuppressive niche
- MPR: TAM_like_IFN, dual role hypothesis (CXCL9/10 pro-immunogenic vs suppressive)
- TAM_like_resident_M2 MPR > NMPR, discordant with GSE243013, may reflect NSCLC vs LUAD composition differences

**Biological interpretation:**
- NMPR immunosuppressive TAM niche confirmed, supports H2
- MPR as biologically unstable intermediate, two divergent trajectories:
  1. Relapse: residual immunosuppressive TAMs + progressive CD8 exhaustion
  2. Deepening toward pCR: pro-immunogenic TAMs + CD8 TPEX plasticity

## Preliminary observations : TAMs UCell (Bloc 4, GSE243013)

**UCell scores pairwise Wilcoxon Bonferroni (p<0.0167):**
- M2_immunosuppressive: non-MPR > pCR > MPR, all pairs p<0.005 
- SPP1_signature: non-MPR > MPR > pCR , all pairs p<0.0001 
- IFN_response: non-MPR dominant; MPR ≈ pCR (p=0.783 ns) 
- M1_inflammatory: non-MPR > MPR > pCR , MPR/pCR ns Bonferroni

**Key observations:**
- Gradient non-MPR > MPR > pCR confirmed all signatures , MPR as biologically unstable intermediate state 
- IFN_response exception: discriminates responders from non-responders, not MPR from pCR
- IFN-stimulated TAMs: hybrid M1/immunosuppressive profile, TME biology, 
not artifact

**Biological interpretation:**
- Non-MPR: triple immunosuppression (M2 + SPP1 + IFN TAMs)
- pCR more heterogeneous than MPR at patient level, multiple pathways to complete response

## Preliminary observations : Pseudobulk UCell_TAMs MPR vs pCR (Bloc 4, GSE243013)

- All signatures ns (IFN p=0.512, M1 p=0.152, M2 p=0.654, SPP1 p=0.468)
- pCR more heterogeneous than MPR at patient level, Case 1 (insufficient power) vs CD8 pseudobulk where heterogeneity similar, Case 2 (biological similarity)
- TAM signatures alone insufficient to separate MPR from pCR
- Integration of CD8 + TAM signatures required planned via CellChat/NicheNet

## Preliminary observations : CollecTRI TAMs (Bloc 4)

**Two complementary visualizations:**
- **Violin plots (Top 6 by variance, absolute scores)** : global MPR vs non-MPR comparison
- **Heatmaps (scale=row, Top 20 by variance)** , subtype-specific TF patterns per condition

**Violin plots , cross-dataset consensus (objective Top 6):**
- MHC II program (RFXAP/RFXANK/RFX5/CIITA) enriched responders both datasets 
- HSF1 enriched non-responders both datasets : chronic stress program
- RELA/REL similar between groups — NMPR wider distribution (inter-patient heterogeneity)

**Heatmap observations , subtype-specific patterns:**

GSE243013 :
- Resident M2 MPR : ELK4 dominant : pro-immunogenic signal 
- Resident M2 non-MPR : HIF1A, NFKB : hypoxic/inflammatory stress
- Resident M2 pCR : ELK4 + DMTF1
- Monocyte FCN1+ non-MPR : HIF1A dominant : hypoxic adaptation
- Stress-response non-MPR : MYC + HSF1/HSF2 : chronic proliferation 
- Proliferating pCR : MYC very strong
- Classical-Mono MPR : ELK4 + MHC II dominant 
- Classical-Mono pCR : ELK4 + MHC II present but weaker
- IFN-stimulated : STAT1 constitutive all conditions, not discriminating

GSE207422 :
- IFN-stimulated MPR : IRF1, STAT1, IRF5, REL, strong IFN program
- MRC1+ M2-like MPR : MHC II dominant , consistent with GSE243013
- SPP1+ immunosuppressive NMPR : highly active, heterogeneous : complex immunosuppressive program
- Monocyte-derived NMPR : HIF1A dominant , consistent with GSE243013
- Stress-response : active both conditions, constitutive

**Cross-dataset consensus (violin + heatmap combined):**
- MHC II enriched MPR TAMs : Resident M2 GSE243013, MRC1+ M2-like GSE207422 
- ELK4 enriched MPR Resident M2 and Classical-Mono GSE243013 : absent GSE207422
  (possible LUAD vs NSCLC histological difference)
- HIF1A enriched non-responder Monocyte TAMs both datasets 
- MYC/HSF enriched non-responder Stress-response TAMs GSE243013 
- SPP1+ immunosuppressive NMPR most transcriptionally active GSE207422 
- IFN-stimulated STAT1 constitutive both datasets : not discriminating

## Preliminary observations : Epithelial compartment CNV (Bloc 4B, GSE207422)

- Epithelial composition: Basal dominant in MPR (~62%), EMT and Tumor_epithelial exclusive to NMPR, transcriptional rather than genomic drivers of non-response
- CNV clonal profiles (SCEVAN) identical between MPR and NMPR: treatment response not determined by chromosomal architecture differences
- Both SCEVAN and CopyKAT confirm higher tumor/aneuploid fraction in NMPR (~47%/~75%) vs MPR (~32%/~38%) : directional concordance between tools
- EMT cells: maximal discordance between tools (SCEVAN 0% tumor vs CopyKAT 38% aneuploid); whether EMT resistance is driven by subtle genomic alterations or transcriptional reprogramming alone remains unresolved
- Tumor_epithelial subtype: strong cross-tool concordance (99% tumor/aneuploid both tools)
- Next steps: CytoTRACE, UCell, CollecTRI on confirmed tumor vs normal epithelial cells per condition to characterize transcriptional state differences driving non-response

## Preliminary observations : CytoTRACE2 epithelial differentiation (Bloc 4B, GSE207422 — SCEVAN labels)

**Violin plot observations (CytoTRACE2 Score, 0=differentiated, 1=totipotent) :**

- Ciliated : very low scores (~0.05-0.10), terminally differentiated, expected for mature bronchial cells
- EMT_NMPR : low scores (~0.10-0.25), majority converging toward differentiation, supports SCEVAN classification (0% tumor); EMT reprogramming appears terminal rather than stem-like
- normal_MPR : pear-shaped distribution with two subpopulations coexist; one highly differentiated and one less differentiated, suggesting epithelial plasticity in the normal compartment under MPR conditions
- normal_NMPR : mass concentrated at low scores, uniformly differentiated, transcriptionally static normal epithelial cells in non-responding patients
- tumor_MPR : compact distribution, high scores (~0.40-0.65), homogeneously stem-like, low differentiation, consistent with a targetable tumor population
- tumor_NMPR : broad bimodal distribution (~0.15-0.80), high intra-tumoral differentiation heterogeneity, consistent with treatment resistance mechanisms

**Biological interpretation :**

In MPR, a dissociation is observed between tumor and normal compartments : tumor cells are homogeneously stem-like while normal epithelial cells display plasticity, suggesting active epithelial remodeling in the treatment-responsive microenvironment. In NMPR, tumor cells show greater heterogeneity while normal cells are uniformly static, a pattern potentially reflecting impaired epithelial dynamics contributing to non-response. EMT_NMPR low differentiation scores support SCEVAN over CopyKAT for this subtype classification, pending CopyKAT-based CytoTRACE2 cross-validation.

**EMT reclassification (Script 09c):**
EMT_NMPR cells were reclassified as normal_NMPR based on three converging lines of evidence:
1. SCEVAN highest specificity (0.75) in published benchmarks vs CopyKAT tendency to overestimate tumor fractions (Lanucara et al., Biomedicines 2024)
2. CytoTRACE2 low stemness scores for EMT cells under both SCEVAN and CopyKAT labels, consistent with terminally differentiated rather than stem-like state
3. CopyKAT cross-validation confirming identical low-score distribution for EMT_NMPR

**Final groups for downstream analyses (UCell, CollecTRI):**
tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, Ciliated

**Literature support for normal epithelial compartment hypothesis:**
- Hu et al. 2023 (Genome Medicine, GSE207422): normal epithelial cells expand in MPR TME after neoadjuvant PD-1 + chemotherapy
- Cui et al. 2025 (Molecular Cancer): neoadjuvant chemoimmunotherapy induces immunosuppressive microenvironment in normal epithelial cells
- Supports CellChat/NicheNet hypothesis: normal epithelial plasticity in MPR may associate with pro-immunogenic signaling to CD8 and TAMs

## Preliminary observations : UCell epithelial compartment (Bloc 4B, GSE207422)

**Signatures applied (13 total):**
- MSigDB Hallmark (Liberzon et al., Cell Systems 2015): Proliferation_E2F, Proliferation_G2M, Apoptosis, EMT, IFN_gamma_response, IL6_JAK_STAT3, TNFA_NFkB, WNT_beta_catenin, Notch_signaling, Unfolded_protein, Hypoxia
- Custom: HSF1_targets (Mendillo et al., Cell 2012), Antigen_presentation (Hu et al., Genome Medicine 2023)

**Dotplot observations (ranked by expression per group):**

tumor_NMPR: Proliferation_E2F > Proliferation_G2M > Apoptosis > TNFA_NFkB > Unfolded_protein > HSF1_targets > Hypoxia = Aggressive survival program: massive proliferation, chronic proteotoxic stress, inflammatory resistance, absent antigen presentation

tumor_MPR: IL6_JAK_STAT3 > WNT_beta_catenin > Notch_signaling > Hypoxia > Unfolded_protein > Proliferation_E2F > Proliferation_G2M = Plasticity program: developmental signaling dominant, moderate proliferation, communicative tumor phenotype potentially sensitive to immunotherapy

normal_NMPR: EMT > IFN_gamma_response > Hypoxia > TNFA_NFkB > Apoptosis = Dysfunctional inflammatory compartment: EMT residual signature (reclassified cells) + sterile inflammation amplified by IFN without coordinated antigen presentation

normal_MPR: Antigen_presentation > IFN_gamma_response > IL6_JAK_STAT3 = Pro-immunogenic compartment: normal epithelial cells actively present antigens and respond to IFN in a coordinated manner, supporting immune activation

Ciliated: Antigen_presentation > Notch_signaling = Stable reference: constitutive antigen presentation, identity maintenance

**Key biological interpretations:**

- IFN_gamma_response higher in normal_NMPR than normal_MPR but coupled to EMT and TNFA_NFkB rather than antigen presentation — IFN acts as a contextual amplifier: pro-immunogenic in MPR, pro-inflammatory/immunosuppressive in NMPR
- HSF1_targets confirmed in tumor_NMPR (compact, high, uniform distribution) > tumor_MPR : completes the HSF1 feedback loop across three compartments: tumor_NMPR → TAMs → CD8 exhaustion
- Unfolded_protein and HSF1_targets concordant in tumor_NMPR, independent signatures converging on same chronic proteotoxic stress program (internal validation)
- Antigen_presentation quasi-absent in tumor_NMPR → complete immune evasion through loss of antigen presentation
- Two distinct tumor survival programs identified:
  tumor_MPR: IL6/JAK/STAT3 + WNT + Notch (plasticity, communicative)
  tumor_NMPR: HSF1 + Proliferation + TNFA (stress resistance, autonomous)
- Apoptosis slightly higher in NMPR groups, tumor cells resist death signals (BIRC5/BCL2 program) while some normal NMPR cells succumb, interpretation limited by absence of pCR group in GSE207422
- Hypoxia higher in tumor_MPR than tumor_NMPR, consistent with active immune remodeling increasing local oxygen consumption rather than hypoxic escape

**Cross-compartment HSF1 conclusion:**
HSF1_targets enriched in tumor_NMPR (epithelial) + non-responder TAMs (Bloc 4A) + non-responder CD8 (Bloc 3), suggesting a HSF1-driven immunosuppressive feedback loop originating from the tumor epithelial compartment propagating through the TME.

## Preliminary observations : CollecTRI TF activity epithelial compartment (Bloc 4B, Script 11)

**Top 20 TFs by cross-cell variance:**
MYC, RFXAP, RFXANK, HSF2, RFX5, CIITA, NKX2-1, SP1, FOXJ1, DMTF1, NFE2L2, E2F4,
HSF1, ELK4, TP53, DACH1, HIF1A, E2F1, NFKB, JUN

**TF programs per group (heatmap + violin observations) :**

tumor_MPR: TP53, MYC, NFE2L2, HIF1A, SP1, NFKB, ELK4, JUN, E2F4, CIITA = Mixed plasticity/stress program with pro-immunogenic window (ELK4 + CIITA)

tumor_NMPR: E2F4, HSF2, E2F1, HSF1, DMTF1, TP53, ELK4, SP1, JUN, NFE2L2, HIF1A, NFKB, MYC = Proliferation (E2F4/E2F1) + chronic proteotoxic stress (HSF1/HSF2) , HSF1 confirmed 

normal_MPR: NKX2-1, RFXAP, RFX5, CIITA, RFXANK, DACH1 = MHC II dominant , pro-immunogenic compartment confirmed

normal_NMPR: DACH1, HIF1A, NFKB, JUN, NKX2-1, signal globally weak = Transcriptionally quiescent, uncoordinated inflammatory state

Ciliated: FOXJ1, RFXANK, RFXAP, RFX5, DACH1, CIITA = Ciliary identity (FOXJ1) + constitutive MHC II : stable reference

**Key cross-compartiment findings:**

**HSF1 feedback loop confirmed (3rd compartment):**
HSF1 enriched in tumor_NMPR epithelial (CollecTRI + UCell) : ORIGIN
+ non-responder TAMs (both datasets) + non-responder CD8 (both datasets)
= HSF1-driven immunosuppressive cascade originating in tumor epithelial compartment

**ELK4 contextual modulator hypothesis confirmed:**
- tumor_MPR: ELK4 + CIITA/TP53 = pro-immunogenic context
- tumor_NMPR: ELK4 + HSF1/E2F = immunosuppressive context
- CD8 MPR: ELK4 + TBX21 = functional cytotoxic program
- CD8 NMPR: ELK4 without TBX21 = abortive cytotoxic program
- TAMs MPR: ELK4 Resident M2 + Classical-Mono (GSE243013)
   - ELK4 functional output determined by co-regulatory partners, not expression level
   - Supporting references: Yao 2013 (BMC Genomics), Xue et al. 2023 (Advanced Science)

**MHC II program cross-compartment (MPR):**
RFXAP/RFXANK/RFX5/CIITA enriched in:
- normal_MPR epithelial (dominant)
- tumor_MPR epithelial (CIITA present)
- Ciliated (constitutive)
- TAMs MPR (Resident M2, Classical-Mono, MRC1+) = Global pro-immunogenic MHC II program in MPR TME

**Literature support:**
- Lanucara et al., Biomedicines 2024, SCEVAN/CopyKAT benchmark
- Hu et al., Genome Medicine 2023, normal epithelial expansion in MPR
- Cui et al., Molecular Cancer 2025, normal epithelial microenvironment after immunotherapy
- Yao, BMC Genomics 2013, ELK4 contextual modulator in macrophages
- Xue et al., Advanced Science 2023, ELK4 co-factor dependency in cancer

### Preliminary observations for LIANA+ GSE207422

**CD8 to TAMs :**
MPR: CCL5-CCR1 (recruitment, 4 subtypes), TNFSF9-HLA-DPA1 (MHC II co-stimulation), CCL5-SDC4 (resolutive remodeling), CD52-SIGLEC10 (IFN-stimulated only). Monocyte-derived: no interaction in MPR.
NMPR: CCL5-CCRL2 (decoy diversion cross-dataset), GZMB-IGF2R (misdirected cytotoxicity), HLA-LILRB1/2 + B2M-LILRB2 (double trap Monocyte-derived), ENTPD1-TMIGD3 (adenosine pathway SPP1+).

**TAMs to CD8 :**
MPR: CLEC4G-LAG3 + HLA-DRB3-LAG3 (modulatory LAG3 TPEX+TEX), GRN-NTRK1 (TPEX survival), HLA-CD8A (antigen presentation). Zero shared interactions with NMPR.
NMPR: HLA-DQA1/DQB1-LAG3 (canonical inhibitory checkpoint TPEX+TEX), S100A8-CD69 (TPEX retention), CCL13-CXCR3 (TEX hostile recruitment).

Bidirectional CD8-TAMs: MPR = pro-immunogenic self-reinforcing loop. 
NMPR = immunosuppressive self-reinforcing loop. Key: LAG3 via CLEC4G (modulatory) in MPR vs HLA-DQA1/DQB1 (canonical inhibitory) in NMPR.

**CD8 to Epithelial :**
Shared: CCL5-SDC1/SDC4, XCL1-ADGRV1.
MPR: CRTAM-CADM1 (functional immune-tumor contact), HMGB1-SDC1, CXCL13-GRM7 (TLS).
NMPR: COL6A3-ITGA/ITGB (fibrosis expanded to Normal+Ciliated), GZMA-PARD3 (misdirected cytotoxicity), HSP90AA1-EGFR (HSF1 target toward Normal epithelial).

**Epithelial to CD8 :**
Shared: CD59-CD2, NECTIN4-TIGIT (TEX both, expanded to TPEX in NMPR).
MPR: GPI/COPA-NTRK1 (TPEX survival), HLA-A-CD8A (TEX presentation).
NMPR: NECTIN4-TIGIT expanded to TPEX (locks reactivable precursors), COL21A1/COL28A1-ITGB1 (fibrosis toward TEX).

Bidirectional CD8-Epithelial: MPR = functional dialogue. NMPR = mutual exclusion, each compartment amplifies the other's dysfunction.

**TAMs to Epithelial :**
Shared: EREG-ERBB2/ERBB4 (Ciliated only).
MPR: HGF-SDC1/MET (regeneration, strongest resolutive signal), EREG-EGFR/ERBB, LAMA4-DAG1/ITGA6, SPINK1-EGFR.
NMPR: CXCL8-SDC1 (HSF1 in Tumor), IL6-F3 (HSF1 in Normal), FN1-SDC1 (fibrosis), F13A1-ITGB1 (coagulation), VCAN-EGFR (pro-tumoral).

**Epithelial to TAMs :**
Shared: SLPI-CD4/PLSCR1 (Lipid-associated + IFN-stimulated).
MPR: MEGF10-ABCA1 (phagocytosis), CD59-STAB1 (complement), GPC3-CD81. Monocyte-derived: no interaction, convergence with CD8 to TAMs MPR.
NMPR: SAA1-FPR1/FPR2/CD36 (alarmin activates Monocyte-derived), SCGB3A1-MARCO (macrophage recruitment), EFNB3-EPHB2 (pro-invasive, 3 subtypes).

Bidirectional TAMs-Epithelial: MPR = regenerative TME (HGF dominant). NMPR = hostile self-amplifying network (HSF1 targets + fibrosis + SAA1 alarmin).
Key: Monocyte-derived ignored in MPR in BOTH CD8→TAMs AND Epithelial→TAMs 

**HSF1 cross-compartment network:**
CXCL8 and IL6 (HSF1 direct targets CollecTRI) in TAMs→Epithelial NMPR. 
HSP90AA1 (HSF1 target) from CD8.TEX toward Normal epithelial and TAM_like_stress. 
Cross-compartment propagation across CD8, TAMs, and epithelial simultaneously,consistent with self-reinforcing network rather than linear cascade. Spatial transcriptomics required to resolve directionality.

### Preliminary observations for CellChat GSE207422

**CD8 → TAMs (MPR vs NMPR) :**

Shared interactions: CCL5-CCR1 (Lipid-associated + IFN-stimulated both conditions),
MT-RNR2-FPRL2 (all subtypes both conditions), PTPRC-MRC1 (both conditions)

MPR dominant:
- HLA-DRB1-CD4 → Resident M2 + SPP1+ + IFN-stimulated : CD4 recruitment signal via MHC II — MPR exclusive. Resident M2 receives double signal HLA-DRA + HLA-DRB1-CD4 = strongest CD4 recruitment target in MPR
- CCL5-CCR1 : higher Comm.Prob. in MPR

NMPR: same interactions, globally attenuated. No qualitative exclusive interaction.
MT-RNR2-FPRL2 and CCL5-CCR1 stronger in NMPR quantitatively.

Conclusion: mainly quantitative difference + one qualitative signal = HLA-DRB1-CD4 exclusive MPR = CD4 helper recruitment toward pro-immunogenic TAMs.

**TAMs → CD8 (MPR vs NMPR) :**

Both conditions: HLA-A-CD8A + HLA-C-CD8A → TPEX + TEX
MPR: higher Comm.Prob. on HLA-A/C-CD8A
NMPR: same interactions, attenuated

Conclusion: purely quantitative difference. MHC I antigen presentation conserved both conditions, stronger in MPR. No qualitative exclusive interaction detected by CellChat on this axis. Biological richness (LAG3 signals) from LIANA+ only.

**CD8 → Epithelial (MPR vs NMPR) :**

MPR dominant:
- CRTAM-CADM1 → Normal epithelial : functional immune recognition = cross-tool confirmed
- CD96-NECTIN1, CD99-CD99 → Tumor epithelial
- HLA-DRA/DRB1-CD4 → Ciliated : CD4 recruitment via 2 alleles
- PPIA-BSG → Tumor + Normal epithelial (both conditions, stronger MPR)
- SEMA4D-PLXNB2 → Normal epithelial

NMPR dominant:
- GZMA-PARD3 → Tumor + Normal epithelial = cross-tool confirmed
- PPIA-BSG amplified → Tumor + Normal epithelial
- HLA-DPA1/DPB1/DQA1/DRA/DRB1-CD4 → Ciliated : pan-HLA CD4 recruitment attempt (5 alleles vs 2 in MPR), unsuccessful compensation
- IFNG-IFNGR1/2 → Normal epithelial : displaced inflammatory signal
- SIRPG-CD47 → Normal epithelial

Conclusion: MPR = functional immune-tumor contact (CRTAM-CADM1) + targeted CD4 recruitment (2 HLA alleles). NMPR = misdirected cytotoxicity (GZMA-PARD3) + pan-HLA compensation attempt (5 alleles) disconnected from tumor compartment.

**Epithelial → CD8 (MPR vs NMPR) :**

Both conditions: HLA-A-CD8A → TEX (stronger MPR) = cross-tool
NMPR exclusive: HLA-B-CD8A + HLA-C-CD8A → TPEX + TEX = pan-HLA engagement in NMPR, attenuated vs MPR

Conclusion: quantitative difference , epithelium maintains MHC I presentation in both conditions but attenuated in NMPR, compensated by HLA-B allele addition.

**TAMs → Epithelial (MPR vs NMPR) :**

MPR dominant:
- LGALS9-CD44 → Tumor + Normal : adhesion/plasticity
- HLA-DPA1/DRA-CD4 → Ciliated : CD4 recruitment from TAMs
- APP-CD74 → Normal : MHC II/CLIP pro-immunogenic
- PPIA-BSG → Tumor + Normal (both conditions)
- LGALS9-P4HB → Tumor + Normal (both conditions)

NMPR dominant:
- FN1-SDC1 → Tumor + Normal epithelial, cross-tool (Tumor only confirmed)
- FN1-SDC4 → Ciliated : fibrotic signal extended to airway
- SPP1-CD44 → Tumor epithelial : pro-tumoral osteopontin
- PPIA-BSG amplified → Tumor + Normal
- LGALS9-P4HB → stress signal amplified

Conclusion: MPR = resolutive remodeling + CD4 recruitment + antigen presentation. 
NMPR = active fibrosis (FN1-SDC1) + pro-tumoral (SPP1-CD44) + stress amplified.

**Epithelial → TAMs (MPR vs NMPR) :**

Both conditions: APP-CD74 + MT-RNR2-FPRL2 → all 4 TAM subtypes
MPR exclusive: HLA-DRA-CD4 → SPP1+ + IFN-stimulated : epithelium sends additional CD4 recruitment signal toward pro-immunogenic TAMs in MPR
Quantitative: APP-CD74 stronger MPR, MT-RNR2-FPRL2 stronger NMPR

Conclusion: shared program + MPR adds HLA-DRA-CD4 signal. MT-RNR2-FPRL2 amplified in NMPR = mitochondrial stress signal from epithelium toward TAMs.

**Cross-tool convergence LIANA+ / CellChat :**

Confirmed by both tools :
- CCL5-CCR1 → Lipid-associated + IFN-stimulated (CD8→TAMs, both conditions)
- HLA-A/C-CD8A → TPEX + TEX (TAMs→CD8, both conditions)
- GZMA-PARD3 → Normal + Tumor epithelial NMPR (CD8→Epithelial)
- CRTAM-CADM1 → Normal epithelial MPR (CD8→Epithelial)
- HLA-A-CD8A → TEX MPR (Epithelial→CD8)
- FN1-SDC1 → Tumor epithelial NMPR (TAMs→Epithelial)

Tool-specific, no cross-tool convergence:
- Epithelial→TAMs axis : complete divergence between tools
- LAG3 signals (TAMs→CD8) : LIANA+ only
- HGF-SDC1, CXCL8-SDC1, IL6-F3 : LIANA+ only
- SPP1-CD44, PPIA-BSG, LGALS9 : CellChat only

**HSF1 axis :**
CXCL8 and IL6 (direct HSF1 targets CollecTRI) detected in NMPR by LIANA+ only (TAMs→Epithelial). HSP90AA1 detected by LIANA+ only (CD8→Epithelial + supplementary figure). Not confirmed by CellChat , database resource difference.SPP1 and PPIA not registered as direct HSF1 targets in CollecTRI.

**CD4 recruitment perspective :**
HLA/MHC II signals toward CD4 detected across three compartments in NMPR (CD8→Epithelial Ciliated, TAMs→Epithelial Ciliated, Epithelial→TAMs) ,consistent across both tools for HLA-DRA-CD4 signal. Pattern suggests unsuccessful CD4 helper recruitment in non-responders, disconnected from 
tumor compartment (Kagamu et al., Nat Commun 2026; Sade-Feldman et al., Cell 2018).

### Preliminary observations for LIANA+ GSE243013

**CD8 → TAMs (MPR vs non-MPR) :**

Shared interactions: CD99-PILRA (Stress-response + Monocyte FCN1+), CIRBP-TREM1 (Stress-response + Monocyte FCN1+), HMGB1-THBD (different 
targets per condition), RPS19-C5AR1 (different targets per condition),TNFSF9-HLA-DPA1 (present both conditions, expanded in non-MPR)

MPR dominant interactions:
- CCL5-ACKR1 + CXCL13-ACKR1 toward Stress-response : TLS chemokine signal = confirmed cross-dataset
- HMGB1-CD163 toward Monocyte FCN1+ + Lipid-associated : resolutive remodeling
- ANXA1-FPR1 toward IFN-stimulated : anti-inflammatory resolution signal
- B2M-LILRB2 + HLA-A/C-LILRB2 toward IFN-stimulated : inhibitory checkpoint paradoxically in MPR

non-MPR dominant interactions:
- CCL5-CCRL2 toward Stress-response + Monocyte FCN1+ : decoy receptor diversion = confirmed cross-dataset
- HMGB1-TLR2 toward Resident M2 : danger/pro-inflammatory signal
- IFNG-IFNGR1/2 toward Resident M2 + Lipid-associated : diffuse IFN, non-coordinated
- HLA-A/C-LILRB2 toward IFN-stimulated : inhibitory MHC checkpoint
- TIGIT-NECTIN2 toward IFN-stimulated : inhibitory checkpoint — double locking of IFN-stimulated TAMs (HLA-LILRB2 + TIGIT-NECTIN2) = confirmed cross-dataset
- TNFSF9-HLA-DPA1 expanded to Lipid-associated + IFN-stimulated : paradoxically amplified but uncoordinated

Conclusion: MPR = TLS-oriented coordinated program (CCL5/CXCL13-ACKR1) + resolutive signals. non-MPR = CCL5 decoy diversion + double inhibitory checkpoint on IFN-stimulated (most pro-immunogenic TAM subtype specifically locked) + diffuse IFNG. Cross-dataset confirmation.

**TAMs → CD8 (MPR vs non-MPR) :**

Shared interactions: CXCL10-CXCR3 + CXCL11-CXCR3 toward TPEX + TEX (both conditions), HLA-DRA-LAG3 toward TPEX + TEX (both conditions)

MPR dominant interactions:
- CXCL9-CXCR3 toward TEX : additional chemokine recruitment signal = MPR exclusive

non-MPR dominant interactions:
- C3-IFITM1 toward TPEX : complement/antiviral signal
- SECTM1-CD7 toward TPEX : co-stimulation signal
- CXCL11-CCR5 toward TEX : retention signal via CCR5 (vs CXCR3 = recruitment), same ligand CXCL11 activates both recruitment (CXCR3) and retention (CCR5) simultaneously in non-MPR

Bidirectional loop conclusion CD8-TAMs:
non-MPR TAMs deploy a dual-signal trapping mechanism toward CD8 TEX, chemokine attraction (CXCL9/10/11-CXCR3) paired with checkpoint activation (HLA-DRA-LAG3) and retention (CXCL11-CCR5). CD8 recruitment itself is co-opted as an immunosuppressive mechanism in non-responders.
MPR = coordinated pro-immunogenic loop with additional CXCL9-CXCR3 signal.

**pCR observations (MPR vs pCR) :**
Note: interpret with caution, single dataset, pseudobulk MPR/pCR non-significant in Blocs 3-4, low pCR n.

CD8 → TAMs:
- Shared: CIRBP-TREM1, RPS19-C5AR1, HMGB1-CD163, TNFSF9-HLA-DPA1 → Resident M2
- pCR shares inhibitory signals with non-MPR: HLA-A/C-LILRB2, TIGIT-NECTIN2, HMGB1-TLR2 toward IFN-stimulated/Monocyte FCN1+
- MPR exclusive vs pCR: CCL5/CXCL13-ACKR1 (TLS signal lost in pCR), ANXA1-FPR1, HMGB1-CD163 toward Lipid-associated
- pCR specific: COPA-P2RY6, HSPA1A-TLR4, TIGIT-PVR (weak)

TAMs → CD8:
- Shared: CXCL9/10/11-CXCR3 toward TPEX + TEX
- MPR exclusive: HLA-DRA-LAG3 → TPEX + TEX — absent pCR. Exploratory 
  hypothesis: residual LAG3 signaling may limit response magnitude, preventing complete tumor elimination
- pCR specific: CXCL11-CCR5 (retention, same as non-MPR), NECTIN2-CD226 (activating — DNAM-1), NECTIN2-TIGIT (inhibitory), same NECTIN2 ligand activates opposing receptors simultaneously in pCR

Conclusion: pCR is biologically distinct from both MPR and non-MPR. 
Shares inhibitory signals with non-MPR but without the decoy diversion (no CCL5-CCRL2). TLS signals lost vs MPR. NECTIN2 contradictory signal 
(CD226 activating + TIGIT inhibitory) unique to pCR. Complete pathological response does not reflect fully resolved immune TME.

### Preliminary observations for CellChat GSE243013** 

**CD8 → TAMs (MPR vs non-MPR) :**

Shared interactions: MIF-(CD74+CD44) toward Stress-response, Resident M2, IFN-stimulated (both conditions, similar intensity); PTPRC-MRC1 toward Stress-response, Resident M2, Lipid-associated (both conditions); PPIA-BSG toward Lipid-associated (both conditions); MIF-(CD74+CXCR4) toward Resident M2 (both conditions)

MPR dominant interactions:
- ANXA1-FPR1 toward IFN-stimulated : resolution signal = cross-tool (LIANA+)
- CD99-PILRA toward IFN-stimulated : immune recognition = cross-tool (LIANA+)
- HLA-DRB1-CD4 toward Stress-response : CD4 recruitment — MPR exclusive
- MIF-(CD74+CD44) toward Lipid-associated : extended MIF signal

non-MPR dominant interactions:
- PPIA-BSG expanded toward Stress-response + IFN-stimulated : stress/immunosuppression toward pro-immunogenic subtypes

**CellChat TAMs → CD8 (MPR vs non-MPR) :**

Shared interactions: HLA-A/B-CD8A → TEX (both conditions)

MPR dominant interactions:
- HLA-A/B/C-CD8A toward TPEX (weak) : MHC I toward precursors
- HLA-C-CD8A toward TEX : third HLA allele MPR exclusive
- MIF-(CD74+CD44) + SPP1-CD44 toward TPEX : survival/retention signals

non-MPR dominant interactions:
- HLA-A/B-CD8A toward TPEX (strong) : enriched toward precursors, attempted antigen presentation in hostile context

Bidirectional loop conclusion CD8↔TAMs CellChat:
MPR = MIF basal + CD4 recruitment (HLA-DRB1-CD4) + resolution signals toward IFN-stimulated. 
non-MPR = PPIA-BSG stress expanded + HLA-A/B enriched toward TPEX without functional recovery.

**CellChat pCR observations (MPR vs pCR) :**
Note: single dataset, pseudobulk MPR/pCR non-significant in Blocs 3-4.

CD8 → TAMs:
- Shared: MIF-(CD74+CD44), PTPRC-MRC1 (4 subtypes), ANXA1-FPR1, CD99-PILRA toward IFN-stimulated, MIF-(CD74+CXCR4) toward Resident M2
- MPR exclusive: HLA-DRB1-CD4 (CD4 recruitment lost in pCR), MIF-(CD74+CD44) toward Lipid-associated, PPIA-BSG toward Lipid-associated
- pCR exclusive: ANXA1-FPR1 extended to Stress-response, PGE2-PTGES3-PTGER4 toward Resident M2 (prostaglandin resolution signal)

TAMs → CD8:
- Shared: HLA-A/B-CD8A toward TEX, MIF-(CD74+CXCR4) toward TPEX (MPR) and TPEX+TEX (pCR)
- MPR exclusive: HLA-C-CD8A toward TEX, MIF-(CD74+CD44) + SPP1-CD44 toward TPEX (weak)
- pCR exclusive: APP-CD74, CXCL16-CXCR6, LGALS9-CD45 toward TPEX+TEX, pan-HLA strong toward TEX, MIF-(CD74+CD44) strong toward TPEX+TEX

**Cross-tool convergence CellChat/LIANA+ GSE243013 :**
Confirmed by both tools:
- ANXA1-FPR1 → IFN-stimulated MPR
- CD99-PILRA → IFN-stimulated MPR 
- HLA-A/B-CD8A → TEX (both conditions) 
- HLA-A/B-CD8A enriched → TPEX non-MPR 
LIANA+ only: CCL5-ACKR1/CXCL13-ACKR1, CCL5-CCRL2, TIGIT-NECTIN2, HLA-DRA-LAG3, CXCL9/10/11-CXCR3
CellChat only: MIF signaling, LGALS9-CD45, CXCL16-CXCR6, pCR-specific signals

### General conclusion : Bloc 5 intercellular communication

CellChat v2 constitutes the primary differential communication analysis tool for both datasets, given its established probabilistic framework and curated ligand-receptor-cofactor database. However, LIANA+ provides critical complementary signals across all biological axes — including the CD8-TAMs narrative (CCL5-CCRL2 decoy diversion, GZMB-IGF2R misdirected cytotoxicity, HLA-LILRB checkpoints, LAG3 differential engagement) and the epithelial axes in GSE207422 (HSF1 targets, NECTIN4-TIGIT), that are not detected by CellChat due to database resource differences.

Cross-tool confirmed interactions, most robust findings:
GSE207422 (6): CCL5-CCR1 (CD8→TAMs), HLA-A/C-CD8A (TAMs→CD8),GZMA-PARD3 + CRTAM-CADM1 (CD8→Epithelial), HLA-A-CD8A (Epithelial→CD8), FN1-SDC1 (TAMs→Epithelial)
GSE243013 (4): ANXA1-FPR1 + CD99-PILRA (CD8→TAMs MPR IFN-stimulated), HLA-A/B-CD8A toward TEX (TAMs→CD8 both conditions), HLA-A/B-CD8A enriched toward TPEX (TAMs→CD8 non-MPR)

HSF1 cross-dataset observation:
HSF1 transcriptional activity is enriched in non-MPR/NMPR CD8 and TAM compartments (CollecTRI, Blocs 3-4) but no HSF1 direct target is detected
as a ligand in CD8↔TAMs interactions in either dataset. HSF1 influences the CD8↔TAMs axis at the intracellular transcriptional level but not at
the intercellular ligand-receptor level. The intercellular link between HSF1 intracellular activity and the communication programs observed remainsto be resolved a question amenable to spatial transcriptomics.

### Bloc 6 — Bulk validation — not performed

A systematic search for publicly available bulk RNA-seq datasets of neoadjuvant anti-PD-1 treated NSCLC with MPR/non-MPR pathological 
response classification was conducted. The following options were evaluated and excluded:

1. **GSE126044** (Kim et al., Genome Medicine 2020, nivolumab, n=16) : RECIST-based classification (CR/PR vs SD/PD), not pathological response (MPR/non-MPR). RECIST underestimates pathological response, in Checkmate 159, 45% achieved MPR while only 10% showed partial response on CT/RECIST. Classification mismatch prevents valid comparison with scRNA-seq findings.

2. **IMvigor210** : urothelial carcinoma cohort (atezolizumab/anti-PD-L1) , wrong cancer type and wrong drug class, not applicable.

3. **Hu et al. Genome Medicine 2023 bulk cohort** (n=21) : pre-treatment biopsies collected before neoadjuvant ICI, baseline context, not post-treatment surgical resection. Response prediction context differs from post-neoadjuvant MPR classification.

4. **HSF1 signature** : preliminary analysis in both scRNA-seq datasets showed HSF1-target genes are not sufficiently discriminant between MPR and non-MPR in CD8 or TAM compartments to construct a robust ssGSEA signature (majority of HSF1-targets expressed in both conditions; only SESN3 cross-dataset exclusive non-MPR in CD8).

**Conclusion:** No publicly available bulk RNA-seq dataset combining NSCLC, neoadjuvant anti-PD-1, and MPR/non-MPR pathological response 
classification was identified. Bulk validation remains a priority for future work, contingent on availability of harmonized datasets 
or institutional cohorts. Signatures retained for future bulk validation: CD8 Exhaustion, CD8 Cytotoxicity, M2_immunosuppressive, SPP1, IFN_response — constructed from UCell modules applied in Blocs 3 and 4.

## Methodological Notes

### TME annotation : iterative approach
Three automated methods tested for global TME annotation:

1. **Azimuth** (Hao et al., Nature Biotechnology 2024) — installation failed on 
   Windows (heavy genomic dependencies). Attempted via WSL/Ubuntu 24.04 — 
   RunAzimuth() blocked by RAM constraints (16GB). Server execution pending.

2. **SingleR** (HumanPrimaryCellAtlas reference) — applied successfully, 
   used for cross-validation of manual annotation.

3. **scType** (Ianevski et al., Nature Communications 2022) — lung tissue 
   reference insufficient granularity for immune subtypes. 
   Used for myeloid cluster validation only.

**Final annotation:** manual (top 50 markers/cluster),cross-validated by SingleR. 
Azimuth validation pending.

---

## Repository structure

```
TumorImmune_Crosstalk_LUAD/
├── README.md
├── Scripts/
│   ├── BLOC 0_Data_acquisition/
│   ├── BLOC 1_QC_and_preprocessing/
│   ├── BLOC 2_Global_TME_Annotation/
│   ├── BLOC 3_CD8_Tcells_Focus/
│   └── BLOC 4_Immunosuppressive_TME_Compartment/
│       ├── BLOC 4A_GSE243013/
│       └── BLOC 4B_GSE207422/
├── Data/
│   ├── collectri_network.csv
│   └── metadata_LUAD.csv
├── Objects/
├── Results/
│   ├── Figures/
│   └── Tables/
├── session_info.txt
└──session_info_linux.txt
 
```

---

## Relationship to prior work

This project is a direct continuation of [CD8_NSCLC_scRNAseq](https://github.com/yasmina-bioinfo/CD8_NSCLC_scRNAseq), which characterized CD8 T cell heterogeneity across two NSCLC datasets (GSE131907, GSE207422) and identified a CD8_Exhausted_Terminal enrichment in MPR patients (OR=3.36, p_adj<0.001) and a STAT2-high exhaustion program in non-MPR patients using DoRothEA.

The present project extends this work by:
- Validating CD8 exhaustion findings in a larger independent cohort (GSE243013, n=63 LUAD vs n=13)
- Expanding the analysis to the TAM and malignant epithelial compartments
- Incorporating paired TCR sequencing to validate CD8 state tumor-reactivity
- Upgrading TF inference from DoRothEA to CollecTRI for improved regulon coverage
- Objectively identifying ELK4 as the most robust cross-dataset TF signal in MPR CD8.TEX by variance-based selection, STAT2 identified by DoRothEA was not reproduced by CollecTRI objective ranking, suggesting tool-specific sensitivity differences
- Adding TAM and epithelial compartment analyses absent from the prior work

---

## Author

**Myriam Yasmina Soumahoro**
MSc Biology, University of Geneva | Mandela Washington Fellow 2025, Arizona State University
[GitHub](https://github.com/yasmina-bioinfo)


 
 