# TumorImmune_Crosstalk_LUAD

> **How does tumor cell plasticity shape CD8 T cell functional states and their response to anti-PD-1 in lung adenocarcinoma?**

This project investigates the crosstalk between malignant epithelial cells and CD8 T cells in the tumor microenvironment (TME) of lung adenocarcinoma (LUAD), using neoadjuvant chemo-immunotherapy data with matched pathological response (MPR vs. non-MPR; MPR vs pCR and non-MPR vs pCR).

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
- Manual annotation + automated validation (`SingleR`)

### Block 3 : CD8 T cell focus
- Robust state annotation on reference atlas (`ProjecTILs`)
- Functional module scoring: exhaustion, cytotoxicity, memory (`UCell`)
- TCR integration: clonotype expansion, repertoire diversity (`scRepertoire`)
- MPR vs. non-MPR, pCR vs MPR, and non_MPR vs pCR comparisons of CD8 state composition

### Block 4 : Malignant epithelial focus
- Isolation and sub-clustering of epithelial population
- Malignancy confirmation (`CopyKAT`)
- Functional scoring: EMT, IFN stress, proliferation (`UCell`)
- MPR vs. non-MPR comparison of tumor cell programs

### Block 5 : TME intercellular communication
- CellChat object stratified by MPR vs. non-MPR (`CellChat`)
- Global signaling network comparison
- Focus on epithelial ↔ CD8 ligand-receptor axes

### Block 6 : Transcription factor activity
- TF activity inference on CD8 and epithelial compartments separately
- Network: CollecTRI (successor to DoRothEA, 12-source curated GRN)
- Method: `run_ulm` via `decoupleR`
- Differential TF activity: MPR vs. non-MPR

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
| CellChat | Cell-cell communication | Jin et al., *Nat Comm* 2021 |
| decoupleR | TF activity inference | Badia-i-Mompel et al., *Bioinformatics Advances* 2022 |
| CollecTRI | TF regulon network | Müller-Dott et al., *Nucleic Acids Research* 2023 |
| SingleR | Automated annotation | Aran et al., *Nat Immunol* 2019 |
| sctype | Automated cell type annotation | Ianevski et al., *Nature Communications* 2022 |
| Azimuth | Reference-based annotation (not installable on Windows — replaced by sctype) |  Hao et al., *bioRxiv* 2022 |

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
- Script 01: FindAllMarkers complete with 20 clusters annotated manually (top50 markers)
    - 18/20 clusters annotated with high confidence
    - 2/20 clusters with medium confidence (Cluster 5: CD4 Tfh/exhausted, Cluster 9: NKT/γδ)
- Script 02: SingleR complete (HumanPrimaryCellAtlas, label.main + label.fine)
- Script 03: sctype complete (Lung tissue reference; limited resolution for immune subtypes)
- Script 04: Azimuth, installed in WSL (Ubuntu 24.04, R 4.6.0) after 3h installation
  - lungref reference (584,884 cells) loaded successfully
  - RunAzimuth blocked by RAM constraint (16GB), server execution pending
- Script 05: Consensus annotation complete with 20 clusters annotated, 3 pending Azimuth confirmation
  (Clusters 5, 9, 13, 18)
- Script 06: TME visualization complete
  - UMAP final annotation + split by response
  - Barplot proportions (Chi-2 p < 2.2e-16)
  - Fisher post-hoc: MPR vs non-MPR, pCR vs non-MPR, MPR vs pCR; post-hoc pairwise comparisons saved in Results/Tables/Bloc2_fisher_posthoc.csv

### Bloc 3 : T cells analysis
- Script 01: T cells subsetting (clusters 1,2,3,4,5,9,10,11,17) : 172,110 cells
  ElbowPlot inspection = 30 PCs retained
- Script 02: Harmony + clustering (resolution=0.4) + UMAP = 16 T cell clusters
- Script 03: ProjecTILs annotation (human CD8 reference)
  - 108,456/172,110 cells (63%) filtered by scGate (non-CD8 pure)
  - 63,654 cells projected = 7 CD8 states identified:
    CD8.EM (25,575), CD8.CM (18,181), CD8.TEX (15,115), CD8.TPEX (2,011),
    CD8.NaiveLike (1,342), CD8.MAIT (1,273), CD8.TEMRA (157)
  - STACAS alignment failed (RAM) = direct projection used
  - CD8.TEX enriched in non-MPR, consistent with portfolio narrative

  - **METHODOLOGICAL NOTE**: ProjecTILs was run before canonical marker-based annotation.
  The recommended workflow is to first identify CD8 clusters by canonical markers, then project only confirmed CD8 cells onto the reference atlas.
  Running ProjecTILs on the full T cell subset (including CD4, Tregs, NKT, MAIT) caused scGate to incompletely filter non-CD8 populations, resulting in artifactual CD8 labels for non-CD8 clusters.
  **For future analyses: subset confirmed CD8 clusters first, then run ProjecTILs.**
  Final annotation relies on canonical markers for non-CD8 clusters.

- Script 07 : CD8 T cells subsetting (corrected workflow)
  - **METHODOLOGICAL CORRECTION**: Script 03 ran ProjecTILs on full T cell subset (172,110 cells including CD4, Tregs, NKT) before canonical marker annotation.
  This caused scGate to incompletely filter non-CD8 populations (63% removed) and generated artifactual CD8 labels for non-CD8 clusters.
  - Corrected workflow: canonical marker annotation (Script 04) performed first to identify confirmed CD8 clusters (1, 2, 5, 6, 10), then ProjecTILs applied to CD8 pure subset only.
  - CD8 confirmed clusters: 75,622 cells
  - ElbowPlot inspection = 25 PCs retained (vs 30 for full T cell subset)

- Script 08 : ProjecTILs on confirmed CD8 cells (corrected workflow)
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

- Script 09 : CD8 ProjecTILs barplot (proportions by response)
  - Chi-2 p < 2.2e-16 : CD8 state distribution highly significantly different across groups
  - CD8.TEX: non-MPR (~30%) > MPR (~25%) > pCR (~15%): gradient confirmed statistically
  - CD8.CM: pCR (~40%) = MPR (~35%) > non-MPR (~30%) : memory enriched in responders
  - CD8.EM: similar across groups (~35-40%) : functional status to be assessed by UCell
  - CD8.TPEX: visible in pCR only, consistent with reactivation hypothesis
  - CD8.MAIT: 47 cells (~0.08%), negligible, not visible on barplot. Reduction from 1,273 (Script 03) to 47 (Script 08) confirms cleaner CD8 population
  - Fisher post-hoc results saved in Results/Tables/Bloc3_CD8_fisher_posthoc.csv

- Script 10 : UCell scoring on CD8 T cells
  **Exhaustion scores (mean), CD8.TEX:**
    - non-MPR = 0.347 > pCR = 0.249 > MPR = 0.220

  **Exhaustion scores (mean), CD8.TPEX:**
    - non-MPR = 0.268 > pCR = 0.225 > MPR = 0.197

- Script 11: scRepertoire TCR analysis (WSL, R 4.6.0)
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

- Script 12 : CollecTRI TF activity on CD8 T cells
  - Tool: decoupleR run_ulm + CollecTRI network (43,159 interactions, 1,186 TFs)
  - RAM constraint: full CD8 object (57,587 cells) requires 14.6 GiB
  Solution: restricted to CD8.TEX, CD8.TPEX, CD8.EM, max 10,000 cells/state
  Total analyzed: ~21,821 cells, set.seed(42) for reproducibility
  - CollecTRI network saved locally from RStudio (OmnipathR issues in WSL)
  - Top 20 TFs by variance: HSF1, HSF2, RFXAP, RFXANK, ELK4, MYC, RFX5, NFYC,
  NFKB, RLF, DMTF1, RELA, CIITA, NFYB, MLXIP, JUN, HOPX, NFKB1, DAXX, TBX21
  - See TF_biological_roles.md for full biological annotation of each TF

  **TFs common with portfolio (GSE207422 DoRothEA analysis):**
  - ELK4 : enriched in MPR CD8.TEX (confirms portfolio finding)
  - TBX21 : present in top 20 (discordance with portfolio, see below)
  - STAT2, STAT1, ELK1, IRF7 : present in CollecTRI but not in top 20 by variance


### Preliminary observations : TME composition (Bloc 2 barplot) and CD8 states analysis (in progress, Bloc 3)

- CD8.TEX visually more abundant in non-MPR and pCR than MPR on UMAP split
- IMPORTANT: non-MPR has 42 patients vs 10 MPR and 11 pCR, absolute cell counts are not directly comparable. Proportional analysis (UCell, next scripts) required before biological conclusions.
- Apparent higher CD8.TEX in pCR vs MPR suggests pCR may harbor more reactivatable exhausted CD8, consistent with complete tumor clearance (pCR = 0% residual tumor)
- Portfolio GSE207422: CD8_Exhausted_Terminal enriched in MPR (OR=3.36): different tool (manual annotation vs ProjecTILs), different dataset (n=8 patients vs n=63), different patient proportions (GSE207422: MPR > non-MPR cells; GSE243013: non-MPR >> MPR), and different annotation granularity. Visual impression ≠ statistical enrichment.
Discordance to be resolved by proportional analysis and UCell scoring.

### Preliminary observations : CD8 T cell states (Bloc 3, Script 08, corrected workflow)
- CD8.TEX enriched in non-MPR (visual), consistent with portfolio (OR=3.36, GSE207422)
- CD8.TPEX enriched in pCR vs MPR, suggests reactivation of exhausted precursors
- CD8.TEX > in MPR than pCR; MPR cells have progressed further toward exhaustion but retain partial function (residual tumor ≤ 10%)
- Proposed gradient: non-MPR (TEX dominant) to MPR (TEX + TPEX) to pCR (TPEX dominant)
  Hypothesis: anti-PD-1 reactivated a fraction of TEX toward TPEX state in responders
  This connects with portfolio finding: CD8_Exhausted_Terminal enriched in MPR (OR=3.36)
- CD8.EM numerically dominant in non-MPR; NOTE: high EM count does not imply functional efficacy. EM cells in non-MPR may be dysfunctional or evaded by tumor. Functional scoring (UCell, next scripts) required to assess EM functional state.
- IMPORTANT: all observations are visual/preliminary; statistical proportional
  analysis and UCell scoring required before biological conclusions.
- Connection with portfolio (GSE207422): STAT2-high program in non-MPR associated with differentiation blockade toward cytotoxic effector state. High EM count in non-MPR may reflect this blockade; cells stalled in EM state, unable to fully differentiate. Immunosuppressive TME (TREM2+ TAMs, CCR8+ Tregs) compounds CD8 dysfunction. Anti-PD-1 partially relieves this blockade in MPR/pCR.

### Preliminary observations : UCell CD8 scoring (Script 10)
- pCR shows higher exhaustion scores than MPR in both TEX and TPEX = pCR cells were exposed to stronger antigenic pressure initially; but their reactivation capacity was superior, enabling complete tumor clearance
- CD8.EM in non-MPR more exhausted (0.128) than MPR (0.107) and pCR (0.106) = confirms EM dysfunction in non-responders, consistent with STAT2-mediated differentiation blockade identified in portfolio (GSE207422)
- TPEX retain higher cytotoxicity score than TEX across all groups = TPEX preserve more residual function, consistent with reactivable precursor state
- Memory score highest in CD8.NaiveLike

**Hypotheses to test in next blocs:**
- H1 (CollecTRI): pCR T cells have superior intrinsic reactivation capacity
- H2 (CellChat): immunosuppressive TME composition determines reactivation failure
- H3: interaction between intrinsic T cell state and TME context

### Preliminary observations : scRepertoire (Script 11)
- CD8.TEX high expansion (85%) confirms antigen-specific tumor-reactive identity despite terminal exhaustion: these cells have been activated and expanded
- pCR shows more non-expanded clones than MPR/non-MPR,consistent with higher clonal diversity, polyclonal response may enable complete tumor clearance
- Clonal diversity gradient: MPR > pCR > non-MPR: responders maintain broader TCR repertoire

### Preliminary observations : CollecTRI (Script 12)

Heatmap observations (TFs ordered by activity strength per state x condition):

| State | pCR | MPR | non-MPR |
|---|---|---|---|
| CD8.EM | NFYB, NFYC, HSF2, DAXX, CIITA, DMTF1, ELK4, RLF | DAXX, ELK4 (~2), RLF | CIITA, RFXAP, DAXX, RFXANK, RFX5, NFYB |
| CD8.TEX | MLXIP (~1), HOPX (~1.5) | TBX21 (~1.5), RFXANK, RFX5, RFXAP | TBX21, RFXANK, RFX5, CIITA, RFXAP, NFKB1 |
| CD8.TPEX | HSF2, NFYC, MLXIP, HOPX, HSF1, DMTF1, JUN, RELA | JUN, NFKB, ELK4, RLF, DMTF1, RELA, NFKB1, MYC | MYC, NFKB1, NFKB, RELA, HSF2, HSF1, NFYB, NFYC, DMTF1 |

Violin observations (key TFs):
- STAT2: non-MPR highest in CD8.TEX and CD8.TPEX, confirms STAT2-high IFN program (same as portfolio)
- STAT1: non-MPR > MPR ≈ pCR in TEX and TPEX, co-activated with STAT2
- ELK4: MPR dominant in CD8.TEX, confirms portfolio finding
- ELK1: pCR > MPR in TEX; pCR > non-MPR in TPEX
- TBX21: non-MPR dominant in TEX, discordance with portfolio (MPR enriched)
  Possible explanations: different tool (DoRothEA vs CollecTRI), different dataset, or TBX21 activates different programs depending on TME context
- IRF7: similar across all three groups in TEX and TPEX

**Key biological interpretation:**
- STAT2/STAT1-high program in non-MPR CD8.TEX and CD8.TPEX confirmed across two independent cohorts and two different TF inference tools, robust signal
- ELK4 enriched in MPR CD8.TEX confirmed, cytotoxic effector engagement
- pCR CD8.TPEX shows richest transcriptional program (HSF2, NFYC, JUN, RELA) consistent with superior reactivation capacity
- TBX21 discordance to be resolved by CollecTRI analysis on full CD8 object on compute server
The apparent discordance in TBX21 activity between cohorts is likely explained by differences in clinical group granularity rather than a biological contradiction. In GSE207422, the MPR category (≤10% residual tumor) may have included "near-pCR" patients with near-complete responses, in whom TBX21 co-activation with ELK4 represented a coordinated cytotoxic effector program. In GSE243013, where pCR is separated from MPR, the MPR category is more homogeneous. TBX21 activity in non-MPR CD8.TEX without ELK4 co-activation, may reflect an abortive cytotoxic program: TBX21 activation insufficient to drive full effector differentiation in the absence of its co-activator. This interpretation suggests that ELK4 may be the key discriminating TF between functional and dysfunctional cytotoxic programs, with TBX21 as a necessary but insufficient partner.
Additionally, the smaller and imbalanced patient cohort in GSE207422 (MPR n=3, non-MPR n=10) may have introduced sampling bias in TF activity estimates. With only 3 MPR patients, the TBX21 signal may have been driven by one or two outlier patients with atypically high TBX21 activity, rather than reflecting a true MPR-specific program.

## Methodological Notes

### Automated annotation : iterative approach
Three automated annotation methods were tested for global TME annotation:

1. **Azimuth** (Hao et al., *bioRxiv* 2022), initially planned as primary tool. 
Installation failed on Windows due to heavy genomic dependencies (BSgenome.Hsapiens.UCSC.hg38, EnsDb.Hsapiens.v86). Documented as Windows limitation.

2. **sctype** (Ianevski et al., *Nat Commun* 2022), selected as lightweight alternative. 
Lung tissue reference insufficiently granular for TME immune subtype resolution. 
Useful for macrophage/myeloid validation only.

3. **Azimuth** via WSL/VS Code, reinstalled under Ubuntu 24.04 LTS to bypass. 
Windows dependency constraints. lungref reference (584,884 cells) loaded successfully. 
RunAzimuth() blocked by insufficient local RAM (16GB), server execution pending.

Manual annotation (top50 markers per cluster) remains the primary reference, validated by SingleR (HumanPrimaryCellAtlas) and partially by sctype (myeloid clusters). 
Azimuth validation pending server execution.

---

## Repository structure

```
TumorImmune_Crosstalk_LUAD/
├── README.md
├── Scripts/
│   ├── BLOC 0_Data acquisition
│   ├── BLOC 1_QC and preprocessing
│   ├── BLOC 2_Global_TME_Annotation
│   ├── BLOC 3_CD8_Tcells_Focus
│   ├── 
│   ├── 
│   └── 
├── Data/
│   └── metadata_LUAD.csv
├── Figures/
└── Results/
```

---

## Relationship to prior work

This project is a direct continuation of [CD8_NSCLC_scRNAseq](https://github.com/yasmina-bioinfo/CD8_NSCLC_scRNAseq), which characterized CD8 T cell heterogeneity across two NSCLC datasets (GSE131907, GSE207422) and identified a CD8_Exhausted_Terminal enrichment in MPR patients (OR = 3.36, p_adj < 0.001) with  and a STAT2-high exhaustion program in non-MPR patients.

The present project extends this work by:
- Validating CD8 exhaustion findings in a larger independent cohort (n = 61 LUAD vs. n = 8)
- Adding the malignant epithelial compartment to investigate tumor-immune crosstalk
- Incorporating paired TCR sequencing to validate T cell state annotations
- Upgrading TF inference from DoRothEA to CollecTRI for improved regulon coverage
- Testing whether the STAT2-high (non-MPR) and ELK4/ELK1/TBX21-high (MPR) transcriptional programs identified in GSE207422 are reproducible in GSE243013

---

## Author

**Myriam Yasmina Soumahoro**   
[GitHub](https://github.com/yasmina-bioinfo)
