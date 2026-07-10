# TumorImmune_Crosstalk_LUAD

## Publication status 
> ***Preprint in preparation** : Multi-cohort single-cell analysis of CD8-TAM crosstalk and anti-PD-1 response in NSCLC.
> Targer: bioRxiv : manuscript in preparation.

> ***How do CD8 T cell exhaustion programs and tumor-associated macrophage transcriptional diversity and plasticity jointly shape response to neoadjuvant anti-PD-1 therapy in NSCLC?**

This project investigates transcription factor activity, functional programs and intercellular communication across CD8 T cells and tumor-associated macrophages (TAMs) in the tumor microenvironment (TME) of lung adenocarcinoma (LUAD), using two independent neoadjuvant chemo-immunotherapy cohorts with matched pathological response (MPR vs. non-MPR; MPR vs pCR and non-MPR vs pCR in GSE243013).

Although histological subtype was not explicitly reported by Hu et al. 2023, the presence of AT2 epithelial cells (Tumor_epithelial_AT2) in the TME annotation is consistent with a predominantly LUAD composition, as AT2 cells represent the recognized cell of origin of lung adenocarcinoma (Sainz de Aja et al. 2021; Xing et al. 2021). KRT5 and TP63 expression, canonical squamous markers, was restricted to the basal epithelial 
cluster and absent in AT2 cells, further supporting an adenocarcinoma histology. Cross-cohort comparison was therefore restricted to LUAD patients in GSE243013 (n=63), excluding LUSC and other NSCLC subtypes, to minimize histological confounding. LUAD and LUSC exhibit distinct tumor immune microenvironments, with LUAD showing higher CD8 T cell infiltration and exhaustion signatures, ensuring that the observed CD8 
exhaustion and TAM immunosuppressive signals reflect LUAD-specific biology.

This work builds on independent prior single-cell analyses of CD8 T cell exhaustion in NSCLC ([CD8_NSCLC_scRNAseq](https://github.com/yasmina-bioinfo/CD8_NSCLC_scRNAseq)), extended here into a full cross-dataset analysis integrating transcription factor activity inference (CollecTRI), functional signature scoring (UCell) and intercellular communication (LIANA+/CellChat) across two independent NSCLC neoadjuvant cohorts.

---

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
- Functional program scoring: unbiased Hallmark MSigDB gene set scoring (`UCell`), 22 biologically relevant discriminant signatures identified
- TCR integration: clonotype expansion, repertoire diversity (`scRepertoire`)
- TF activity inference: Top 20 TFs by between-group variance (subtype × response means) across CD8 states (`CollecTRI/decoupleR`), with Wilcoxon significance testing per state
- MPR vs. non-MPR comparisons of CD8 state composition (both cohorts); MPR vs pCR patient-level comparison performed in GSE243013 only, reserved for perspectives
- Cross-cohort validation: ProjecTILs, UCell Hallmark and CollecTRI applied to both GSE243013 and GSE207422 CD8 T cells (completed)

### Block 4 : Immunosuppressive TME Compartment 
- TAM subtype annotation and reclustering (GSE243013 + GSE207422)
- Functional program scoring: unbiased Hallmark MSigDB gene set scoring (`UCell`), 25 (GSE243013) and 23 (GSE207422) biologically relevant discriminant signatures identified
- Patient-level pseudobulk analysis MPR vs pCR (`UCell`, Hallmark signatures)
- TF activity inference on TAMs: Top 20 TFs by between-group variance (subtype × response means) across TAM subtypes (`CollecTRI/decoupleR`), with Wilcoxon significance testing per subtype, both cohorts (completed)
- Malignant epithelial analysis: CopyKAT, SCEVAN, CytoTRACE 2, UCell, CollecTRI (GSE207422 only), deferred to perspectives/limitations, not part of the main narrative

### Block 5 : TME intercellular communication
- Intercellular communication inference across full TME using LIANA+ (consensus across CellChatDB, OmniPath and other LR databases)
- Differential cell-cell communication analysis across MPR vs. non-MPR vs. pCR using CellChat v2 (MultiNicheNet evaluated but abandoned: n=3 MPR patients insufficient for pseudobulk-based differential testing)
- Focus on TAM ↔ CD8 ligand-receptor axes, prioritized on TF and functional programs established in Blocks 3 and 4
- Full TME opened to detect unexpected interactions; TAM ↔ Epithelial and CD8 ↔ Epithelial axes deferred to perspectives, since epithelial malignant cell analysis was performed in GSE207422 only, precluding a cross-dataset comparison consistent with the rest of the main narrative

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
| MSigDB Hallmark | Standardized gene set signatures for UCell functional scoring | Liberzon et al., *Cell Systems* 2015 |
---


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

#### Script 10 : UCell scoring on CD8 T cells (custom signatures)
  **Exhaustion scores (mean), CD8.TEX:**
    - non-MPR = 0.347 > pCR = 0.249 > MPR = 0.220

  **Exhaustion scores (mean), CD8.TPEX:**
    - non-MPR = 0.268 > pCR = 0.225 > MPR = 0.197

#### Script 10b : UCell Hallmark scoring on CD8 T cells (all 50 signatures)
  - Unbiased approach — all 50 Hallmark MSigDB signatures tested
  - Wilcoxon MPR vs non-MPR with BH correction
  - 21 biologically relevant discriminant signatures retained (p_adj < 0.05)
  - Full results : Results/Tables/Bloc3_GSE243013_UCell_Hallmark_CD8_wilcox_results.csv

#### Script 10c : UCell Hallmark preprint figures
  - Main figure : heatmap 21 signatures (scale=row, p-value asterisks per signature)
  - Biologically filtered : non-relevant signatures excluded (ADIPOGENESIS, HEDGEHOG_SIGNALING, ESTROGEN_RESPONSE, ANDROGEN_RESPONSE, PANCREAS_BETA_CELLS etc.)
  - Figure : Results/Figures/CD8/Preprint/Bloc3_GSE243013_Hallmark_CD8_heatmap_main.png

#### Script 11: scRepertoire TCR analysis (WSL, R 4.6.0)
  - NOTE: scRepertoire requires gsl >= R 4.5.0, not available on Windows R 4.4 Installed and run in WSL environment
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
**Superseded:** this Top 20 selection (per-cell variance, with downsampling) and the observations above, including the TBX21+ELK4 "functional cytotoxic program" narrative, RLF, and DMTF1, are superseded by Script 12b (below), which uses the corrected between-group variance selection method and adds formal Wilcoxon significance testing. See "Preliminary observations : CollecTRI CD8 (Bloc 3)" for the current, statistically verified findings.

#### Script 13 : Pseudobulk Hallmark UCell scores MPR vs pCR (patient-level), CD8 T cells

- Input: Objects/Bloc3_GSE243013_12c_seu_CD8_Hallmark.rds
- Signatures: the same 21 Hallmark MSigDB signatures established as discriminant between MPR and non-MPR (not the original custom UCell signatures), restricted to CD8.TEX/CD8.TPEX
- Subset: MPR (n=10 patients) and pCR (n=11 patients); non-MPR excluded
- Pseudobulk: mean UCell score per patient per signature
- Wilcoxon rank-sum test MPR vs pCR (patient-level), BH-corrected across the 21 signatures
  References: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)

**Superseded:** earlier figures used four custom UCell signatures (Cytotoxicity, Exhaustion, Memory, TPEX), all non-significant (p>0.05). This approach has been replaced to align with the Hallmark-only, BH-corrected methodology used throughout the project.

#### Script 14 : CollecTRI pseudobulk MPR vs pCR, CD8, GSE243013

- Purpose: extend the non-MPR vs MPR CollecTRI discriminant TF panel (Script 12b) to a direct MPR vs pCR comparison, to characterize pCR relative to MPR at the transcription factor level
- Input: Objects/Bloc3_12_tf_acts_raw.rds (raw per-cell TF scores, already computed)
          Objects/Bloc3_12b_seu_CD8_TF.rds (metadata: sampleID, PathResponse)
- Same Top 20 TFs as Script 12b, restricted to MPR and pCR patients
- Patient-level pseudobulk (mean TF score per patient per CD8 state), Wilcoxon rank-sum MPR vs pCR, BH-corrected within each CD8 state
- Output: Results/Tables/Bloc3_14_CollecTRI_pseudobulk_MPR_vs_pCR.csv

### Script 12b : CollecTRI CD8 T cells GSE243013 (no downsampling)
- **Input:** Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
- **Subset:** CD8.TEX and CD8.TPEX only (CD8.EM excluded, secondary to narrative, RAM constraint)
- **Downsampling:** None, full subset used (TEX + TPEX, ~25,000 cells)
- **RAM note:** Script 12 required downsampling to 10,000 cells/state due to WSL RAM constraints. Script 12b runs on 12Gi WSL allocation with TEX+TPEX only.
- **TF selection method:** corrected from per-cell variance to between-group variance (variance of state × response means), consistent with the project's comparative design, see Methodological correction note above.
- **Output:**
  - Results/Tables/Bloc3_CollecTRI_TF_activity.csv (Top 20 TF means per state × condition)
  - Results/Tables/Bloc3_CollecTRI_TF_wilcox_bystate.csv (Wilcoxon MPR vs non-MPR, per state, BH-corrected)
  - Results/Figures/CD8/Preprint/Bloc3_CollecTRI_heatmap_MPR_nonMPR.png (2-condition heatmap, preprint)
  - Results/Figures/CD8/Preprint/Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR.png (2-condition violin, preprint)
  - Results/Figures/CD8/Exploratory/Bloc3_CollecTRI_heatmap_MPR_nonMPR_pCR.png (3-condition heatmap, exploratory only, not in preprint)
  - Results/Figures/CD8/Exploratory/Bloc3_CollecTRI_key_TFs_violin_MPR_nonMPR_pCR.png (3-condition violin, exploratory only, not in preprint)
  - Objects/Bloc3_12_tf_acts_raw.rds (raw per-cell TF scores, saved to avoid re-running decoupleR if further adjustments are needed)
- **Top 20 TFs by between-group variance:** MTF1, ELK4, HOPX, HSF4, HSF2, FOXA1, HSF1, NFAT5, NFYB, NFYC, FOXF2, KLF4, MYC, MLX, RFX5, RFXANK, RFXAP, MLXIP, DAXX, FOXO3
- **Top 6 TFs by between-group variance (violin plots):** HSF1, HSF2, NFYC, MLXIP, FOXO3, HOPX

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

**Superseded:** this analysis used four custom UCell signatures (Exhaustion, Cytotoxicity, Memory, TPEX) rather than the Hallmark-only, unbiased approach adopted for the main narrative. Superseded by Script 03b/03c (below), which use the same 22 Hallmark MSigDB signatures established as discriminant for this dataset.

#### Script 03b : UCell Hallmark scoring on CD8 T cells (all 50 signatures)
  - Unbiased approach, all 50 Hallmark MSigDB signatures tested
  - Wilcoxon MPR vs NMPR with BH correction
  - 22 biologically relevant discriminant signatures retained (p_adj < 0.05)
  - Full results : Results/Tables/Bloc3_GSE207422_UCell_Hallmark_CD8_wilcox_results.csv

#### Script 03c : UCell Hallmark preprint figures
  - Main figure : heatmap 22 signatures (scale=row, p-value asterisks per signature)
  - Biologically filtered : non-relevant signatures excluded (ADIPOGENESIS, HEDGEHOG_SIGNALING, ESTROGEN_RESPONSE, ANDROGEN_RESPONSE, PANCREAS_BETA_CELLS, SPERMATOGENESIS, MYOGENESIS, NOTCH_SIGNALING etc.)
  - Figure : Results/Figures/CD8/Preprint/Bloc3_GSE207422_Hallmark_CD8_heatmap_main.png  

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

#### Script 04b : CollecTRI TF activity on CD8 T cells (GSE207422)
  - Tool: decoupleR run_ulm + CollecTRI network (local CSV)
  - Subset: CD8.TEX, CD8.TPEX only (CD8.EM excluded, consistent with GSE243013 narrative scope)
  - Downsampling: none, full subset used (corrected from earlier Script 04, which downsampled to 10,000 cells/state)
  - TF selection method: between-group variance (variance of state × response means), corrected from per-cell variance, consistent with GSE243013 methodology
  - Top 20 TFs by between-group variance: MYC, RFXAP, RFXANK, RFX5, SP1, NFKB, ELK4, HIF1A, E2F4, IRF6, SRSF2, TP53, HOPX, E2F1, PAWR, NAB2, STAT3, NFKB1, JUN, TFDP1
  - Top 6 TFs by between-group variance (violin plots): MYC, RFXAP, RFXANK, RFX5, SP1, NFKB
  - Statistical test: Wilcoxon MPR vs NMPR, per CD8 state, BH-corrected within state
  - Output: Results/Tables/Bloc3_GSE207422_CollecTRI_TF_activity.csv (Top 20 means per state × condition)
            Results/Tables/Bloc3_GSE207422_CollecTRI_TF_wilcox_bystate.csv (Wilcoxon results)
            Objects/Bloc3_GSE207422_04_tf_acts_raw.rds (raw per-cell TF scores, saved for future re-analysis)

**Superseded:** the original Script 04 (with downsampling and per-cell variance selection, Top 20: MYC, NFKB, JUN, SP1, HSF1, E2F4, E2F1, STAT1, RFXAP, RFXANK, ELK4, RFX5, HSF2, ABL1, SRSF2, TFDP1, NFYC, ZBTB4, CIITA, DOT1L) is retained in version history but no longer used for the main narrative.  

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

#### Script 05c : UCell Hallmark scoring on TAMs (all 50 signatures)
  - Unbiased approach : all 50 Hallmark MSigDB signatures tested
  - Wilcoxon MPR vs non-MPR with BH correction
  - 47 signatures discriminant (p_adj < 0.05)
  - 25 biologically relevant signatures retained after filtering
  - Full results : Results/Tables/GSE243013_UCell_Hallmark_TAMs_wilcox_results.csv
  - Object saved : Objects/GSE243013_seu_TAMs_Hallmark.rds

#### Script 05d : UCell Hallmark TAMs preprint figures
  - Main figure : heatmap 25 signatures (scale=row, p-value asterisks per signature)
  - Biologically irrelevant signatures excluded
  - Figure : Results/Figures/BLOC4A_TAMs/Preprint/GSE243013_Hallmark_TAMs_heatmap_main.png

#### Script 06 : Pseudobulk Hallmark UCell scores MPR vs pCR (patient-level)
- Input: Bloc4A_05_seu_TAMs_UCell.rds
- Signatures: the same 25 Hallmark MSigDB signatures established as discriminant between MPR and non-MPR (not the original custom UCell signatures)
- Subset: MPR (n=10 patients) and pCR (n=11 patients)
  NOTE: non-MPR excluded, focus on MPR vs pCR to test intermediate state hypothesis
- Pseudobulk: mean UCell score per patient per signature
  NOTE: aggregation to patient-level avoids pseudo-replication
- Wilcoxon rank-sum test MPR vs pCR (patient-level), BH-corrected across the 25 signatures
  References: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)

**Superseded:** earlier figures used four custom UCell signatures (M2_immunosuppressive, M1_inflammatory, SPP1_signature, IFN_response) across three groups (non-MPR, MPR, pCR), with pairwise Wilcoxon tests, Bonferroni-corrected. This approach has been replaced to align with the Hallmark-only, BH-corrected methodology used throughout the project.

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

**Superseded:** this Top 20 selection (per-cell variance) and the observations above, including ELK4 in Resident M2, AEBP1, DMTF1, and KAT6B, are superseded by the corrected between-group variance method applied later in this project (see "Preliminary observations : CollecTRI TAMs (Bloc 4-GSE243013)" for the current, statistically verified findings, including the added Wilcoxon test producing Bloc4A_CollecTRI_TF_wilcox_bysubtype.csv).

#### Script 08 : CollecTRI pseudobulk MPR vs pCR, TAMs, GSE243013

- Purpose: extend the non-MPR vs MPR CollecTRI discriminant TF panel (Script 07) to a direct MPR vs pCR comparison, to characterize pCR relative to MPR at the transcription factor level
- Input: Objects/Bloc4A_07_tf_acts_raw.rds (raw per-cell TF scores, already computed)
          Objects/Bloc4A_07_seu_TAMs_TF.rds (metadata: sampleID, pathological_response, tam_short)
- Same Top 20 TFs as Script 07, restricted to MPR and pCR patients
- Patient-level pseudobulk (mean TF score per patient per TAM subtype), Wilcoxon rank-sum MPR vs pCR, BH-corrected within each subtype
- Output: Results/Tables/Bloc4A_08_CollecTRI_pseudobulk_MPR_vs_pCR.csv

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

#### Script 05b : UCell Hallmark scoring on TAMs (all 50 signatures)
  - Unbiased approach : all 50 Hallmark MSigDB signatures tested
  - Wilcoxon MPR vs NMPR with BH correction
  - 43 signatures discriminant (p_adj < 0.05)
  - 23 biologically relevant signatures retained after filtering
  - Full results : Results/Tables/GSE207422_UCell_Hallmark_TAMs_wilcox_results.csv
  - Object saved : Objects/GSE207422_seu_TAMs_Hallmark.rds

#### Script 05c : UCell Hallmark TAMs preprint figures
  - Main figure : heatmap 23 signatures (scale=row, p-value asterisks per signature)
  - Biologically irrelevant signatures excluded
  - Figure : Results/Figures/BLOC4B_Epithelial_TAMs/Preprint/GSE207422_Hallmark_TAMs_heatmap_main.png

#### Script 06 : CollecTRI TF activity on TAMs (GSE207422)
  - Tool: decoupleR run_ulm + CollecTRI network (local CSV)
  - All 9 TAM subtypes included in full, no downsampling (corrected from an earlier version that downsampled to max 10,000 cells/subtype)
  - TF selection method: between-group variance (variance of subtype × response means), corrected from per-cell variance — consistent with GSE243013 methodology
  - Top 20 TFs (between-group variance): RFXAP, RFXANK, RFX5, CIITA, RELA, HSF1, HSF2, NFKB1, NFKBIB, XBP1, REL, CREB1, EGR1, KLF6, IRF1, HIF1A, JUN, JUND, NFKB, IRF5
  - Top 6 (violin plots): RFXAP, RFXANK, HSF1, RFX5, CIITA, RELA
  - Statistical test: Wilcoxon MPR vs NMPR, per TAM subtype, BH-corrected within subtype (Monocyte-derived excluded, 0 MPR cells available)
  - Output: Results/Tables/Bloc4B_GSE207422_CollecTRI_TF_activity.csv (Top 20 means per subtype × condition)
            Results/Tables/Bloc4B_GSE207422_CollecTRI_TF_wilcox_bysubtype.csv (Wilcoxon results)
            Objects/Bloc4B_06_tf_acts_raw.rds (raw per-cell TF scores, saved for future re-analysis)
            2 heatmaps (MPR/NMPR) + key TFs violin
  - Objects saved: Bloc4B_06_seu_TAMs_TF.rds

**Superseded:** the heatmap state-specific observations previously listed here (IFN-stimulated MPR: IRF1/STAT1/IRF5/REL; SPP1+ NMPR: heterogeneous immunosuppressive program; Monocyte-derived NMPR: HIF1A; Stress-response: active in both conditions) were based on a per-cell-variance Top 20 selection and z-score-only visual reading, without formal statistical testing. These are replaced by the statistically verified observations in the "Preliminary observations" section above.

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

#### Script 10 : UCell Hallmark scoring on epithelial cells

- Tool: UCell v2.x (Andreatta & Carmona, Computational and Structural Biotechnology Journal 2021)
- Input: Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds (n=8,944 cells)
- 5 groups: tumor_MPR (n=506), normal_MPR (n=493), tumor_NMPR (n=3,386), normal_NMPR (n=3,574, including reclassified EMT cells), Ciliated (n=565)
- Unbiased approach: all 50 Hallmark MSigDB signatures tested (Liberzon et al., Cell Systems 2015), imported via msigdbr, consistent with CD8/TAM methodology
- **Superseded:** the two custom signatures (HSF1_targets, Antigen_presentation) used in an earlier version of this script are no longer used, replaced by the unbiased Hallmark-only approach
- 3 comparisons tested by Wilcoxon rank-sum, BH-corrected within each: tumor_MPR vs tumor_NMPR; normal_MPR vs normal_NMPR; Ciliated_MPR vs Ciliated_NMPR (Ciliated tested rather than assumed stable, given its original role as a CopyKAT/SCEVAN reference population only)
- Discriminant signatures after biological relevance filtering (12 irrelevant gene sets excluded): 33/50 for tumor_MPR vs NMPR, 33/50 for normal_MPR vs NMPR, 13/50 for Ciliated_MPR vs NMPR
- Output: Results/Tables/Bloc4B_10_UCell_Epithelial_scores.csv (all 50 signatures, per-cell)
          Results/Tables/Bloc4B_10_UCell_Epithelial_wilcox_results.csv (all comparisons, all signatures)
          Results/Tables/Bloc4B_10_UCell_Epithelial_wilcox_discriminant.csv (79 discriminant results)
          Objects/Bloc4B_10_seu_Epithelial_UCell.rds
- Figures: 1 dotplot (all 50 signatures x 5 groups, descriptive) + violin plots for discriminant signatures only

#### Script 11 : CollecTRI TF activity on epithelial cells

- Tool: decoupleR run_ulm (Badia-i-Mompel et al., Bioinformatics Advances 2022)
- Network: CollecTRI (Müller-Dott et al., Nucleic Acids Research 2023), 43,159 interactions, 1,186 TFs
- Input: Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds (n=8,944 cells)
- 5 final groups: tumor_MPR (n=506), normal_MPR (n=493), tumor_NMPR (n=3,386), normal_NMPR (n=3,574), Ciliated (n=565)
- No downsampling applied (manageable dataset size)
- TF selection method: between-group variance (variance of final_group means), corrected from per-cell variance, consistent with CD8/TAM methodology
- Top 20 TFs (between-group variance): RFXANK, RFX5, RFXAP, HSF2, FOXJ1, CREB1, TP53, DMTF1, MYC, EPAS1, ZBTB4, MEIS2, ATF1, ELK4, SP1, NFE2L2, ESR1, HIF1A, NKX2-1, DACH1
- Top 6 TFs (violin plots): RFXANK, RFX5, RFXAP, HSF2, FOXJ1, CREB1
- Statistical test: Wilcoxon, 3 comparisons (tumor_MPR vs NMPR, normal_MPR vs NMPR, Ciliated_MPR vs NMPR), BH-corrected within each comparison
- **Superseded:** the earlier Top 20 (per-cell variance: MYC, RFXAP, RFXANK, HSF2, RFX5, CIITA, NKX2-1, SP1, FOXJ1, DMTF1, NFE2L2, E2F4, HSF1, ELK4, TP53, DACH1, HIF1A, E2F1, NFKB, JUN) and all descriptive (non-statistically-tested) observations built on it are replaced by the between-group variance selection and Wilcoxon results above
- Output: Results/Tables/Bloc4B_11_CollecTRI_Epithelial_TF_activity.csv (Top 20 means per group)
          Results/Tables/Bloc4B_11_CollecTRI_Epithelial_wilcox_results.csv (Wilcoxon results, all 3 comparisons)
          Objects/Bloc4B_11_tf_acts_raw.rds (raw per-cell TF scores, saved for future re-analysis)
          Objects/Bloc4B_11_seu_Epithelial_TF.rds
- Figures: 1 heatmap (descriptive, all 5 groups) + violin plots split tumor/normal (Ciliated included in normal group figure)

## Bloc 5 : TME Intercellular Communication

### Note on tool selection:
Both LIANA+ and CellChat v2 were run and are fully documented below for transparency and cross-validation purposes, across both GSE207422 and GSE243013. LIANA+ was selected as the primary framework for the preprint. Its consensus score is computed as a simple product of mean ligand and receptor expression, which could be recomputed directly from the raw expression matrix on a per-patient basis; this allowed a patient-level pseudobulk Wilcoxon test to be added on top of LIANA's own prioritization score, treating each patient as an independent observation and avoiding pseudo-replication at the single-cell level, consistent with the pseudobulk methodology used throughout this project.

CellChat's communication probability, by contrast, is derived from a mass-action model with Hill-function saturation terms, computed internally by the tool itself rather than as a simple, externally reproducible formula. Recomputing this score on a per-patient basis would require re-running CellChat's full algorithm separately for each patient, computationally demanding, and particularly impractical given the smaller response groups (e.g., n=3 MPR patients in GSE207422) typical of neoadjuvant immunotherapy cohorts. CellChat's native permutation test instead pools all cells within each condition without distinguishing between patients, a design vulnerable to pseudo-replication, where a single patient contributing disproportionately more cells could dominate the result; which does not directly support the patient-level comparison this project's methodology requires. CellChat results are retained here as a complementary, qualitative cross-check but are not reported in the preprint.

### Biological question
Does the transcriptional divergence between MPR and NMPR/non-MPR compartments (CD8, TAMs) reflect distinct intercellular communication programs that reinforce or undermine anti-PD-1 response? Specifically, do responders and non-responders differ in the directionality, specificity, and coordination of bidirectional communication loops between CD8 T cells and TAMs, and does this crosstalk constitute a self-reinforcing immunosuppressive network in non-responders?

### Cells of interest : guided by Blocs 3 and 4 findings (statistically verified, between-group variance + Wilcoxon)

Cell subtypes of interest were selected based on the number of significantly discriminant transcription factors (CollecTRI, Wilcoxon MPR vs non-MPR/NMPR, BH-corrected within each subtype/state, out of the Top 20 tested), using the natural drop-off in each dataset's ranking rather than a fixed cutoff. Full TME was opened for inference; interpretation is prioritized on these axes without excluding unexpected interactions.

**CD8 T cells - both datasets:**
- CD8.TEX and CD8.TPEX: the only two states analyzed, established as the primary exhausted CD8 populations responsive to PD-1 blockade (Liu B et al., 2022). ELK4 is significantly MPR-enriched in both states, in both datasets — the most robust cross-dataset signal identified in this project.

**TAMs - GSE243013 (6 of 7 subtypes retained, natural drop-off after Resident M2):**
LAMs (20/20 significant TFs), Proliferating (20/20), Monocyte FCN1+ (19/20), Resident M2 (19/20), IFN-stimulated (16/20), Stress-response (16/20). Classical-Mono excluded (9/20, the least discriminant subtype).

**TAMs - GSE207422 (5 of 9 subtypes retained, natural drop-off after Resident M2):**
MRC1+ M2-like (19/20 significant TFs), SPP1+ immunosuppressive (19/20), IFN-stimulated (16/20), M2-SIGLEC8+ (12/20), Resident M2 (10/20). Stress-response (1/20), Lipid-associated (0/20) and Regulatory (0/20) excluded. Monocyte-derived not tested (no MPR cells available in this dataset).

**Epithelial: GSE207422 only:**
Excluded from the main preprint narrative (single dataset, no cross-dataset validation possible). Retained in perspectives/supplementary: tumor and normal epithelial compartments each show extensive discriminant TF and Hallmark programs (see Bloc 4B Preliminary observations, Scripts 10-11), with the RFX complex consistently MPR-enriched across tumor, normal, and Ciliated populations.

### Dataset coverage
- GSE207422 : CD8 + TAMs + Epithelial (tumor_MPR, normal_MPR, tumor_NMPR, normal_NMPR, Ciliated) , full bidirectional analysis across all three compartments
- GSE243013 : CD8 + TAMs (no epithelial compartment available) , bidirectional CD8 ↔ TAMs analysis; pCR retained in analysis, excluded from main narrative

### Tools
- LIANA+ v0.1.14 (Dimitrov et al., Nature Communications 2022), consensus LR inference (sca + natmi + connectome)
- CellChat v2 (Jin et al., Nature Communications 2021), differential communication analysis

### Design rationale
Full TME opened for inference: interpretation prioritized on bidirectional CD8 ↔ TAMs, CD8 ↔ Epithelial, TAMs↔Epithelial axes established in Blocs 3 and 4. Both directions of each axis are analyzed to characterize communication loops, not unidirectional signals. Unexpected interactions retained only if they corroborate findings from prior blocs.

### Script 01 : LIANA+ GSE207422, intercellular communication inference

- Input: Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds + Bloc4B_04_seu_TAMs_combined.rds + 
  Bloc4B_09c_seu_Epithelial_FinalGroups.rds
- CD8 restricted to CD8.TEX/CD8.TPEX only (CM, EM, TEMRA, NaiveLike, MAIT excluded, consistent with CD8/TAM CollecTRI methodology)
- Merged TME object: 27,307 cells (CD8 n=5,353 | TAMs n=13,430 | Epithelial n=8,524)
- Methods: sca + natmi + connectome aggregated via liana_aggregate()
- Resource: Consensus (CellChatDB + OmniPath + others)
- Conversion: SingleCellExperiment (SCE) for Seurat v5 compatibility
- Per condition runs: MPR and NMPR (cell counts to confirm via table(seu_TME$cell_type, seu_TME$PathResponse))
- Statistical test: patient-level pseudobulk Wilcoxon, BH-corrected, on interactions prioritized by aggregate_rank <= 0.05 (Dimitrov et al., 2022)
- Output: Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_aggregated.csv
          Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv
          Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_NMPR_aggregated.csv
          Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_differential_wilcox.csv (CD8<->TAM/Epithelial)
          Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_TAM_CD8_Epithelial_wilcox.csv (TAM/CD8<->tumor/normal axis comparison)
          Objects/Bloc5_01_seu_TME_GSE207422.rds

### Script 02 : LIANA+ GSE207422, main axis figures (exploratory, not in preprint)

- Input: Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv + NMPR_aggregated.csv
- 6 bidirectional figures (MPR vs NMPR comparative dotplots), full TME (all 9 TAM subtypes, CD8.TEX/CD8.TPEX only, all epithelial groups)
- Descriptive only: interactions prioritized by aggregate_rank ≤ 0.05 (Dimitrov et al., 2022), NOT statistically tested for MPR vs NMPR difference in this script
- TAM short labels applied for readability (consistent with Blocs 4A/4B)
- Epithelial labels simplified: tumor_MPR + tumor_NMPR → Tumor epithelial; normal_MPR + normal_NMPR → Normal epithelial; Ciliated kept separate
- Top 5 interactions per target group (aggregate_rank ≤ 0.05)
- Output:
  1. CD8 → TAMs (Bloc5_02_LIANA_CD8_TAMs.png)
  2. TAMs → CD8 (Bloc5_02_LIANA_TAMs_CD8.png)
  3. CD8 → Epithelial (Bloc5_02_LIANA_CD8_Epithelial.png)
  4. Epithelial → CD8 (Bloc5_02_LIANA_Epithelial_CD8.png)
  5. TAMs → Epithelial (Bloc5_02_LIANA_TAMs_Epithelial.png)
  6. Epithelial → TAMs (Bloc5_02_LIANA_Epithelial_TAMs.png)
- **Not used in the preprint main text or supplementary** — kept on disk for reference only. Formal Wilcoxon statistical testing (patient-level pseudobulk, BH-corrected), consistent with the rest of the project, is applied in Script 02b instead.

### Script 02b : LIANA+ GSE207422, preprint figures with cells of interest

- Input: Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv + NMPR_aggregated.csv + Objects/Bloc5_01_seu_TME_GSE207422.rds (for statistical testing)
- Same 6 bidirectional axes as Script 02, filtered on cells of interest guided by Blocs 3 and 4 CollecTRI findings
- CD8 priority: CD8.TEX, CD8.TPEX (restricted at Script 01, merged object level)
- TAMs priority (based on CollecTRI significant TF count, natural drop-off after Resident M2): MRC1+ M2-like, SPP1+ immunosuppressive, IFN-stimulated, M2-SIGLEC8+, Resident M2
- Epithelial: Tumor epithelial + Normal epithelial + Ciliated (all kept)
- Top 5 interactions per target group prioritized by aggregate_rank ≤ 0.05 (descriptive), each additionally tested for MPR vs NMPR difference by patient-level pseudobulk Wilcoxon, BH-corrected (statistical, not just prioritization)
- On figures: non-significant interactions (Wilcoxon p_adj ≥ 0.05) shown with reduced opacity, to visually distinguish LIANA prioritization from statistical significance
- Output:
  1. Bloc5_02b_LIANA_CD8_TAMs_preprint.png + Bloc5_02b_LIANA_GSE207422_wilcox_CD8_to_TAMs.csv
  2. Bloc5_02b_LIANA_TAMs_CD8_preprint.png + Bloc5_02b_LIANA_GSE207422_wilcox_TAMs_to_CD8.csv
  3. Bloc5_02b_LIANA_CD8_Epithelial_preprint.png + Bloc5_02b_LIANA_GSE207422_wilcox_CD8_to_Epithelial.csv
  4. Bloc5_02b_LIANA_Epithelial_CD8_preprint.png + Bloc5_02b_LIANA_GSE207422_wilcox_Epithelial_to_CD8.csv
  5. Bloc5_02b_LIANA_TAMs_Epithelial_preprint.png + Bloc5_02b_LIANA_GSE207422_wilcox_TAMs_to_Epithelial.csv
  6. Bloc5_02b_LIANA_Epithelial_TAMs_preprint.png + Bloc5_02b_LIANA_GSE207422_wilcox_Epithelial_to_TAMs.csv

### Script 03 : CellChat v2 GSE207422, differential communication analysis

- Input: Bloc5_01_seu_TME_GSE207422.rds
- CD8 restricted to CD8.TEX/CD8.TPEX only (consistent with Bloc 5 methodology)
- CellChatDB.human full database
- Communication probability computed via trimean method; significance via built-in permutation test (nboot = 100)
- Per condition objects: cellchat_MPR.rds + cellchat_NMPR.rds (each computed independently on cells from that condition only)
- min.cells = 5
- Merged object: Bloc5_03_cellchat_merged_GSE207422.rds
- Output: Objects/Bloc5_03_cellchat_MPR.rds
          Objects/Bloc5_03_cellchat_NMPR.rds
          Objects/Bloc5_03_cellchat_merged_GSE207422.rds
          Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_MPR_interactions.csv
          Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_NMPR_interactions.csv

### Script 04 : CellChat v2 GSE207422, custom figures (exploratory)

- Input: Bloc5_03_CellChat_GSE207422_MPR_interactions.csv + NMPR_interactions.csv
- Custom ggplot2 figures built from exported CellChat CSV (not netVisual_bubble)
- 6 bidirectional figures, full TME:
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
- TAMs priority (based on CollecTRI significant TF count, natural drop-off after Resident M2): MRC1+ M2-like, SPP1+ immunosuppressive, IFN-stimulated, M2-SIGLEC8+, Resident M2
- Epithelial: Tumor epithelial + Normal epithelial + Ciliated (all kept)
- MPR vs NMPR crossing for all 3 axes (CD8<->TAM, CD8<->Epithelial, TAM<->Epithelial): interactions matched on (source, target, ligand, receptor); presence flagged as Both/MPR_only/NMPR_only. For Tumor/Normal epithelial (condition encoded in cell type name), the MPR-arm and NMPR-arm of the same axis are compared instead of a direct match
- Only MPR_only/NMPR_only interactions retained (both rest on CellChat's own native permutation test, applied independently within each condition); "Both" excluded, since CellChat provides no native test to compare an interaction's strength directly between conditions
- Output: Results/Tables/BLOC5/Bloc5_04b_CellChat_GSE207422_CD8_TAM_crossed.csv
          Results/Tables/BLOC5/Bloc5_04b_CellChat_GSE207422_CD8_Epithelial_crossed.csv
          Results/Tables/BLOC5/Bloc5_04b_CellChat_GSE207422_TAM_Epithelial_crossed.csv
- Preprint figures (CD8 ↔ TAMs, main narrative):
  1. CD8 → TAMs (Bloc5_04b_CellChat_CD8_TAMs_preprint.png)
  2. TAMs → CD8 (Bloc5_04b_CellChat_TAMs_CD8_preprint.png)
- Discussion figures (epithelial axes — GSE207422 only):
  3. CD8 → Epithelial (Bloc5_04b_CellChat_CD8_Epithelial_discussion.png)
  4. Epithelial → CD8 (Bloc5_04b_CellChat_Epithelial_CD8_discussion.png)
  5. TAMs → Epithelial (Bloc5_04b_CellChat_TAMs_Epithelial_discussion.png)
  6. Epithelial → TAMs (Bloc5_04b_CellChat_Epithelial_TAMs_discussion.png)
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

#### Script 09 : Individual HLA-D gene expression by TAM subtype, GSE207422

- Purpose: resolve an apparent discrepancy between CollecTRI's aggregate RFX/CIITA regulon activity and the direction of the HLA-D-LAG3 ligand-receptor communication signal, by testing individual HLA-D gene expression rather than the combined ligand-receptor product score
- Input: Objects/Bloc4B_04_seu_TAMs_combined.rds
- Restricted to the 3 TAM subtypes engaging LAG3 in communication analysis: M2-SIGLEC8+, MRC1+ M2-like, SPP1+ immunosuppressive
- Genes tested: HLA-DQB1, HLA-DRB3
- Patient-level pseudobulk (mean expression per patient per subtype), Wilcoxon rank-sum MPR vs NMPR, BH-corrected
- Output: Results/Tables/BLOC5/Bloc5_09_HLA_D_individual_genes_wilcox.csv

#### Script 09b : LAG3 expression, CD8 side, GSE207422

- Purpose: companion test to Script 09, verifying whether LAG3 itself (CD8 side of the HLA-D-LAG3 axis) shows a significant standalone MPR vs NMPR difference
- Input: Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds, restricted to CD8.TEX/CD8.TPEX
- Gene tested: LAG3
- Patient-level pseudobulk, Wilcoxon rank-sum MPR vs NMPR, BH-corrected, tested separately for CD8.TEX and CD8.TPEX
- Output: Results/Tables/BLOC5/Bloc5_09b_LAG3_CD8_wilcox.csv


---

# Preliminary Observations: Reorganized by Compartment

---

# CD8 T CELLS

## Preliminary observations : TME composition and CD8 states (Bloc 2-3)

- CD8.TEX gradient: non-MPR > MPR > pCR (visual UMAP), confirmed statistically by ProjecTILs barplot (Chi-2 p<2.2e-16)
- non-MPR dominance (n=42) requires proportional analysis, absolute counts not comparable
- GSE207422 preliminary analysis identified CD8_Exhausted_Terminal enriched in MPR (OR=3.36, different annotation tool). Current harmonized analysis applies identical tools (ProjecTILs, UCell, CollecTRI) to both datasets for direct cross-dataset comparison.

## Preliminary observations : ProjecTILs, CD8 T cell states (Bloc 3, GSE243013)

- CD8.TEX gradient confirmed: non-MPR (32.3%) > MPR (20.2%); pCR (11.5%) retained for reference
- CD8.CM inverse pattern: MPR (35.1%) > non-MPR (25.9%)
- CD8.TPEX enriched in pCR, consistent with reactivation of precursor-exhausted cells
- CD8.EM similar across groups, functional status assessed by UCell (Script 10)
- Proposed gradient: non-MPR (TEX dominant) → MPR (TEX + TPEX) → pCR (TPEX dominant)
Hypothesis: anti-PD-1 reactivated a fraction of TEX toward TPEX in responders

## Preliminary observations : UCell custom signatures, CD8 (Bloc 3, GSE243013)

- Exhaustion gradient confirmed: non-MPR (0.347) > pCR (0.249) > MPR (0.220) in CD8.TEX
- CD8.EM exhaustion: non-MPR (0.128) > MPR (0.107) ≈ pCR (0.106), EM dysfunction in non-responders
- CD8.TPEX retains higher cytotoxicity than CD8.TEX across all groups, preserved precursor function
- Memory score highest in CD8.NaiveLike = expected

**Hypotheses:**
- H1 (CollecTRI): superior intrinsic CD8 reactivation capacity in responders
- H2 (CellChat/NicheNet): immunosuppressive TME determines reactivation failure
- H3: interaction between CD8 intrinsic state and TME context

## Preliminary observations : UCell Hallmark CD8 (Bloc 3, GSE243013)

Unbiased scoring of all 50 Hallmark MSigDB signatures on CD8.TEX and CD8.TPEX (MPR vs non-MPR). Discriminant signatures identified by Wilcoxon test, BH-corrected separately within each CD8 state (not on the pooled TEX+TPEX population), to avoid one state's signal masking or inflating the other's.

**CD8.TEX - non-MPR enriched (significant):**
FATTY_ACID_METABOLISM, GLYCOLYSIS, G2M_CHECKPOINT, OXIDATIVE_PHOSPHORYLATION, INFLAMMATORY_RESPONSE, COMPLEMENT, APOPTOSIS, E2F_TARGETS, MYC_TARGETS_V2, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE, UNFOLDED_PROTEIN_RESPONSE, MTORC1_SIGNALING, consistent with broad metabolic, proliferative and inflammatory dysregulation.

**CD8.TEX - MPR enriched (significant):**
HYPOXIA, KRAS_SIGNALING_UP, P53_PATHWAY, TNFA_SIGNALING_VIA_NFKB, EPITHELIAL_MESENCHYMAL_TRANSITION, TGF_BETA_SIGNALING.

**CD8.TPEX - non-MPR enriched (significant):**
FATTY_ACID_METABOLISM, COMPLEMENT, E2F_TARGETS, MYC_TARGETS_V2, INTERFERON_ALPHA_RESPONSE, INTERFERON_GAMMA_RESPONSE, UNFOLDED_PROTEIN_RESPONSE, MTORC1_SIGNALING.

**CD8.TPEX - MPR enriched (significant):**
HYPOXIA, TNFA_SIGNALING_VIA_NFKB, EPITHELIAL_MESENCHYMAL_TRANSITION.

**CD8.TPEX - not significant despite pooled-test significance (state-specific dilution):**
INFLAMMATORY_RESPONSE, APOPTOSIS, KRAS_SIGNALING_UP, P53_PATHWAY, GLYCOLYSIS, OXIDATIVE_PHOSPHORYLATION, G2M_CHECKPOINT, TGF_BETA_SIGNALING. These reached significance in the pooled TEX+TPEX test (and in TEX alone for most), but not when tested within TPEX specifically, a reminder that the global test can be driven by one state alone.

**Overall pattern:**
HYPOXIA, TNFA_SIGNALING_VIA_NFKB, and EPITHELIAL_MESENCHYMAL_TRANSITION are the only signatures confirmed MPR-enriched in both CD8.TEX and CD8.TPEX, a more modest but state-consistent MPR signal than previously described. Non-MPR shows a broader, largely shared program across both states (FATTY_ACID_METABOLISM, COMPLEMENT, E2F_TARGETS, MYC_TARGETS_V2, INTERFERON_ALPHA/GAMMA_RESPONSE, UNFOLDED_PROTEIN_RESPONSE, MTORC1_SIGNALING), consistent with broad transcriptional dysregulation in non-responders across the exhausted CD8 lineage, rather than two fully distinct state-specific programs as previously suggested.

## Preliminary observations : TCR analysis (Bloc 3, GSE243013)

- CD8.TEX 85.4% expanded, confirms tumor-reactive identity
- Clonal diversity: MPR (0.490) > pCR (0.458) > non-MPR (0.438), responders maintain broader TCR repertoire
- pCR more polyclonal, broader anti-tumor response

## Preliminary observations : ProjecTILs, CD8 T cell states (Bloc 3, GSE207422)

**ProjecTILs:**
- 7 CD8 states identified, consistent with GSE243013
- CD8.TEX dominant MPR 51.2% vs NMPR 47.0%, less pronounced gradient than GSE243013
- CD8.TEMRA enriched NMPR (4.0% vs 2.2%), consistent with terminal differentiation
- CD8.TPEX minor population both conditions (MPR 1.3%, NMPR 0.9%)

*Note: proportions consistent with Bloc3_GSE207422_Barplot_CD8_ProjecTILs.png.*

## Preliminary observations : UCell Hallmark CD8 (Bloc 3, GSE207422)

Unbiased scoring of all 50 Hallmark MSigDB signatures on CD8.TEX and CD8.TPEX (MPR vs NMPR). 22 biologically relevant discriminant signatures retained after filtering (p_adj < 0.05 in at least one CD8 state, Wilcoxon BH-corrected within each state).

**Methodological note on cohort imbalance:** MPR (n=3) vs NMPR (n=10) in this dataset is severely imbalanced. With unequal group sizes, statistical power is asymmetric and can inflate apparent significance on the larger group's side; the pattern below should be interpreted with this caveat and cross-checked against GSE243013 (more balanced, MPR n=10 vs non-MPR n=40) before treating it as a robust biological signal.

**CD8.TEX, near-uniform NMPR enrichment:**
All 22 tested signatures are significantly enriched in NMPR relative to MPR (all p_adj < 0.05), with no exception. This includes metabolic (OXIDATIVE_PHOSPHORYLATION, GLYCOLYSIS, FATTY_ACID_METABOLISM), proliferative (MYC_TARGETS_V1/V2, E2F_TARGETS, G2M_CHECKPOINT, MTORC1_SIGNALING), inflammatory (TNFA_SIGNALING_VIA_NFKB, INFLAMMATORY_RESPONSE, IL6_JAK_STAT3_SIGNALING, INTERFERON_ALPHA/GAMMA_RESPONSE, ALLOGRAFT_REJECTION), and stress/other programs (P53_PATHWAY, APOPTOSIS, UNFOLDED_PROTEIN_RESPONSE, HYPOXIA, KRAS_SIGNALING_UP, PI3K_AKT_MTOR_SIGNALING, IL2_STAT5_SIGNALING, EPITHELIAL_MESENCHYMAL_TRANSITION).

**CD8.TPEX, largely consistent with TEX, weaker statistical power:**
17/22 signatures significant, all direction NMPR > MPR (no signature reaches significance in the MPR direction). Significant: G2M_CHECKPOINT, P53_PATHWAY, MYC_TARGETS_V2, PI3K_AKT_MTOR_SIGNALING, E2F_TARGETS, KRAS_SIGNALING_UP, APOPTOSIS, TNFA_SIGNALING_VIA_NFKB, GLYCOLYSIS, OXIDATIVE_PHOSPHORYLATION, MYC_TARGETS_V1, HYPOXIA, INFLAMMATORY_RESPONSE, MTORC1_SIGNALING, EPITHELIAL_MESENCHYMAL_TRANSITION. Not significant (but same directional trend): ALLOGRAFT_REJECTION, UNFOLDED_PROTEIN_RESPONSE, IL2_STAT5_SIGNALING, INTERFERON_GAMMA_RESPONSE, INTERFERON_ALPHA_RESPONSE, FATTY_ACID_METABOLISM. The single exception is IL6_JAK_STAT3_SIGNALING, which trends toward MPR but does not reach significance (p_adj = 0.15).

**Overall pattern:**
Contrary to a nuanced state-specific split, both CD8.TEX and CD8.TPEX show a broadly NMPR-dominant transcriptional landscape across nearly all tested Hallmark programs in this dataset; no signature shows robust, significant MPR enrichment in either state. This near-uniform direction, combined with the severe MPR/NMPR imbalance (n=3 vs n=10), warrants caution: it may reflect a genuine biological signal, a power asymmetry artifact, or both. Cross-dataset comparison with GSE243013 is needed before drawing conclusions about CD8 Hallmark-level MPR/NMPR divergence.

## Preliminary observations : CollecTRI CD8 (Bloc 3, both datasets)

**Methodological correction:** Initial violin plots used manually selected TFs based on prior DoRothEA hypotheses, confirmation bias risk. Corrected to Top 6 by cross-cell variance (objective). Top 10 and Top 15 tested, conclusions unchanged.

**Second methodological correction:** TF selection (Top 6 and Top 20) was further corrected from per-cell variance to between-group variance (variance of state × response means), consistent with comparing response groups rather than generic feature selection. This changed the composition of the Top 20 substantially in both datasets, see below. Applied uniformly to GSE243013 and GSE207422.

**Visualization approach:** Heatmaps row-scaled (z-score) show state-specific TF patterns within each CD8 state per condition. Violin plots (Top 6 by variance, absolute scores) show global MPR vs non-MPR differences. The two are complementary.

### GSE243013 (Script 12b, no downsampling, CD8.TEX and CD8.TPEX)

**Two complementary visualizations:**
- **Violin plots (Top 6 by between-group variance, absolute scores)**: global MPR vs non-MPR comparison across CD8.TEX and CD8.TPEX
- **Heatmaps (scale=row, Top 20 by between-group variance)**: state-specific TF patterns within each CD8 state per condition

**Violin plots, Top 6 by between-group variance (GSE243013):**
HSF1, HSF2, NFYC, MLXIP, FOXO3, HOPX. All six show the same directional pattern: lowest in MPR, intermediate in non-MPR, highest in pCR (FOXO3 being the exception, with MPR ≈ non-MPR and only pCR clearly elevated). This Top 6 is distinct from ELK4/MTF1, which, despite showing the most statistically significant MPR-enrichment in the Wilcoxon test (see below), do not rank among the six TFs with the highest variance across all six state × condition groups (including pCR).

*Note: cross-dataset consensus for the violin Top 6 cannot yet be confirmed, GSE207422's Top 6 (below) uses a different TF list under the same corrected method.*

**Heatmap observations, statistically verified (Wilcoxon MPR vs non-MPR, BH-corrected within each state; supersedes all earlier observations, both the downsampled Script 12 and the initial no-downsampling per-cell-variance Top 20):**

New Top 20 by between-group variance: MTF1, ELK4, HOPX, HSF4, HSF2, FOXA1, HSF1, NFAT5, NFYB, NFYC, FOXF2, KLF4, MYC, MLX, RFX5, RFXANK, RFXAP, MLXIP, DAXX, FOXO3.

*Note: RLF, ZBTB4, CIITA, TBX21, E2F1 and E2F4, all featured in earlier observations, are no longer part of the Top 20 under this corrected selection method and are not discussed below.*

- **CD8.TEX MPR-enriched (significant):** MTF1, ELK4, MLX
- **CD8.TEX non-MPR-enriched (significant):** HOPX, HSF4, HSF2, FOXA1, HSF1, NFAT5, NFYB, NFYC, FOXF2, KLF4, MYC, RFX5, RFXANK, RFXAP
- **CD8.TEX not significant:** MLXIP, DAXX, FOXO3
- **CD8.TPEX MPR-enriched (significant):** MTF1, ELK4
- **CD8.TPEX non-MPR-enriched (significant):** FOXA1, HSF2, NFAT5, NFYB, HSF4, NFYC, HSF1, FOXF2, HOPX, KLF4, MYC, DAXX, MLXIP, RFX5
- **CD8.TPEX not significant:** RFXANK, RFXAP, MLX, FOXO3

**Key points:**
- ELK4 remains the most robust MPR-enriched TF, now with formal statistical confirmation in both states (CD8.TEX: p_adj ≈ 0; CD8.TPEX: p_adj = 0.021), not just a descriptive z-score elevation as before.
- MTF1 emerges as a new, strongly significant MPR-enriched TF in both states (CD8.TEX: p_adj ≈ 0; CD8.TPEX: p_adj = 3.7×10⁻¹⁵), not previously identified; its role has not yet been characterized in this context and warrants literature review before being discussed in the Discussion section.
- The MHC II signature in non-MPR CD8.TEX is now supported by all three RFX subunits (RFX5, RFXANK, RFXAP), but **not by CIITA**, which ranked 42nd out of 752 TFs by between-group variance (same direction as RFX, non-MPR > MPR, but well below the Top 20 cutoff). Since CIITA is the non-DNA-binding coactivator required for RFX-bound enhanceosomes to drive MHC II transcription (RFX alone cannot transactivate MHC II; Gobin et al., 1998; Masternak et al., 2000), this may indicate an incomplete or only partially assembled MHC II program rather than a fully functional antigen-presentation shift.
- In CD8.TPEX, only RFX5 reaches significance (non-MPR); RFXANK and RFXAP fall just short (p_adj = 0.076 and 0.080), the MHC II signal is weaker and less complete in TPEX than in TEX.
- pCR and CD8.EM: not included in this analysis (Script 12b restricted to CD8.TEX/CD8.TPEX, MPR vs non-MPR); pCR-specific observations are retained separately in the perspectives/limitations section only.

### GSE207422 (Script 04b, no downsampling, CD8.TEX and CD8.TPEX)

**Top 20 by between-group variance:** MYC, RFXAP, RFXANK, RFX5, SP1, NFKB, ELK4, HIF1A, E2F4, IRF6, SRSF2, TP53, HOPX, E2F1, PAWR, NAB2, STAT3, NFKB1, JUN, TFDP1.
**Top 6 (violin plots):** MYC, RFXAP, RFXANK, RFX5, SP1, NFKB.

*Note: RLF, TBX21, NFYC, ABL1, DOT1L, ZBTB4, STAT1 and CIITA, all featured in earlier observations from the downsampled per-cell-variance Script 04, are no longer part of the Top 20 under this corrected selection method and are not discussed below.*

**Heatmap observations, statistically verified (Wilcoxon MPR vs NMPR, BH-corrected within each state):**

- **CD8.TEX MPR-enriched (significant):** ELK4, PAWR
- **CD8.TEX NMPR-enriched (significant):** RFXAP, HOPX, RFX5, HIF1A, RFXANK, NFKB, TP53, SP1, MYC, NFKB1, NAB2, IRF6, JUN, TFDP1, STAT3, SRSF2, E2F4, E2F1 (all remaining Top 20 TFs)
- **CD8.TPEX MPR-enriched (significant):** PAWR, ELK4
- **CD8.TPEX NMPR-enriched (significant):** SP1, HIF1A, STAT3, MYC, NFKB, IRF6, TP53, JUN, SRSF2, TFDP1, E2F1, NFKB1, HOPX, E2F4
- **CD8.TPEX not significant:** RFXANK, RFXAP, RFX5, NAB2

**Key points:**
- ELK4 significantly MPR-enriched in both states (CD8.TEX: p_adj ≈ 0; CD8.TPEX: p_adj = 7.1×10⁻³), confirming the cross-dataset signal.
- PAWR, a pro-apoptotic tumor suppressor not identified in GSE243013's Top 20, emerges as an additional significant MPR-enriched TF in both states (CD8.TEX: p_adj ≈ 0; CD8.TPEX: p_adj = 2.2×10⁻⁵).
- The RFX complex (RFX5/RFXANK/RFXAP) is significant in CD8.TEX (all NMPR-enriched) but **none of the three subunits reach significance in CD8.TPEX** (all p_adj > 0.10), the same TEX-only MHC-II pattern seen in GSE243013, here even more pronounced.
- CIITA ranked 22nd out of ~752 TFs by between-group variance, narrowly missing the Top 20. In CD8.TEX, same direction as RFX (higher NMPR: 12.4 vs 11.0), consistent with an incomplete MHC II program. In CD8.TPEX, opposite direction (higher MPR: 11.2 vs 10.8), diverging from the (non-significant) RFX trend in this state, unlike GSE243013 where CIITA's direction matched RFX in both states.

### Cross-dataset conclusion (statistically verified, between-group variance + Wilcoxon by state)

- **ELK4 is the most robust MPR-associated transcription factor**, significantly enriched in both CD8.TEX and CD8.TPEX, in both datasets, the only TF with this level of cross-state, cross-cohort consistency.
- **The RFX complex is a consistent non-MPR/NMPR signal restricted mainly to CD8.TEX.** In CD8.TPEX, the RFX signal is markedly weaker: only RFX5 reaches significance in GSE243013 (p_adj = 7.7×10⁻³); none of the three subunits reach significance in GSE207422. This TEX-specific pattern replicates across both cohorts.
- **CIITA is absent from the Top 20 in both datasets** (ranked 42nd in GSE243013, 22nd in GSE207422), despite RFX being a strong discriminant signal. Since CIITA is required as the non-DNA-binding coactivator for RFX-bound enhanceosomes to transactivate MHC II genes, this points to a possibly incomplete or only partially assembled MHC II program in non-MPR/NMPR CD8.TEX, rather than a fully functional antigen-presentation shift. This interpretation holds more consistently in GSE243013 (CIITA's direction matches RFX in both states) than in GSE207422, where CIITA diverges from RFX specifically in CD8.TPEX.
- **ELK4's co-elevated partner TF differs between datasets**: MTF1 in GSE243013 (both states, highly significant, uncharacterized in this context) versus PAWR in GSE207422 (both states, highly significant, a known pro-apoptotic tumor suppressor). Neither partner was tested in the other dataset (absent from its Top 20), so cross-dataset consensus on this specific companion signal cannot yet be established, treated as an open question, not a discrepancy to resolve.
- **Six transcription factors were independently selected into the Top 20 of both datasets** despite differing cohort composition: ELK4, HOPX, MYC, RFX5, RFXANK, RFXAP.
- non-MPR/NMPR CD8.TPEX shows a chronic proliferative/stress program (MYC, E2F-family, HSF-family members, depending on dataset) without cytotoxic coordination, confirmed cross-dataset, though the exact TF composition differs (HSF4/HSF2/FOXA1/NFAT5/NFYB/NFYC/FOXF2/KLF4/DAXX/MLXIP in GSE243013 vs SP1/HIF1A/STAT3/MYC/NFKB/IRF6/TP53/JUN/SRSF2/TFDP1/E2F1/NFKB1/HOPX/E2F4 in GSE207422).

*Full interpretive discussion of these findings (ELK4/MTF1/PAWR, CIITA/RFX incompleteness) is developed in the preprint, Discussion Section 1.*

## Preliminary observations : Pseudobulk CD8 MPR vs pCR (Bloc 3, Script 13, GSE243013)

**Methodological note:** earlier version used four custom UCell signatures (Cytotoxicity, Exhaustion, Memory, TPEX), all non-significant. Corrected to reuse the same 21 Hallmark MSigDB signatures established as discriminant between MPR and non-MPR, tested here between MPR and pCR (patient-level pseudobulk, BH-corrected).

No signature reached significance after BH correction (all p_adj > 0.05). Three signatures approached the threshold without crossing it: OXIDATIVE_PHOSPHORYLATION, INFLAMMATORY_RESPONSE and UNFOLDED_PROTEIN_RESPONSE (all p_adj = 0.056), a weak trend rather than a null result, in contrast with TAMs, where all 25 signatures were far from significance (all p_adj > 0.89, Bloc4A Script 06). This raises the possibility that CD8 T cell state, rather than TAM composition alone, may carry more of the residual biological signal distinguishing MPR from pCR, though this interpretation requires caution given the small sample size (n=10 vs n=11) and the fact that no signature reached formal significance. This weak CD8-level trend is examined alongside the TAM compartment below.

**Perspectives:**
- TAM-CD8 interaction analysis (CellChat/NicheNet) may help resolve whether the immunosuppressive TAM context differentially constrains CD8 function in MPR vs pCR; planned in subsequent blocs
- Larger balanced cohort with sufficient statistical power would be required to formally test CD8-level differences between MPR and pCR at patient level
- Spatial transcriptomics could resolve TAM/CD8 co-localization patterns in the TME to complement these transcriptomic observations

## Preliminary observations : CollecTRI pseudobulk MPR vs pCR, CD8 (Bloc 3, Script 14, GSE243013)

**Methodological note:** this analysis extends the existing UCell-based MPR vs pCR pseudobulk comparison to the transcription factor level, testing the same 20 CollecTRI TFs already established as discriminant between non-MPR and MPR (Script 12b) between MPR and pCR specifically. Patient-level pseudobulk, Wilcoxon rank-sum, BH-corrected within each CD8 state.

**Result:** 11 of 40 tests (20 TFs × 2 states) reached significance. In CD8.TEX, 9 TFs were significant, all higher in pCR than MPR (HOPX, KLF4, NFAT5, FOXA1, HSF1, HSF2, HSF4, NFYC, FOXO3). In CD8.TPEX, only 2 TFs were significant: KLF4 (higher in pCR) and MTF1 (higher in MPR, opposite direction).

**Cross-referencing with the non-MPR vs MPR comparison (Script 12b):** of the 9 CD8.TEX-significant TFs here, 8 (HOPX, KLF4, NFAT5, FOXA1, HSF1, HSF2, HSF4, NFYC) were also independently significant and non-MPR-enriched in the earlier non-MPR vs MPR comparison. For these 8 TFs, MPR shows the lowest activity of the three groups (non-MPR > MPR and pCR > MPR), rather than sitting at an intermediate point on a monotonic gradient. This pattern is specific to CD8.TEX; CD8.TPEX shows no comparable pattern, with only KLF4 shared and MTF1 showing the opposite relationship (MPR > pCR, consistent with MTF1's earlier MPR-enriched status vs non-MPR). FOXO3 was significant only in this MPR vs pCR comparison and was not part of the original non-MPR vs MPR discriminant panel, so its full three-group behavior is unconfirmed.

**Interpretation, held provisional:** this analysis does not establish whether non-MPR and pCR are similar to or different from each other; only their relationship to MPR was tested. A direct non-MPR vs pCR comparison was not performed here given the substantial imbalance between these groups (non-MPR n=42 vs pCR n=11), which would limit the reliability of such a test in this already exploratory three-way comparison. Whether MPR represents a genuine outlier state relative to both non-MPR and pCR, or whether this pattern reflects some other structure in the data, remains an open question, to be examined further once communication-level analysis (LIANA+) is extended to the MPR vs pCR comparison for this dataset.

**Output:** Results/Tables/Bloc3_14_CollecTRI_pseudobulk_MPR_vs_pCR.csv

---

# TUMOR-ASSOCIATED MACROPHAGES (TAMs)

## TAM subtype annotation reference (GSE207422)

## TAM subtype annotation reference

### GSE207422

| Short name (used throughout) | Original annotation | Marker basis |
|---|---|---|
| MRC1+ M2-like | TAM_like_MRC1 | MRC1 (CD206) |
| SPP1+ immunosuppressive | TAM_like_SPP1 | SPP1 (osteopontin) |
| Resident M2 | TAM_like_resident_M2 (iron metabolism/anti-inflammatory) | Iron metabolism, anti-inflammatory program |
| IFN-stimulated | TAM_like_IFN (PD-L1+/IDO1+/CXCL9+) | Hybrid: CXCL9 (pro-inflammatory/M1) + PD-L1, IDO1 (immunosuppressive) |
| Monocyte-derived | TAM_like_monocyte (classical inflammatory) | Classical monocyte-derived, inflammatory |
| Lipid-associated | TAM_like_lipid (CCL18+/AKR+) | CCL18 (M2 marker), AKR (lipid metabolism) |
| Stress-response | TAM_like_stress (HSP-high/M1-like) | Heat shock proteins, M1-like |
| Regulatory | TAM_like_regulatory (glucocorticoid-responsive) | Glucocorticoid-responsive program |
| M2-SIGLEC8+ | TAM_like_M2 (SIGLEC8+/CCL18+) | SIGLEC8, CCL18 (M2 marker) |

### GSE243013

| Short name (tam_short) | Original annotation | Marker basis |
|---|---|---|
| Resident M2 | Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like) | Anti-inflammatory, M2-like program |
| LAMs | TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs) | TREM2, APOE (lipid-associated TAM markers) |
| Monocyte FCN1+ | Inflammatory monocyte-derived TAMs (FCN1+/S100A8+) | FCN1, S100A8 (inflammatory monocyte markers) |
| Stress-response | Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high) | MARCO, PPARG, heat shock proteins |
| Proliferating | Proliferating TAMs (cycling/MKI67+) | MKI67 (proliferation marker) |
| IFN-stimulated | IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+) | Interferon-stimulated genes, PD-L1, IDO1 (immunomodulatory) |
| Classical-Mono | Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+) | S100A8, S100A9, CCR2 (classical monocyte markers) |

*Note: "M1/M2" labels reflect original marker-based annotation and should be read as descriptive shorthand rather than strict functional classification. Subtype names differ between the two datasets (no universal TAM nomenclature exists in the field; Zhu et al., 2023), with three exceptions carrying the same short label in both cohorts: Resident M2, Stress-response, and IFN-stimulated, though their full marker basis differs slightly (see above), so direct cross-dataset equivalence should not be assumed without checking marker overlap.*

## Preliminary observations : TAMs UCell custom signatures (Bloc 4, GSE207422)

**UCell scores (9 subtypes combined):**
- M2_immunosuppressive: NMPR > MPR, p<0.0001
- SPP1_signature: NMPR > MPR, p<0.0001
- IFN_response: MPR > NMPR, p=3e-04
- M1_inflammatory: ns (p=0.2848), limited power (MPR n=3)

**Key subtype enrichments:**
- NMPR: TAM_like_lipid, monocyte, regulatory, SPP1, stress-immunosuppressive niche
- MPR: TAM_like_IFN, dual role hypothesis (CXCL9/10 pro-immunogenic vs suppressive)
- TAM_like_resident_M2 MPR > NMPR, discordant with GSE243013, may reflect NSCLC vs LUAD composition differences

**Biological interpretation:**
- NMPR immunosuppressive TAM niche confirmed, supports H2
- MPR as biologically unstable intermediate, two divergent trajectories:
  1. Relapse: residual immunosuppressive TAMs + progressive CD8 exhaustion
  2. Deepening toward pCR: pro-immunogenic TAMs + CD8 TPEX plasticity

## Preliminary observations : UCell Hallmark TAMs (Bloc 4B, GSE207422)

Unbiased scoring of all 50 Hallmark MSigDB signatures on TAMs (MPR vs NMPR). 23 biologically relevant discriminant signatures retained after filtering (p_adj < 0.05, Wilcoxon rank-sum test, BH-corrected).

**NMPR, enriched signatures:**
HYPOXIA, ANGIOGENESIS, KRAS_SIGNALING_UP, GLYCOLYSIS, EPITHELIAL_MESENCHYMAL_TRANSITION, TNFA_SIGNALING_VIA_NFKB, TGF_BETA_SIGNALING, INFLAMMATORY_RESPONSE, IL2_STAT5_SIGNALING, P53_PATHWAY, G2M_CHECKPOINT, UNFOLDED_PROTEIN_RESPONSE, IL6_JAK_STAT3_SIGNALING, MYC_TARGETS_V2.

**MPR, enriched signatures:**
COMPLEMENT, REACTIVE_OXYGEN_SPECIES_PATHWAY, MYC_TARGETS_V1, FATTY_ACID_METABOLISM, OXIDATIVE_PHOSPHORYLATION, INTERFERON_ALPHA_RESPONSE, ALLOGRAFT_REJECTION, PI3K_AKT_MTOR_SIGNALING, NOTCH_SIGNALING.

**Inter-dataset asymmetries:**
HYPOXIA, enriched non-MPR GSE243013 but enriched MPR GSE207422.
MYC_TARGETS_V2, enriched non-MPR GSE243013 but enriched MPR GSE207422.
MYC_TARGETS_V1, enriched non-MPR in both datasets.
TNFA_SIGNALING_VIA_NFKB, enriched non-MPR in both datasets.

**Overall pattern:**
More nuanced than GSE243013: two distinct clusters with signatures enriched in both MPR and NMPR, reflecting greater TAM transcriptional heterogeneity in this dataset, consistent with the known context-dependency of macrophage programs.

## Preliminary observations : TAMs UCell custom signatures (Bloc 4, GSE243013)

**UCell scores pairwise Wilcoxon Bonferroni (p<0.0167):**
- M2_immunosuppressive: non-MPR > pCR > MPR, all pairs p<0.005
- SPP1_signature: non-MPR > MPR > pCR, all pairs p<0.0001
- IFN_response: non-MPR dominant; MPR ≈ pCR (p=0.783 ns)
- M1_inflammatory: non-MPR > MPR > pCR, MPR/pCR ns Bonferroni

**Key observations:**
- Gradient non-MPR > MPR > pCR confirmed all signatures, MPR as biologically unstable intermediate state
- IFN_response exception: discriminates responders from non-responders, not MPR from pCR
- IFN-stimulated TAMs: hybrid M1/immunosuppressive profile, TME biology, not artifact

**Biological interpretation:**
- Non-MPR: triple immunosuppression (M2 + SPP1 + IFN TAMs)
- pCR more heterogeneous than MPR at patient level, multiple pathways to complete response

## Preliminary observations : UCell Hallmark TAMs (Bloc 4A, GSE243013)

Unbiased scoring of all 50 Hallmark MSigDB signatures on TAMs (MPR vs non-MPR, pCR excluded). 25 biologically relevant discriminant signatures retained after filtering (p_adj < 0.05, Wilcoxon BH-corrected).

**MPR-enriched signature:**
G2M_CHECKPOINT, the only signature enriched in MPR TAMs, though with a modest effect size and the weakest significance of the 25 retained signatures (p_adj = 1.7×10⁻³, versus p_adj < 10⁻⁹ for all other signatures in this panel).

**non-MPR-enriched signatures (24/25):**
Nearly all retained signatures enriched in non-MPR, including metabolic programs (OXIDATIVE_PHOSPHORYLATION, GLYCOLYSIS, FATTY_ACID_METABOLISM), stress programs (UNFOLDED_PROTEIN_RESPONSE, TGF_BETA_SIGNALING, HYPOXIA), oncogenic programs (MYC_TARGETS_V1/V2, KRAS_SIGNALING_UP, MTORC1_SIGNALING, PI3K_AKT_MTOR_SIGNALING), inflammatory/immune programs (TNFA_SIGNALING_VIA_NFKB, INFLAMMATORY_RESPONSE, IL6_JAK_STAT3_SIGNALING, IL2_STAT5_SIGNALING, INTERFERON_ALPHA/GAMMA_RESPONSE, COMPLEMENT, ALLOGRAFT_REJECTION), and other programs (P53_PATHWAY, APOPTOSIS, EPITHELIAL_MESENCHYMAL_TRANSITION, REACTIVE_OXYGEN_SPECIES_PATHWAY, NOTCH_SIGNALING).

**Overall pattern:**
Highly polarized landscape: TAMs in non-MPR show broad and intense activation of stress, metabolic, oncogenic and inflammatory programs. MPR TAMs are transcriptionally quieter, with only G2M_CHECKPOINT showing a modest, borderline enrichment.

## Preliminary observations : TAMs Pseudobulk MPR vs pCR (Bloc 4A, GSE243013)

**Methodological note:** Earlier figures compared non-MPR, MPR and pCR using four custom UCell signatures (M2_immunosuppressive, M1_inflammatory, SPP1_signature, IFN_response) with pairwise Wilcoxon tests, Bonferroni-corrected. This has been superseded: the current analysis uses the same 25 Hallmark signatures already established as discriminant between MPR and non-MPR, restricted to MPR vs pCR only (n=10 vs n=11 patients), tested at the pseudobulk patient level (mean UCell score per patient, avoiding pseudo-replication), Wilcoxon rank-sum test, BH-corrected, consistent with the methodology used throughout this project.

**Result:**
None of the 25 Hallmark signatures discriminant between MPR and non-MPR reached significance between MPR and pCR at the patient level (all p_adj > 0.89). No signature showed even a trend toward separation between the two groups.

**Perspective:**
At the pseudobulk, patient level, the TAM transcriptional programs that distinguish responders from non-responders (Hallmark signatures) do not further separate partial (MPR) from complete (pCR) pathological responders. This suggests that, within the TAM compartment alone, MPR and pCR are not transcriptionally distinct states along the same axis that separates response from non-response. By contrast, exhausted CD8 T cells showed a weak (non-significant) trend toward separating MPR from pCR on a subset of metabolic and stress-related programs (see CD8 pseudobulk analysis above), raising the possibility that residual biological differences between partial and complete pathological response are more likely to reside in the CD8 compartment than in TAM composition, though neither signal reaches formal significance in the current cohorts.

## Preliminary observations : CollecTRI pseudobulk MPR vs pCR, TAMs (Bloc 4A, Script 08, GSE243013)

**Methodological note:** same approach as the CD8 CollecTRI pseudobulk above, testing the same 20 CollecTRI TFs already established as discriminant between non-MPR and MPR (Script 07) between MPR and pCR specifically, within each of the 7 TAM subtypes.

**Result:** none of the 140 tests (20 TFs × 7 subtypes) reached significance (all p_adj > 0.15). This mirrors the earlier UCell Hallmark pseudobulk finding for TAMs (all p_adj > 0.89, Bloc4A Script 06): across two independent analytical layers (functional signature and transcription factor activity), the TAM compartment shows no detectable difference between MPR and pCR.

**Interpretation, held provisional:** the convergence of two independent methods on a null result strengthens the possibility that, within the TAM compartment specifically, MPR and pCR are not transcriptionally distinct states. Combined with the CD8 finding above (8 TFs distinguishing MPR from both non-MPR and pCR in CD8.TEX), this suggests any residual biological signal separating partial from complete pathological response is more likely to reside in the CD8 compartment than in TAM composition or transcriptional activity, consistent with the same suggestion made from the UCell-level pseudobulk analysis.

**Output:** Results/Tables/Bloc4A_08_CollecTRI_pseudobulk_MPR_vs_pCR.csv

## Preliminary observations : CollecTRI TAMs (Bloc 4, both datasets)

### GSE243013

**Violin plots, Top 6 by between-group variance (GSE243013):**
RFXAP, HSF1, RELA, RFXANK, NFKB, STAT1. HSF1, RELA, NFKB and STAT1 show a clear, consistent pattern: lowest in MPR, highest in non-MPR, intermediate in pCR, matching the broad non-MPR-associated program described below. RFXAP and RFXANK show a much weaker version of the same trend at this pooled, all-subtypes level (MPR slightly lower than non-MPR/pCR), which does not on its own reveal the sharp subtype-specific divergence seen in the heatmap (MPR-enriched in Monocyte FCN1+ and Resident M2, non-MPR-enriched in LAMs and Proliferating), a reminder that the pooled violin view can mask opposing subtype-level patterns that partly cancel out.

**Heatmap observations, statistically verified (Wilcoxon MPR vs non-MPR, BH-corrected within each subtype):**

Transcription factor activity inference revealed distinct, statistically verified program-specific patterns across all seven TAM subtypes (Wilcoxon MPR vs non-MPR, BH-corrected within each subtype).

**MPR-enriched (significant), by subtype:**
- Resident M2: RFX5, RFXANK, RFXAP (CIITA not significant, p_adj = 0.895, trending non-MPR, an incomplete MHC II signal, since CIITA is required as coactivator for RFX-bound complexes to transactivate MHC II genes)
- Monocyte FCN1+: RFXANK, RFXAP, RFX5, CIITA (the complete MHC II module, the only subtype where all four align in the same direction)
- Classical-Mono, IFN-stimulated, Stress-response, LAMs, Proliferating: no significant MPR-enriched TF among the Top 20

**non-MPR-enriched (significant), by subtype:**
- LAMs (20/20): STAT1, IRF1, HSF2, REL, HIF1A, SPI1, RELA, NFKB, HSF1, NFKB1, CREB1, DDIT3, CIITA, STAT3, JUN, EGR1, AP1, RFX5, RFXAP, RFXANK, the complete Top 20, including the full MHC II module
- Proliferating (19/20): STAT1, HSF2, HSF1, NFKB, EGR1, SPI1, RELA, NFKB1, DDIT3, HIF1A, IRF1, REL, JUN, CREB1, CIITA, AP1, RFXAP, RFX5, RFXANK, nearly the complete Top 20, including the full MHC II module
- Resident M2 (16/20): HIF1A, NFKB1, RELA, HSF2, NFKB, HSF1, STAT1, EGR1, STAT3, CREB1, JUN, AP1, IRF1, REL, SPI1, DDIT3, everything except the MHC II cluster, which instead trends MPR (see above)
- Monocyte FCN1+ (13/20): HIF1A, HSF1, HSF2, NFKB1, RELA, STAT1, STAT3, NFKB, JUN, AP1, CREB1, SPI1, IRF1, everything except the MHC II cluster, which instead trends MPR (see above)
- IFN-stimulated (16/20): HSF2, HSF1, STAT1, DDIT3, RELA, NFKB1, IRF1, REL, STAT3, HIF1A, CREB1, NFKB, EGR1, SPI1, JUN, AP1 (RFX/CIITA not significant)
- Stress-response (16/20): HIF1A, HSF1, NFKB1, HSF2, RELA, NFKB, AP1, STAT3, IRF1, STAT1, JUN, CREB1, EGR1, SPI1, REL, DDIT3 (RFX/CIITA not significant)
- Classical-Mono (9/20, the least discriminant subtype): HIF1A, HSF2, STAT3, HSF1, NFKB1, IRF1, RELA, NFKB, SPI1

Beyond the MHC II axis, the same broad set of transcription factors (the NF-κB family, NFKB, NFKB1, RELA, REL; STAT1/STAT3; HSF1/HSF2; IRF1; JUN; AP1; CREB1; SPI1; EGR1; DDIT3; HIF1A) is consistently non-MPR-enriched across nearly every subtype, with the degree of completeness scaling with overall subtype discriminance (complete in LAMs/Proliferating, near-complete in Resident M2/Monocyte FCN1+/IFN-stimulated/Stress-response, partial in Classical-Mono).

The clearest response-associated switch in this dataset occurs in Resident M2 and Monocyte FCN1+: both show the broad non-MPR-associated program (HIF1A, NF-κB family, STAT1/3, HSF1/2, etc.) alongside an MHC II signal that runs in the opposite direction from the rest of the subtype's profile, RFX-only (incomplete) toward MPR in Resident M2, and the complete RFX+CIITA module toward MPR in Monocyte FCN1+. LAMs and Proliferating show no such divergence: the MHC II module aligns with the broader non-MPR program in both. This subtype-specific decoupling of MHC II direction from the rest of the transcriptional program will be examined further alongside ligand-receptor communication signals in the CD8-TAM crosstalk section.

### GSE207422

**Two complementary visualizations:**
- **Violin plots (Top 6 by between-group variance, absolute scores)**: global MPR vs non-MPR comparison
- **Heatmaps (scale=row, Top 20 by between-group variance)**, subtype-specific TF patterns per condition

*TF selection performed using the same 20 transcription factors identified by between-group variance for this dataset: RFXAP, RFXANK, RFX5, CIITA, RELA, HSF1, HSF2, NFKB1, NFKBIB, XBP1, REL, CREB1, EGR1, KLF6, IRF1, HIF1A, JUN, JUND, NFKB, IRF5.

**Violin plots, Top 6 by between-group variance:**
RFXAP, RFXANK, HSF1, RFX5, CIITA, RELA, the complete MHC II module (RFX5, RFXANK, RFXAP, CIITA) plus HSF1 and RELA. This composition already hints at what the by-subtype breakdown makes explicit: MHC II activity is not uniformly directional in this dataset, unlike in GSE243013 where it leaned consistently toward non-MPR in most subtypes.

**Heatmap observations, statistically verified (Wilcoxon MPR vs NMPR, BH-corrected within each subtype; Monocyte-derived excluded, 0 MPR cells available):**

**MRC1+ M2-like (19/20 significant, 17 MPR, 2 NMPR):** the most MPR-dominant subtype in either dataset. The complete MHC II module was strongly MPR-enriched (RFX5, RFXAP, RFXANK, CIITA, all p_adj < 10⁻¹⁵⁸), alongside nearly the entire remaining panel, HSF1, HSF2, RELA, NFKB1, NFKBIB, XBP1, REL, CREB1, EGR1, NFKB, HIF1A, JUN, IRF1. Only KLF6 and IRF5 trended NMPR.

**SPP1+ immunosuppressive (19/20 significant, 4 MPR, 15 NMPR):** a clean decoupling. The complete MHC II module (CIITA, RFX5, RFXAP, RFXANK) was significantly MPR-enriched, while the remaining 15 transcription factors (EGR1, NFKBIB, IRF5, HIF1A, JUND, KLF6, XBP1, JUN, NFKB, CREB1, NFKB1, REL, RELA, HSF1, HSF2), were all significantly NMPR-enriched.

**IFN-stimulated (16/20 significant, 2 MPR, 14 NMPR):** an incomplete MHC II signal in the opposite direction from Resident M2 (below): CIITA and IRF1 were significantly MPR-enriched, but none of the three RFX subunits reached significance (all trending MPR but not significant). The remaining 14 TFs (EGR1, HIF1A, JUND, JUN, IRF5, HSF1, HSF2, NFKB, KLF6, XBP1, NFKB1, CREB1, RELA, REL) were NMPR-enriched.

**Resident M2 (10/20 significant, 4 MPR, 6 NMPR):** an incomplete MHC II signal mirroring GSE243013's own Resident M2 subtype. RFXAP, RFXANK and RFX5 were significantly MPR-enriched, but CIITA was not (p_adj = 0.229, trending MPR). HSF2 was also MPR-enriched. EGR1, IRF5, HIF1A, XBP1, KLF6 and NFKBIB were NMPR-enriched.

**M2-SIGLEC8+ (12/20 significant, 0 MPR, 12 NMPR):** NFKBIB, REL, JUND, JUN, NFKB, RELA, IRF5, EGR1, NFKB1, CREB1, HIF1A and HSF1 were NMPR-enriched. The MHC II module (RFX5, RFXANK, RFXAP, CIITA) was not significant.

**Stress-response (1/20 significant, 0 MPR, 1 NMPR):** only JUND reached significance (NMPR). All other TFs, including the MHC II module, showed a non-significant NMPR-leaning trend.

**Lipid-associated and Regulatory (0/20 significant each):** the two least discriminant subtypes in this dataset, with no transcription factor reaching significance in either direction.

**Cross-dataset note:** unlike GSE243013, where non-MPR dominated nearly every subtype (mirroring the global UCell Hallmark pattern, 24/25 signatures NMPR-enriched), GSE207422 shows a more balanced MPR/NMPR split at both the CollecTRI and UCell Hallmark level (9/23 Hallmark signatures MPR-enriched here, versus 1/25 in GSE243013). The RFX/CIITA decoupling pattern observed in Resident M2 (RFX significant, CIITA not) replicates across both datasets in the same subtype, while MRC1+ M2-like's near-total MPR dominance and SPP1+'s clean MHC II/rest-of-program split have no direct counterpart in GSE243013.

### Cross-dataset observations

While transcription factor networks governing CD8 T cell exhaustion and reactivation have been extensively characterized in the context of immunotherapy response, systematic mapping of TF activity across TAM subtypes in relation to pathological response to neoadjuvant immunotherapy remains largely unexplored in NSCLC. TAM subtype annotations mostly differed between cohorts, reflecting the known absence of universal TAM nomenclature in the field (Zhu et al., 2023; Chen et al., 2024), with three notable exceptions: Resident M2, Stress-response and IFN-stimulated carry the same label in both datasets, allowing a direct subtype-level comparison for these three.

**Resident M2 (direct cross-dataset replication):**
In both GSE243013 and GSE207422, Resident M2 shows the RFX complex (RFX5, RFXANK, RFXAP) significantly MPR-enriched, while CIITA does not reach significance in either dataset (p_adj = 0.895 in GSE243013, 0.229 in GSE207422). Since CIITA is required as coactivator for RFX-bound complexes to transactivate MHC II genes, this represents a consistent, cross-dataset-replicated pattern of incomplete MHC II signal specifically in this subtype.

**IFN-stimulated:**
Both datasets show a broad non-MPR-leaning program of similar magnitude (16/20 significant TFs in each). However, the two diverge on two points. First, MHC II: neither RFX nor CIITA reach significance in GSE243013, while in GSE207422 CIITA is significantly MPR-enriched (though without RFX support), an incomplete MHC II signal present only in one cohort. Second, IRF1 itself is non-MPR-enriched in GSE243013 but MPR-enriched in GSE207422, the only transcription factor identified here that reverses direction between datasets in the same subtype.

**Stress-response:**
Both datasets show a non-MPR-leaning direction, though with very different strength: GSE243013 shows a broad, highly significant non-MPR program (16/20 TFs), while GSE207422 shows only one significant TF (JUND, non-MPR) with the rest trending in the same direction without reaching significance. The MHC II module is not significant in either dataset for this subtype.

**Overall MPR/NMPR balance differs substantially between cohorts:**
GSE243013 shows an overwhelmingly non-MPR-dominant TAM compartment, consistent across both CollecTRI (most subtypes with 0 significant MPR-enriched TFs) and Hallmark UCell (24/25 signatures non-MPR-enriched). GSE207422 shows a more balanced picture at both levels (multiple subtypes with substantial MPR-enriched signal, notably MRC1+ M2-like and SPP1+ immunosuppressive; 9/23 Hallmark signatures MPR-enriched). This asymmetry does not appear to be an artifact of a single analysis layer, since it replicates across two independent methods (CollecTRI and UCell) within each dataset.

**Subtype-specific MHC II decoupling is a recurring theme, but its direction is not fixed:**
Beyond Resident M2's consistent RFX-without-CIITA pattern, other subtypes show the complete MHC II module decoupling from the rest of their transcriptional program in either direction: toward MPR in GSE243013's Monocyte FCN1+ and GSE207422's SPP1+ immunosuppressive (complete module MPR-enriched while the rest of the program is non-MPR-enriched), and toward non-MPR in most other discriminant subtypes (complete module aligned with the broader non-MPR program, as in GSE243013's LAMs and Proliferating). No single direction for MHC II generalizes across the TAM compartment as a whole.

*Note: earlier cross-dataset observations referencing AEBP1, ELK4 (in Resident M2), and DMTF1 relied on TF selections superseded by the current between-group-variance method; none of these three TFs are part of the current Top 20 in either dataset and are no longer discussed.*

---

# EPITHELIAL COMPARTMENT (GSE207422 only)

## Preliminary observations : Epithelial compartment CNV (Bloc 4B, GSE207422)

- Epithelial composition: Basal dominant in MPR (~62%), EMT and Tumor_epithelial exclusive to NMPR, transcriptional rather than genomic drivers of non-response
- CNV clonal profiles (SCEVAN) identical between MPR and NMPR: treatment response not determined by chromosomal architecture differences
- Both SCEVAN and CopyKAT confirm higher tumor/aneuploid fraction in NMPR (~47%/~75%) vs MPR (~32%/~38%): directional concordance between tools
- EMT cells: maximal discordance between tools (SCEVAN 0% tumor vs CopyKAT 38% aneuploid); whether EMT resistance is driven by subtle genomic alterations or transcriptional reprogramming alone remains unresolved
- Tumor_epithelial subtype: strong cross-tool concordance (99% tumor/aneuploid both tools)
- Next steps: CytoTRACE, UCell, CollecTRI on confirmed tumor vs normal epithelial cells per condition to characterize transcriptional state differences driving non-response

## Preliminary observations : CytoTRACE2 epithelial differentiation (Bloc 4B, GSE207422, SCEVAN labels)

**Violin plot observations (CytoTRACE2 Score, 0=differentiated, 1=totipotent):**

- Ciliated: very low scores (~0.05-0.10), terminally differentiated, expected for mature bronchial cells
- EMT_NMPR: low scores (~0.10-0.25), majority converging toward differentiation, supports SCEVAN classification (0% tumor); EMT reprogramming appears terminal rather than stem-like
- normal_MPR: pear-shaped distribution with two subpopulations coexist; one highly differentiated and one less differentiated, suggesting epithelial plasticity in the normal compartment under MPR conditions
- normal_NMPR: mass concentrated at low scores, uniformly differentiated, transcriptionally static normal epithelial cells in non-responding patients
- tumor_MPR: compact distribution, high scores (~0.40-0.65), homogeneously stem-like, low differentiation, consistent with a targetable tumor population
- tumor_NMPR: broad bimodal distribution (~0.15-0.80), high intra-tumoral differentiation heterogeneity, consistent with treatment resistance mechanisms

**Biological interpretation:**

In MPR, a dissociation is observed between tumor and normal compartments: tumor cells are homogeneously stem-like while normal epithelial cells display plasticity, suggesting active epithelial remodeling in the treatment-responsive microenvironment. In NMPR, tumor cells show greater heterogeneity while normal cells are uniformly static, a pattern potentially reflecting impaired epithelial dynamics contributing to non-response. EMT_NMPR low differentiation scores support SCEVAN over CopyKAT for this subtype classification, pending CopyKAT-based CytoTRACE2 cross-validation.

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

## Preliminary observations : UCell Hallmark epithelial compartment (Bloc 4B, GSE207422)

**Methodological note:** earlier version used 11 pre-selected Hallmark signatures plus 2 custom signatures (HSF1_targets, Antigen_presentation). Superseded by the unbiased approach: all 50 Hallmark signatures tested, discriminant signatures identified post-hoc by Wilcoxon test (BH-corrected within each comparison), consistent with CD8/TAM methodology.

**tumor_MPR vs tumor_NMPR (33/50 significant):**
*Higher in tumor_MPR*: IL6_JAK_STAT3_SIGNALING (strongest signal, p_adj=3.0×10⁻¹⁷⁵), ALLOGRAFT_REJECTION, PROTEIN_SECRETION, GLYCOLYSIS, OXIDATIVE_PHOSPHORYLATION, WNT_BETA_CATENIN_SIGNALING, NOTCH_SIGNALING, REACTIVE_OXYGEN_SPECIES_PATHWAY, FATTY_ACID_METABOLISM, HYPOXIA, COMPLEMENT, INTERFERON_GAMMA_RESPONSE, INFLAMMATORY_RESPONSE, HEME_METABOLISM, UV_RESPONSE_UP, INTERFERON_ALPHA_RESPONSE.
*Higher in tumor_NMPR*: APICAL_JUNCTION, MITOTIC_SPINDLE, APOPTOSIS, CHOLESTEROL_HOMEOSTASIS, P53_PATHWAY, TNFA_SIGNALING_VIA_NFKB, KRAS_SIGNALING_DN, G2M_CHECKPOINT, UNFOLDED_PROTEIN_RESPONSE, TGF_BETA_SIGNALING, PI3K_AKT_MTOR_SIGNALING, E2F_TARGETS, DNA_REPAIR, KRAS_SIGNALING_UP, MYC_TARGETS_V2, EPITHELIAL_MESENCHYMAL_TRANSITION, MTORC1_SIGNALING.

**normal_MPR vs normal_NMPR (33/50 significant):**
*Higher in normal_NMPR*: KRAS_SIGNALING_DN, PI3K_AKT_MTOR_SIGNALING, HYPOXIA, EPITHELIAL_MESENCHYMAL_TRANSITION, APOPTOSIS, KRAS_SIGNALING_UP, APICAL_JUNCTION, TNFA_SIGNALING_VIA_NFKB, TGF_BETA_SIGNALING, INFLAMMATORY_RESPONSE, MTORC1_SIGNALING, IL2_STAT5_SIGNALING, GLYCOLYSIS, UNFOLDED_PROTEIN_RESPONSE, COMPLEMENT, ANGIOGENESIS, MYC_TARGETS_V1, OXIDATIVE_PHOSPHORYLATION, UV_RESPONSE_UP, INTERFERON_GAMMA_RESPONSE, P53_PATHWAY, PROTEIN_SECRETION, MITOTIC_SPINDLE, INTERFERON_ALPHA_RESPONSE, MYC_TARGETS_V2, HEME_METABOLISM, DNA_REPAIR, REACTIVE_OXYGEN_SPECIES_PATHWAY, WNT_BETA_CATENIN_SIGNALING, UV_RESPONSE_DN (30 of 33).
*Higher in normal_MPR*: FATTY_ACID_METABOLISM, NOTCH_SIGNALING, CHOLESTEROL_HOMEOSTASIS (3 of 33).

**Ciliated_MPR vs Ciliated_NMPR (13/50 significant):**
*Higher in MPR*: INTERFERON_GAMMA_RESPONSE, MYC_TARGETS_V1, APICAL_SURFACE, ALLOGRAFT_REJECTION, PROTEIN_SECRETION, G2M_CHECKPOINT.
*Higher in NMPR*: DNA_REPAIR, INFLAMMATORY_RESPONSE, REACTIVE_OXYGEN_SPECIES_PATHWAY, WNT_BETA_CATENIN_SIGNALING, NOTCH_SIGNALING, OXIDATIVE_PHOSPHORYLATION, TNFA_SIGNALING_VIA_NFKB.

**Key observations:**
- The tumor compartment shows the most extensive divergence between conditions (33/50 signatures), with IL6_JAK_STAT3_SIGNALING as the single most significant signature in this analysis (p_adj=3.0×10⁻¹⁷⁵).
- The normal epithelial compartment shows a comparably extensive divergence (33/50 signatures), but overwhelmingly in one direction: 30 of 33 signatures are higher in normal_NMPR.
- Ciliated cells, despite their original role as a stable reference population for malignancy calling, show 13/50 significant differences between MPR and NMPR, fewer than tumor or normal compartments, but not zero.
- NOTCH_SIGNALING direction differs by population: higher in tumor_MPR, higher in normal_MPR, but higher in Ciliated_NMPR.
- FATTY_ACID_METABOLISM is higher in MPR in both tumor and normal compartments, the only signature showing this same direction in both.

## Preliminary observations : CollecTRI TF activity epithelial compartment (Bloc 4B, Script 11, GSE207422)

**Methodological note:** TF selection (Top 20/Top 6) corrected from per-cell variance to between-group variance (variance of final_group means), consistent with CD8/TAM methodology. This changed the composition of the Top 20 relative to earlier observations, see below.

**Top 20 TFs by between-group variance:**
RFXANK, RFX5, RFXAP, HSF2, FOXJ1, CREB1, TP53, DMTF1, MYC, EPAS1, ZBTB4, MEIS2, ATF1, ELK4, SP1, NFE2L2, ESR1, HIF1A, NKX2-1, DACH1

*Note: CIITA, E2F1, E2F4, HSF1, JUN and NFKB, all featured in earlier observations, are no longer part of the Top 20 under this corrected selection method and are not discussed below.*

**Wilcoxon results, statistically verified (BH-corrected within each comparison):**

**tumor_MPR vs tumor_NMPR (20/20 significant):**
*Higher in tumor_MPR*: ZBTB4, RFXANK, RFX5, RFXAP, EPAS1, NKX2-1, MYC.
*Higher in tumor_NMPR*: HSF2, DMTF1, CREB1, DACH1, ESR1, ELK4, FOXJ1, MEIS2, TP53, ATF1, SP1, HIF1A, NFE2L2.

**normal_MPR vs normal_NMPR (19/20 significant):**
*Higher in normal_MPR*: RFXAP, RFXANK, RFX5, NKX2-1, HSF2, MEIS2.
*Higher in normal_NMPR*: HIF1A, CREB1, DMTF1, EPAS1, TP53, ATF1, MYC, SP1, DACH1, ESR1, ELK4, FOXJ1, NFE2L2.
Not significant: ZBTB4.

**Ciliated_MPR vs Ciliated_NMPR (10/20 significant):**
*Higher in MPR*: RFXANK, RFX5, RFXAP, HSF2, FOXJ1, CREB1, TP53, MYC.
*Higher in NMPR*: DMTF1, EPAS1.
Not significant: ZBTB4, MEIS2, ATF1, ELK4, SP1, NFE2L2, ESR1, HIF1A, NKX2-1, DACH1.

**Key points:**
- The RFX complex (RFXANK, RFX5, RFXAP) is significantly MPR-enriched in all three comparisons, the only TF module showing this same direction across tumor, normal, and Ciliated populations in this dataset.
- CIITA is absent from the Top 20 (as in the CD8 and TAM compartments in both datasets), despite RFX being consistently discriminant.
- ELK4 is significantly NMPR-enriched in both tumor and normal comparisons, not tested in Ciliated (not significant there, p_adj=0.159).
- tumor_MPR_vs_NMPR is the most extensively discriminant comparison in this analysis (20/20 TFs significant).

## Preliminary observations : ELK4/malignancy axis, epithelial compartment (Bloc 4B, GSE207422)

ELK4 transcription factor activity was consistently elevated in malignant (tumor) epithelial cells relative to their normal counterparts, regardless of pathological response (p_adj < 10⁻³⁷ in both MPR and NMPR arms). This is consistent with ELK4's broadly reported pro-oncogenic role across multiple cancer types, and specifically in NSCLC, where ELK4 has been shown to transcriptionally activate MSI2, which in turn drives NSCLC proliferation, migration and invasion through TGF-β/SMAD3 pathway activation (Shi et al., 2025, Kaohsiung J Med Sci).

This mechanistic chain connects to two independent findings already established in this dataset. First, ELK4 itself is NMPR-enriched (not MPR-enriched) in both tumor and normal epithelial cells when compared by pathological response (CollecTRI, see above), the same direction in both compartments. Second, the EPITHELIAL_MESENCHYMAL_TRANSITION Hallmark signature, a functional readout of the TGF-β/SMAD3-driven process ELK4 activates via MSI2, was independently found NMPR-enriched in both tumor and normal epithelial compartments (UCell Hallmark, see above). TGF-β/SMAD signaling is an established driver of EMT specifically in lung cancer, including via SMAD3/4-mediated N-cadherin induction.

Together, three independent layers, TF activity (ELK4), a documented NSCLC-specific mechanistic pathway (ELK4→MSI2→TGF-β/SMAD3), and a functional signature (EMT), converge on the same NMPR-associated program across both malignant and non-malignant epithelial populations. Whether this malignancy-associated, NMPR-enriched ELK4/EMT axis shares any regulatory relationship with the MPR-associated ELK4 signal observed independently in exhausted CD8 T cells (CollecTRI, Bloc 3) remains an open question, given that these represent distinct cell lineages, likely distinct upstream regulatory contexts, and, notably, opposite response-group directions (NMPR in epithelial cells vs MPR in CD8 T cells) despite sharing the same transcription factor.
This divergence, the same transcription factor associating with opposite response-group directions depending on cell lineage, illustrates a general principle of transcription factor biology: a TF's functional consequence depends on the target gene repertoire accessible within a given cell type's chromatin and cofactor landscape, not on a fixed, cell-type-independent function. ELK4 in epithelial cells (via MSI2/TGF-β/SMAD3, driving EMT) and ELK4 in exhausted CD8 T cells (mechanism uncharacterized in this project) likely represent two distinct regulatory programs sharing the same upstream transcription factor, rather than a single pathway operating consistently across the tumor microenvironment.

---

# INTERCELLULAR COMMUNICATION (Bloc 5)

**Note on tool selection:** Both LIANA+ and CellChat v2 were run and are fully documented below for transparency and cross-validation purposes, across both GSE207422 and GSE243013. LIANA+ was selected as the primary framework for the preprint. Its consensus score is computed as a simple product of mean ligand and receptor expression, which could be recomputed directly from the raw expression matrix on a per-patient basis; this allowed a patient-level pseudobulk Wilcoxon test to be added on top of LIANA's own prioritization score, treating each patient as an independent observation and avoiding pseudo-replication at the single-cell level, consistent with the pseudobulk methodology used throughout this project.

CellChat's communication probability, by contrast, is derived from a mass-action model with Hill-function saturation terms, computed internally by the tool itself rather than as a simple, externally reproducible formula. Recomputing this score on a per-patient basis would require re-running CellChat's full algorithm separately for each patient, computationally demanding, and particularly impractical given the smaller response groups (e.g., n=3 MPR patients in GSE207422) typical of neoadjuvant immunotherapy cohorts. CellChat's native permutation test instead pools all cells within each condition without distinguishing between patients — a design vulnerable to pseudo-replication, where a single patient contributing disproportionately more cells could dominate the result — which does not directly support the patient-level comparison this project's methodology requires. CellChat results are retained here as a complementary, qualitative cross-check but are not reported in the preprint.

## GSE207422

### Preliminary observations : LIANA+ communication, cells of interest (Bloc 5, Script 02b)

**Methodological note:** figures restricted to CD8.TEX/CD8.TPEX and the 5 TAM subtypes identified as most discriminant by CollecTRI (MRC1+ M2-like, SPP1+ immunosuppressive, IFN-stimulated, M2-SIGLEC8+, Resident M2). Top 5 interactions per target group prioritized by aggregate_rank ≤ 0.05 (Dimitrov et al., 2022), each tested for MPR vs NMPR difference by patient-level pseudobulk Wilcoxon (BH-corrected). Multi-subunit complexes scored by geometric mean of subunits.

#### CD8 ↔ TAM communication

**CD8 → TAM (42 interactions tested, 38 significant):**

- *IFN-stimulated:* CCL5-CCR1 MPR-enriched in both CD8.TEX and CD8.TPEX; CD52-SIGLEC10 NMPR-enriched in both; GZMB-IGF2R and IFNG-IFNGR1_IFNGR2 NMPR-enriched (CD8.TEX only); TNFSF9-HLA-DPA1 not significant (CD8.TPEX).
- *M2-SIGLEC8+:* CCL5-CCR1 MPR-enriched in both states; CD52-SIGLEC10 NMPR-enriched (CD8.TEX); ENTPD1-ADORA3, GZMB-IGF2R, SPN-SIGLEC1 all MPR-enriched (CD8.TEX); CRTAM-CADM1 MPR-enriched (CD8.TPEX); TNFSF9-HLA-DPA1 not significant (CD8.TPEX).
- *MRC1+ M2-like:* CCL5-CCR1 and IFNG-IFNGR1_IFNGR2 both MPR/NMPR-enriched in a direction opposite to the other four subtypes (CCL5-CCR1 NMPR-enriched, IFNG-IFNGR1_IFNGR2 MPR-enriched, CD8.TEX). CCL5-SDC4 and GZMB-IGF2R NMPR-enriched (CD8.TEX and/or CD8.TPEX); CIRBP-TREM1 NMPR-enriched (CD8.TPEX); TNFSF9-HLA-DPA1 not significant.
- *Resident M2:* CCL5-CCR1 MPR-enriched in both states; CD52-SIGLEC10 and IFNG-IFNGR1_IFNGR2 both NMPR-enriched in both states; GZMB-IGF2R MPR-enriched (CD8.TEX); CCL5-ACKR1, CXCL13-ACKR1 significant with median = 0 in both conditions (direction not interpretable); TNFSF9-HLA-DPA1 not significant.
- *SPP1+ immunosuppressive:* CCL5-CCR1 MPR-enriched in both states; CCL4-CCR1 and IFNG-IFNGR1_IFNGR2 NMPR-enriched (CD8.TPEX); CCL5-CCRL2 NMPR-enriched, SPN-SIGLEC1 MPR-enriched (CD8.TEX); CD200-CD200R1 and TNFSF9-HLA-DPA1 significant with median = 0 in both conditions (direction not interpretable).

**TAM → CD8 (19 interactions tested, 15 significant):**

- *IFN-stimulated:* HLA-C-CD8A NMPR-enriched toward CD8.TEX but MPR-enriched toward CD8.TPEX (opposite direction by target state); SECTM1-CD7 MPR-enriched (CD8.TPEX).
- *M2-SIGLEC8+:* HLA-DQB1-LAG3 NMPR-enriched toward both CD8.TEX and CD8.TPEX; HLA-DQA1-LAG3 NMPR-enriched toward CD8.TPEX (not significant toward CD8.TEX); HLA-DRB3-LAG3 MPR-enriched toward CD8.TEX (not significant toward CD8.TPEX); HLA-A-CD8A and HLA-F-CD8A NMPR-enriched.
- *MRC1+ M2-like:* HLA-DQB1-LAG3 NMPR-enriched (CD8.TPEX; absent in MPR, median = 0).
- *Resident M2:* S100A8-CD69 NMPR-enriched (CD8.TPEX; absent in MPR, median = 0).
- *SPP1+ immunosuppressive:* CCL13-CXCR3 NMPR-enriched toward both CD8.TEX and CD8.TPEX; HLA-DRB3-LAG3 NMPR-enriched toward CD8.TPEX (not significant toward CD8.TEX).

**Recurring pattern:** the HLA-DQ (A1 and B1 subunits)-LAG3 axis is NMPR-enriched consistently across the TAM subtypes where it appears (M2-SIGLEC8+, MRC1+ M2-like). HLA-DQ is an MHC class II molecule and LAG3 is an inhibitory checkpoint receptor whose canonical ligand is MHC class II (Huard et al., 1997, PNAS), one of several molecular families through which TAMs can transmit inhibitory signals to CD8 T cells. This pattern is discussed further alongside the CollecTRI MHC II findings (RFX/CIITA) in the Discussion.

#### Individual HLA-D and LAG3 gene expression, follow-up test (Bloc 5, Script 09/09b)

**Methodological note:** this follow-up analysis was prompted by an apparent discrepancy between CollecTRI's RFX/CIITA regulon activity (MPR-enriched in MRC1+ M2-like and SPP1+ immunosuppressive TAMs) and the NMPR-enriched direction of the HLA-D-LAG3 ligand-receptor communication signal identified by LIANA+ in the same subtypes (Preliminary observations, LIANA+ communication section, above). Individual gene expression (not the ligand-receptor product score) was tested at the patient-level pseudobulk, Wilcoxon rank-sum, BH-corrected, consistent with the methodology used throughout this project.

**Script 09 (TAM side):** HLA-DQB1 and HLA-DRB3 expression tested individually in the three TAM subtypes engaging LAG3 in communication analysis (M2-SIGLEC8+, MRC1+ M2-like, SPP1+ immunosuppressive). No gene reached significance in any subtype (all p_adj ≥ 0.719; M2-SIGLEC8+ underpowered, n_MPR = 2).

**Script 09b (CD8 side):** LAG3 expression tested individually in CD8.TEX and CD8.TPEX. Neither state reached significance (CD8.TEX: p_adj = 0.6; CD8.TPEX: p_adj = 0.6).

**Conclusion:** neither the TAM-side ligand (HLA-DQB1, HLA-DRB3) nor the CD8-side receptor (LAG3) shows a significant, standalone MPR vs NMPR difference at the individual gene level, despite the combined ligand-receptor product score being significant in the original LIANA+ analysis. This indicates the NMPR-enriched direction of the HLA-D-LAG3 communication signal cannot be attributed to a robust, independent expression shift in either gene alone. This finding does not contradict the MPR-enriched RFX/CIITA regulon activity identified by CollecTRI in MRC1+ M2-like and SPP1+ immunosuppressive TAMs, since regulon-level TF activity (aggregated across many target genes) and single-gene expression operate at different analytical levels that are not directly comparable. The apparent discrepancy is most likely attributable to the limited statistical power of single-gene tests in this cohort's small MPR group (n = 3) rather than a genuine biological conflict.

**Output:** Results/Tables/BLOC5/Bloc5_09_HLA_D_individual_genes_wilcox.csv
          Results/Tables/BLOC5/Bloc5_09b_LAG3_CD8_wilcox.csv

#### CD8 ↔ Epithelial communication

**Methodological note:** tumor_MPR/tumor_NMPR and normal_MPR/normal_NMPR were tested using their original condition-encoded labels (patient-level pseudobulk Wilcoxon, since tumor_MPR only exists in MPR patients), then collapsed to "Tumor" and "Normal" for reporting, since both labels of the same axis necessarily yield identical statistics. Ciliated tested directly (condition not encoded in its label).

**CD8 → Epithelial (21 interactions tested, 18 significant):**

- *Ciliated (7/7 significant):* CCL5-SDC4, CCL5-GRM7, GZMA-PARD3, SEMA4D-PLXNB1 all MPR-enriched, in both CD8.TEX and CD8.TPEX where tested. CXCL13-GRM7 significant with median = 0 in both conditions (direction not interpretable).
- *Normal (5/7 significant):* CCL5-SDC4 MPR-enriched; CCL5-DPP4 NMPR-enriched (the only NMPR-leaning interpretable interaction in this axis). CCL5-SDC1 (both states), COL6A3-SDC1, CRTAM-CADM1 significant with median = 0 in both conditions (direction not interpretable).
- *Tumor (3/7 with interpretable direction, all significant):* CCL5-SDC1, HMGB1-SDC1, GZMA-PARD3 all MPR-enriched (CD8.TEX only; absent in NMPR, median = 0). The remaining 4 pairs are significant with median = 0 in both conditions (direction not interpretable) or not significant (XCL1-ADGRV1).

**Epithelial → CD8 (14 interactions tested, 13 significant):**

- *Ciliated → (6/6 significant):* CD59-CD2 MPR-enriched in both states; OMG-TNFRSF1B NMPR-enriched in both states; PCSK1N-GPR171 direction differs by CD8 state (NMPR-enriched toward CD8.TEX, MPR-enriched toward CD8.TPEX).
- *Normal → (1/1 significant):* CADM1-CRTAM significant with median = 0 in both conditions (direction not interpretable).
- *Tumor → (4/7 interpretable, all significant):* HLA-A-CD8A, NECTIN1-CD96 (CD8.TEX and CD8.TPEX), NECTIN4-TIGIT (CD8.TEX) all MPR-enriched, with very low absolute values (0.0001–0.06). CEACAM5-CD8A and GPI-NTRK1 significant with median = 0 in both conditions; NECTIN4-TIGIT toward CD8.TPEX not significant.

**Overall pattern:** across both directions of the CD8 ↔ Epithelial axis, nearly all interactions with an interpretable direction (i.e. excluding pairs with median = 0 in both conditions) are MPR-enriched, with CCL5-DPP4 (CD8 → Normal) as the sole clear exception. Based on these results alone, MPR patients show a greater volume of detectable CD8-epithelial communication than NMPR patients. Whether this reflects a functionally beneficial (e.g. cytotoxic, activating) or simply more active communication landscape has not been characterized here and requires examining the specific nature of each ligand-receptor pair individually, a question taken up in the Discussion alongside the CD8 ↔ TAM axis, which shows the same general pattern.

#### TAM ↔ Epithelial communication

**Methodological note:** same approach as the CD8 ↔ Epithelial axis with tumor_MPR/tumor_NMPR and normal_MPR/normal_NMPR tested using their condition-encoded labels, then collapsed to "Tumor"/"Normal" for reporting. Ciliated tested directly.

**TAM → Epithelial (28 interactions tested, 25 significant):**

- *IFN-stimulated (4/4 significant, all MPR):* CXCL10-SDC4 (Ciliated, Normal), CXCL9-GRM7 (Ciliated), FN1-SDC1 (Tumor).
- *M2-SIGLEC8+ (5/6 significant):* PSAP-CELSR1 (Ciliated) and GRN-EGFR (Tumor) MPR-enriched; HGF-SDC1 and IGF1-IGF1R significant with median = 0 in both conditions (direction not interpretable); HGF-MET not significant.
- *MRC1+ M2-like (7/8 significant):* ADAM17-ERBB4, ALCAM-CHL1, CXCL5-GRM7 (all → Ciliated), FN1-SDC1 and LPL-SDC1 (→ Tumor), LPL-SDC1 (→ Normal) all MPR-enriched, the sole exception being FN1-SDC1 → Normal, which is NMPR-enriched. GRN-EGFR (→ Tumor) not significant.
- *Resident M2 (2/2 significant, both MPR):* MMP9-CD44 (Tumor), PSAP-CELSR1 (Ciliated).
- *SPP1+ immunosuppressive (6/8 significant, nearly all NMPR):* CD209-CEACAM1, CXCL8-SDC1 (Normal and Tumor), F13A1-ITGB1, HGF-SDC1, IGF1-IGF1R all NMPR-enriched. PLAU-ITGA3 significant with median = 0 in both conditions; IGF1-INSR not significant.

**Epithelial → TAM (43 interactions tested, 37 significant):**

- *Ciliated → (24/25 significant, overwhelmingly NMPR):* SLPI-CD4 and SLPI-PLSCR1 recur across nearly every TAM subtype (IFN-stimulated, M2-SIGLEC8+, MRC1+ M2-like, Resident M2, SPP1+), all NMPR-enriched. CD24-SIGLEC10, SAA1-FPR2/SCARB1, PIP-CD4/NPTN, CD59-STAB1, CP-SLC40A1, OMG-RTN4R (→ M2-SIGLEC8+)/OMG-TNFRSF1B all NMPR-enriched. Exceptions (MPR-enriched): FGF14-FGFR1 (→ MRC1+ M2-like and → SPP1+), OMG-RTN4R (→ MRC1+ M2-like), SCGB3A1/A2-MARCO (→ MRC1+ M2-like).
- *Normal → (5/7 significant, mixed):* SCGB3A1-MARCO, SCGB3A2-MARCO (→ MRC1+ M2-like) and SLPI-CD4 (→ SPP1+) MPR-enriched; SLPI-PLSCR1 (→ IFN-stimulated) NMPR-enriched; SLPI-CD4/PLSCR1 (→ Resident M2) and MXRA5-SIGLEC8 (→ M2-SIGLEC8+) significant/not significant with median = 0 in both conditions.
- *Tumor → (7/11 significant, overwhelmingly MPR):* GPC3-CD81 recurs across 4 of the 5 TAM subtypes of interest (IFN-stimulated, M2-SIGLEC8+, MRC1+ M2-like, Resident M2), all MPR-enriched. ANGPTL4-SDC3, MXRA5-SIGLEC8, NTS-SORT1 also MPR-enriched. Remaining pairs (COL4A1/A5-CD93, CX3CL1-CX3CR1, FGF2-SDC3, all → SPP1+/M2-SIGLEC8+) significant with median = 0 in both conditions or not significant.

**Notable contrast:** the direction of the TAM ↔ Epithelial axis depends strongly on which epithelial population is the source. Ciliated → TAM interactions (SLPI axis in particular) are overwhelmingly NMPR-enriched, while Tumor → TAM interactions (GPC3-CD81 in particular) are overwhelmingly MPR-enriched, opposite directions depending on the epithelial source, rather than a uniform MPR- or NMPR-leaning pattern as observed for the CD8 ↔ TAM and CD8 ↔ Epithelial axes. SPP1+ immunosuppressive is also the one TAM subtype whose outgoing signal (TAM → Epithelial) leans consistently NMPR, opposite to the other four subtypes.

### Preliminary observations : LIANA+ bidirectional loops, CD8/TAM ↔ Epithelial (Bloc 5, GSE207422)

**Methodological note:** loops defined as significant, directionally interpretable signaling (excluding median = 0 in both conditions) in both directions, within the same condition and cell pair, applying the same criterion established for CD8-TAM loops.

#### CD8 ↔ Epithelial

Three bidirectional loops identified, all MPR, no NMPR loop identified with any epithelial population:

- **CD8.TEX ↔ Tumor**: CCL5-SDC1, HMGB1-SDC1, GZMA-PARD3 ↔ NECTIN4-TIGIT (specific, TIGIT-exclusive ligand; Reches et al., 2020, J Immunother Cancer; NECTIN4 relevant in NSCLC, contributing to anti-PD-1 resistance via CD155 stabilization, Springer 2025), NECTIN1-CD96 (lower-affinity than CD96's primary ligand CD155; functional direction context-dependent), HLA-A-CD8A.
- **CD8.TEX ↔ Ciliated**: CCL5-SDC4, CCL5-GRM7, GZMA-PARD3, SEMA4D-PLXNB1 (confirmed T cell semaphorin, CD100) ↔ CD59-CD2 (established second CD2 ligand alongside CD58, T cell adhesion/activation; Deckert et al., 1992, Eur J Immunol).
- **CD8.TPEX ↔ Ciliated**: CCL5-SDC4, SEMA4D-PLXNB1 ↔ CD59-CD2, PCSK1N-GPR171 (GPR171/BigLEN is a genuine T cell checkpoint pathway suppressing TCR signaling; Ma et al., 2021, Nat Commun; direction here is MPR-enriched, opposite to its NMPR-enriched engagement of CD8.TEX noted in the CD8-TAM section, direction differs by CD8 state for this same pair).
- **No loop with Normal epithelial cells**: CD8 signaled outward (CCL5-SDC4 MPR; CCL5-DPP4 NMPR) without any interpretable reciprocal response.

#### TAM ↔ Epithelial

Four MPR loops with Tumor (IFN-stimulated, M2-SIGLEC8+, MRC1+ M2-like, Resident M2, all subtypes tested), one MPR loop with Ciliated (MRC1+ M2-like), one NMPR loop with Ciliated (SPP1+ immunosuppressive), one partial MPR loop with Normal (MRC1+ M2-like).

**Tumor loops (all MPR):**
- IFN-stimulated ↔ Tumor: FN1-SDC1 ↔ GPC3-CD81
- **M2-SIGLEC8+ ↔ Tumor** (most extensively characterized; every tested pair essentially MPR-exclusive, median_NMPR ≈ 0 throughout): GRN-EGFR (progranulin, pro-tumorigenic EGFR-directed signaling, under clinical investigation in NSCLC, anti-progranulin antibody trial NCT05627960; macrophage-derived EGFR ligands documented to drive NSCLC therapy resistance), HGF-SDC1 (HGF/c-MET regulates EMT, cancer stem cell properties and immune evasion in NSCLC TME; connects to the CytoTRACE2-identified MPR-specific tumor/normal stemness dissociation) ↔ MXRA5-SIGLEC8 (not supported by current receptor biology, SIGLEC8's characterized ligands are sialylated glycans, not MXRA5; excluded from interpretation), GPC3-CD81 (specific CD81-GPC3 axis, characterized primarily in HCC; GPC3 shows moderate, clinically relevant expression specifically in squamous NSCLC), ANGPTL4-SDC3 (direct HIF-1α target gene, strongly hypoxia-correlated in NSCLC; converges with EPAS1 and HYPOXIA signature already established as MPR-enriched in this dataset's tumor cells, though hypoxia is otherwise classically linked to NSCLC treatment resistance, Cancers 2021; no mechanistic explanation for the MPR association is proposed).
- MRC1+ M2-like ↔ Tumor: FN1-SDC1, LPL-SDC1 ↔ GPC3-CD81, NTS-SORT1
- Resident M2 ↔ Tumor: MMP9-CD44 ↔ GPC3-CD81

**Ciliated loops:**
- MRC1+ M2-like ↔ Ciliated (MPR): CXCL5-GRM7, ADAM17-ERBB4, ALCAM-CHL1 ↔ SCGB3A1-MARCO, SCGB3A2-MARCO, OMG-RTN4R, FGF14-FGFR1
- SPP1+ immunosuppressive ↔ Ciliated (NMPR, the only NMPR loop identified in this axis): IGF1-IGF1R ↔ SLPI-CD4, CD59-STAB1, CP-SLC40A1

**Normal (partial MPR loop only):**
- MRC1+ M2-like ↔ Normal: LPL-SDC1 ↔ SCGB3A1-MARCO, SCGB3A2-MARCO (FN1-SDC1 from TAM side is NMPR, conflicting direction, not part of the coherent loop)
- SPP1+ immunosuppressive → Normal (NMPR, unidirectional): HGF-SDC1, CXCL8-SDC1, F13A1-ITGB1, no coherent reciprocal signal (Normal's only response, SLPI-CD4, is MPR, opposite direction)

**GPC3-CD81 recurs across 4 of 5 TAM subtypes with Tumor**, the single most shared signal in this axis; CD81 is GPC3's principal binding partner (Hippo/HHEX pathway regulation), characterized mainly in hepatocellular carcinoma, with uncharacterized relevance to tumor-TAM communication specifically.

**Perspectives note:** These epithelial-compartment findings derive from a single dataset (GSE207422; epithelial cells not available in GSE243013) and should be interpreted as hypothesis-generating rather than confirmed. The near-total absence of NMPR loops with Tumor/Normal cells (against a background of frequent unidirectional or condition-discordant signaling) contrasts with the CD8-TAM axis, where NMPR loops were common; whether this reflects a genuine biological asymmetry or a power/detection limitation specific to this cohort's small MPR group (n=3) cannot be resolved here. Because TAMs, rather than CD8 T cells, showed the more extensive and more frequently reciprocated communication with tumor and normal epithelial cells in this dataset (see also the general absence of CD8-Normal loops), TAMs are considered a central axis for future investigation, ideally combined with spatial transcriptomics to resolve physical cell-cell proximity underlying these co-expression-based signals.

### Preliminary observations : CellChat, cells of interest (Bloc 5, Script 04b)

**Methodological note:** CellChat's native permutation-based significance test (nboot=100) is applied independently within each condition (MPR, NMPR); it provides no native statistical test to directly compare an interaction's strength between the two groups. Only interactions classified as MPR_only or NMPR_only (present and significant in one condition, absent from the significant set in the other) are reported below; interactions present in both conditions ("Both") are not further discriminated by CellChat's native output and are excluded here.

#### CD8 ↔ TAM communication

Full listing, including "Both" (416/693 tested), in Bloc5_04b_CellChat_GSE207422_CD8_TAM_crossed.csv.

**MPR_only (74 interactions):**
- CD8.TEX → M2-SIGLEC8+: CCL3-CCR5, CCL4-CCR5, CCL5-CCR5, HLA-DMA-CD4, PTPRC-CD22
- CD8.TEX → MRC1+ M2-like: CD96-PVR, FASLG-FAS, HLA-DMA-CD4, PGE2-PTGES3-PTGER4, PTPRC-CD22
- CD8.TEX → SPP1+ immunosuppressive: CCL4-CCR5, CCL5-CCR5, CD99-CD99, FASLG-FAS, HLA-DMA-CD4
- CD8.TEX → IFN-stimulated: HLA-DMA-CD4
- CD8.TEX → Resident M2: FASLG-FAS
- CD8.TPEX → M2-SIGLEC8+: CCL3-CCR5, CCL4-CCR5, CCL5-CCR5, COL6A3-SDC4, COL6A3-CD44, HLA-DRB3-CD4, PTPRC-CD22
- CD8.TPEX → MRC1+ M2-like: CD96-PVR, COL6A3-ITGAV_ITGB8, COL6A3-CD44, COL6A3-SDC4, FASLG-FAS, HLA-DRB3-CD4, PTPRC-CD22
- CD8.TPEX → SPP1+ immunosuppressive: CCL3-CCR5, CCL4-CCR5, CD99-PILRA, COL6A3-CD44, FASLG-FAS, HLA-DRB3-CD4, HLA-F-LILRB1
- CD8.TPEX → IFN-stimulated: COL6A3-CD44, COL6A3-SDC4, HLA-DRB3-CD4
- CD8.TPEX → Resident M2: COL6A3-SDC4, COL6A3-CD44, FASLG-FAS, HLA-DRB3-CD4, TNFSF9-TNFRSF9
- MRC1+ M2-like → CD8.TEX: CD274-PDCD1, MDK-ITGA4_ITGB1, NECTIN2-CD226, PDCD1LG2-PDCD1, PVR-CD226, PVR-TIGIT, TNF-TNFRSF1B
- MRC1+ M2-like → CD8.TPEX: CD274-PDCD1, ICAM1-ITGAL, ICAM1-ITGAL_ITGB2, LGALS9-HAVCR2, PDCD1LG2-PDCD1, PVR-CD226, PVR-TIGIT, TNF-TNFRSF1B
- Resident M2 → CD8.TEX: CCL5-CCR5, CD274-PDCD1, FN1-ITGA4_ITGB7, FN1-ITGA4_ITGB1, NECTIN2-CD226, PDCD1LG2-PDCD1
- Resident M2 → CD8.TPEX: CD274-PDCD1, PDCD1LG2-PDCD1
- IFN-stimulated → CD8.TEX: NECTIN2-CD226
- IFN-stimulated → CD8.TPEX: LGALS9-HAVCR2
- M2-SIGLEC8+ → CD8.TEX: NECTIN2-CD226
- M2-SIGLEC8+ → CD8.TPEX: LGALS9-HAVCR2
- SPP1+ immunosuppressive → CD8.TEX: CD99-CD99

**NMPR_only (203 interactions):**
- CD8.TEX/TPEX → all 5 TAM subtypes: ANXA1-FPR1, ANXA1-FPR3, CD6-ALCAM (recurring across nearly every subtype)
- CD8.TEX/TPEX → IFN-stimulated: additionally ANXA1-FPR2, ANXA1-FPR2_LXA4, HLA-DRB5-CD4, MT-RNR2-FPR2, SEMA7A-PLXNC1 (TPEX)
- CD8.TEX/TPEX → M2-SIGLEC8+: additionally CD99-CD99L2, CRTAM-CADM1, HLA-DRB5-CD4, PGE2-PTGES3-PTGER2/PTGER4, SEMA7A-PLXNC1 (TPEX)
- CD8.TEX/TPEX → MRC1+ M2-like: additionally HLA-DRB5-CD4, HLA-DPB1-CD4 (TPEX), CD99-CD99 (TPEX), SEMA7A-PLXNC1 (TPEX)
- CD8.TEX/TPEX → Resident M2: additionally CCL3-CCR5, CRTAM-CADM1, HLA-DQA1-CD4, HLA-DPB1-CD4 (TPEX), HLA-DRB5-CD4, PGE2-PTGES3-PTGER2, PTPRC-MRC1, SEMA4D-CD72, SEMA7A-PLXNC1 (TPEX), CD99-PILRA (TPEX)
- CD8.TEX/TPEX → SPP1+ immunosuppressive: additionally CD99-CD99L2, ENTPD1-TMIGD3, HLA-DRB5-CD4, PGE2-PTGES3-PTGER4, SEMA7A-PLXNC1 (TPEX)
- All 5 TAM subtypes → CD8.TEX/TPEX: CCL3-CCR1, CD55-ADGRE5, HLA-A/B/C-CD8B (recurring across nearly every subtype)
- IFN-stimulated → CD8: additionally CCL5-CCR5, CCL7-CCR1, CCL8-CCR1, CXCL9/10/11-CXCR3, ITGA4_ITGB1/ITGB7-VCAM1, SPP1-CD44, SPP1-ITGA4_ITGB1, CDH1-ITGAE_ITGB7 (TPEX), Cholesterol-LIPA-RORA (TPEX), FN1-ITGA4_ITGB7 (TPEX), PGE2-PTGES2-PTGER4 (TPEX)
- M2-SIGLEC8+ → CD8: additionally CCL23-CCR1, CCL5-CCR5, CCL8-CCR1, CXCL9/10-CXCR3, MDK-ITGA4_ITGB1, PGE2-PTGES2/PTGES3-PTGER4, CDH1-ITGAE_ITGB7 (TPEX), Cholesterol-LIPA-RORA (TPEX), F11R-ITGAL_ITGB2 (TPEX), LGALS9-CD44 (TPEX)
- MRC1+ M2-like → CD8: additionally CCL23-CCR1, ITGA4_ITGB1-VCAM1, SPP1-CD44, SPP1-ITGA4_ITGB1, CD99-CD99 (TPEX), CDH1-ITGAE_ITGB7 (TPEX), Cholesterol-DHCR24/LIPA-RORA (TPEX), F11R-ITGAL_ITGB2 (TPEX), FN1-ITGA4_ITGB7 (TPEX), HLA-E-CD8B (TPEX), ITGAV_ITGB1-ADGRE5 (TPEX), PGE2-PTGES2-PTGER4 (TPEX), THBS1-CD47 (TPEX)
- Resident M2 → CD8: additionally CCL8-CCR1, CXCL9/10-CXCR3, F11R-ITGAL_ITGB2, LGALS9-CD44, PGE2-PTGES3-PTGER4, RETN-CAP1, SIGLEC1-SPN, CDH1-ITGAE_ITGB7 (TPEX), Cholesterol-LIPA-RORA (TPEX), ITGAV_ITGB1-ADGRE5 (TPEX), PGE2-PTGES2-PTGER4 (TPEX)
- SPP1+ immunosuppressive → CD8: additionally CCL13-CCR1, CCL3L3-CCR1, CCL5-CCR5, CCL8-CCR1, CD80-CTLA4, COL1A1/COL1A2-CD44/ITGA1_ITGB1, FN1-ITGA4_ITGB1/ITGB7, ITGA4_ITGB1-VCAM1, NECTIN2-TIGIT, SPP1-CD44, SPP1-ITGA4_ITGB1, TNF-TNFRSF1B, TNFSF9-TNFRSF9, CDH1-ITGAE_ITGB7 (TPEX), COL1A1-SDC4 (TPEX), COL1A2-SDC4 (TPEX), Cholesterol-LIPA-RORA (TPEX), NECTIN2-CD226 (TPEX)

#### CD8 ↔ Epithelial communication

Full listing, including "Both" (159/314 tested), in Bloc5_04b_CellChat_GSE207422_CD8_Epithelial_crossed.csv.

**MPR_only (83 interactions):**
- CD8.TEX → Ciliated: CD99-CD99/CD99L2, HLA-DMA/DPA1/DPB1/DRA/DRB1-CD4
- CD8.TEX → normal: CRTAM-CADM1
- CD8.TEX → tumor: CD99-CD99, SEMA4D-PLXNB1
- CD8.TPEX → Ciliated: CD99-CD99L2, COL6A3-ITGA2_ITGB1/ITGA3_ITGB1/ITGAV_ITGB8/SDC4, HLA-DPA1/DRA/DRB1/DRB3-CD4
- CD8.TPEX → normal: COL6A3-ITGA3_ITGB1/CD44/SDC1/SDC4, CRTAM-CADM1
- CD8.TPEX → tumor: COL6A3-ITGA2_ITGB1/ITGAV_ITGB8/CD44/SDC1/SDC4, SEMA4D-PLXNB1
- Ciliated → CD8.TEX: CD99-CD99, COL9A2-ITGA1_ITGB1/CD44, ICAM1-ITGAL/SPN/ITGAL_ITGB2, LAMA5-ITGA1_ITGB1/CD44, NECTIN2-CD226
- Ciliated → CD8.TPEX: COL9A2-ITGA1_ITGB1/CD44/SDC4, LAMA5-ITGA1_ITGB1/CD44
- normal → CD8.TEX: NECTIN2-CD226
- normal → CD8.TPEX: APP-CD74, ICAM1-ITGAL_ITGB2
- tumor → CD8.TEX: CD274-PDCD1, CD99-CD99, COL4A1/COL9A3-ITGA1_ITGB1/CD44, CXCL12-CXCR4, LAMB1/LAMB3-ITGA1_ITGB1/CD44, LGALS9-HAVCR2/PTPRC, LPAR1/LPAR2/LPAR3-ADGRE5, NECTIN2-CD226, NECTIN3-TIGIT
- tumor → CD8.TPEX: ANGPTL4-SDC4, CD274-PDCD1, COL4A1-ITGA1_ITGB1/CD44/SDC4, COL9A3-ITGA1_ITGB1/CD44/SDC4, CXCL12-CXCR4, LAMB1/LAMB3-ITGA1_ITGB1/CD44, LGALS9-PTPRC, LPAR1/LPAR2/LPAR3-ADGRE5, NECTIN3-TIGIT

**NMPR_only (72 interactions):**
- CD8.TEX → Ciliated: CD6-ALCAM, HLA-DRB5-CD4
- CD8.TEX → normal: CD6-ALCAM, GZMA-PARD3, IFNG-IFNGR1_IFNGR2
- CD8.TEX → tumor: CD6-ALCAM, TNFSF10-TNFRSF10A
- CD8.TPEX → Ciliated: CD6-ALCAM, Cholesterol-LIPA-RORA, SIRPG-CD47
- CD8.TPEX → normal: CD6-ALCAM, GZMA-PARD3, IFNG-IFNGR1_IFNGR2, SIRPG-CD47
- CD8.TPEX → tumor: CD6-ALCAM, Cholesterol-LIPA-RORA, SIRPG-CD47, TNFSF10-TNFRSF10A
- Ciliated → CD8.TPEX: CDH1-ITGAE_ITGB7, Cholesterol-DHCR24/LIPA-RORA, HLA-A/B/C/E/F-CD8B, MDK-SDC4, PGE2-PTGES3-PTGER4
- normal → CD8.TEX: ITGAV_ITGB1-ADGRE5, LAMA5/LAMB3/LAMC1/LAMC2-ITGA1_ITGB1/CD44, PGE2-PTGES-PTGER4
- normal → CD8.TPEX: CD55-ADGRE5, CDH1-ITGAE_ITGB7, Cholesterol-DHCR24/LIPA-RORA, HLA-A/B/C/E-CD8B, HLA-F-CD8A/CD8B, LAMA5/LAMB3/LAMC1/LAMC2-ITGA1_ITGB1/CD44, PGE2-PTGES-PTGER4
- tumor → CD8.TEX: CD55-ADGRE5, KLK6-F2R, SPP1-ITGA4_ITGB1
- tumor → CD8.TPEX: CDH1-ITGAE_ITGB7, Cholesterol-DHCR24/DHCR7/LIPA-RORA, HLA-A/B/C/E/F-CD8B, PGE2-PTGES2/PTGES3-PTGER4, SPP1-ITGA4_ITGB1

#### TAM ↔ Epithelial communication

**Methodological note:** only MPR_only/NMPR_only interactions reported (792/1491 tested); "Both" (699/1491) excluded. Full listing in Bloc5_04b_CellChat_GSE207422_TAM_Epithelial_crossed.csv.

**MPR_only (467 interactions):**
- Ciliated → all 5 TAM subtypes: HLA-DMA/DMB/DOA/DPA1/DPB1/DQA1/DQB1/DRB5-CD4 (recurring MHC-II module, all subtypes); COL9A2-CD44/SDC4, ICAM1 variants, MPZL1-MPZL1, JAG2-NOTCH2, LAMA5-CD44 also recurring
- IFN-stimulated/M2-SIGLEC8+/MRC1+ M2-like/Resident M2/SPP1+ → Ciliated: APP-CD74, CD99-CD99/CD99L2, HBEGF-EGFR/EGFR_ERBB2, MPZL1-MPZL1, PPIA-BSG, SEMA4A/SEMA4D-PLXNB1/PLXNB2 (recurring across subtypes)
- Same 5 TAM subtypes → normal: APP-CD74/TNFRSF21, HGF-MET, LGALS9-CD44, IL1B-IL1R2 (recurring)
- Same 5 TAM subtypes → tumor: APP-SORL1, PECAM1-CD38, SEMA4A/SEMA4D-PLXNB1, LTB4-LTA4H-LTB4R, NAMPT-ITGA5_ITGB1 (recurring); OSM-LIFR_IL6ST/OSMR_IL6ST (M2-SIGLEC8+, SPP1+)
- normal → all 5 TAM subtypes: HLA-DMA/DMB/DOA/DPA1/DPB1/DQA1/DQB1/DRA/DRB1/DRB5-CD4 (complete MHC-II module, recurring across all subtypes); SCGB3A2-MARCO (IFN-stimulated, M2-SIGLEC8+, MRC1+ M2-like)
- tumor → all 5 TAM subtypes: ANGPTL4-ITGA5_ITGB1/SDC3/SDC4, COL4A1/COL9A3-CD44/SDC4, CXCL12-CXCR4, IGFBP3-TMEM219, LAMB1/LAMB3-CD44, LPAR1/LPAR2/LPAR3-ADGRE5, NAMPT-ITGA5_ITGB1, NECTIN3-NECTIN2, PLAU-PLAUR, SEMA4C-PLXNB2 (all recurring across the 5 subtypes)

**NMPR_only (325 interactions):**
- Ciliated → TAM subtypes: SEMA3C variants (M2-SIGLEC8+), MDK-SDC2/ITGA4_ITGB1 (recurring), GAS6/PROS1/TUB-MERTK (MRC1+ M2-like)
- IFN-stimulated/M2-SIGLEC8+/MRC1+ M2-like/Resident M2/SPP1+ → Ciliated: Cholesterol-LIPA/DHCR24-RORA (recurring); SPP1-ITGAV_ITGB1/ITGB5 (IFN-stimulated, MRC1+ M2-like)
- Same 5 subtypes → normal/tumor: TNFSF12-TNFRSF12A (recurring across nearly all); SPP1-CD44/ITGAV variants (recurring); PGE2-PTGES2-PTGER4 (recurring); HBEGF-EGFR_ERBB2, CDH1-ITGA2_ITGB1 (recurring)
- SPP1+ immunosuppressive → Ciliated/normal/tumor: extensive collagen signaling (COL1A1/COL1A2-ITGA2_ITGB1/ITGA3_ITGB1/CD44/SDC1/SDC4), FN1-ITGA3_ITGB1/ITGAV_ITGB8/SDC1/SDC4, TNF-TNFRSF1A, CLDN1-CLDN1, TGM2-ADGRG1, CD209-CEACAM1 (recurring, this subtype the most extensively NMPR-skewed toward epithelium)
- normal → all 5 TAM subtypes: ITGAV_ITGB1-ADGRE2/ADGRE5, LAMA5/LAMB3/LAMC1/LAMC2-CD44 (complete laminin module, recurring across all subtypes), IGFBP3-TMEM219, PGE2-PTGES-PTGER2/PTGER4, PLAU-PLAUR, MDK-SDC2 (recurring)
- tumor → all 5 TAM subtypes: GDF15-TGFBR2, MDK-SDC2, SPP1-CD44/ITGA4_ITGB1/ITGA5_ITGB1/ITGAV_ITGB1/ITGB5, SEMA4D-PLXNB2/CD72, ADO-NT5E_SLC29A1-ADORA3, VEGFA/VEGFB-FLT1 (IFN-stimulated), CDH1-CDH1 (recurring)

**Notable pattern:** the complete MHC class II module (HLA-DMA/DMB/DOA/DPA1/DPB1/DQA1/DQB1/DRA/DRB1/DRB5-CD4) recurs as MPR-enriched in both directions (Ciliated↔TAM and normal↔TAM), across nearly all 5 TAM subtypes, and it is the most extensive and consistent MHC-II signal identified across all communication axes analyzed in this dataset.

## GSE243013

> **Status: DRAFT, NOT YET FINALIZED.** The following LIANA+ and CellChat observations for GSE243013 were produced before the current methodology (patient-level Wilcoxon test on LIANA+, CD8 restricted to TEX/TPEX, biological verification of ligand-receptor pairs before literature interpretation, MPR_only/NMPR_only framework for CellChat) was established on GSE207422. They are retained below for reference only and will be reworked once the GSE207422 CD8-TAM communication asymmetry Discussion section is complete. **Do not cite these observations in the preprint as currently written.**

*(GSE243013 LIANA+ and CellChat draft content retained here, unchanged, pending rework — see prior README version for full text)*

---

# BULK VALIDATION

## Preliminary observations : Bulk validation (Bloc 6)

**Status: not performed.**

A systematic search for publicly available bulk RNA-seq datasets of neoadjuvant anti-PD-1 treated NSCLC with MPR/non-MPR pathological response classification was conducted. The following options were evaluated and excluded:

1. **GSE126044** (Kim et al., Genome Medicine 2020, nivolumab, n=16): RECIST-based classification (CR/PR vs SD/PD), not pathological response (MPR/non-MPR). RECIST underestimates pathological response; in Checkmate 159, 45% achieved MPR while only 10% showed partial response on CT/RECIST. Classification mismatch prevents valid comparison with scRNA-seq findings.

2. **IMvigor210**: urothelial carcinoma cohort (atezolizumab/anti-PD-L1), wrong cancer type and wrong drug class, not applicable.

3. **Hu et al. Genome Medicine 2023 bulk cohort** (n=21): pre-treatment biopsies collected before neoadjuvant ICI, baseline context, not post-treatment surgical resection. Response prediction context differs from post-neoadjuvant MPR classification.

4. **TCGA-LUAD**: LUAD cohort available, but without the specific neoadjuvant anti-PD-1 treatment context required. TCGA cohorts largely predate widespread neoadjuvant immunotherapy protocols, so treatment-naive surgical resection data cannot substitute for a post-neoadjuvant MPR/non-MPR classification.

**Conclusion:** No publicly available bulk RNA-seq dataset combining NSCLC, neoadjuvant anti-PD-1, and MPR/non-MPR pathological response classification was identified. Bulk validation remains a priority for future work, contingent on availability of harmonized datasets or institutional cohorts.



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

## Author

**Myriam Yasmina Soumahoro**
MSc Biology, University of Geneva | Mandela Washington Fellow 2025, Arizona State University
[GitHub](https://github.com/yasmina-bioinfo)


 
 