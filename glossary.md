# R Functions Glossary: TumorImmune_Crosstalk_LUAD

## base R
- `length()`: number of elements in a vector
- `dim()` : dimensions of a matrix or dataframe (rows x cols)
- `head()` :show first 6 rows of an object
- `readLines()` : read raw lines of a file as plain text, without parsing. 
  Used to inspect file structure before loading (headers, format, separators).
  Example: readLines(file, n = 5) reads only the first 5 lines.
- `dir.create()` : create a directory
- `file.path()` : build a file path from components
- `table()` : contingency table of counts
- `setdiff()` : return elements in the first vector that are not in the second.
  Used here to identify expected columns that are missing from the metadata.
- `colSums(is.na())` :count the number of NA values per column. Used here to detect missing clinical data per variable.

## data.table
- `fread()` : fast read of large files (csv, gz) into a data.table object
- `uniqueN()` : count number of unique values in a column.
Return unique rows of a data.table or unique values of a vector. When used with `by = "sampleID"`, keeps only the first row per unique sampleID.
  Used here to deduplicate from cell-level to patient-level (336685 rows → 63 rows).
- `..cols_patient` : data.table syntax to select columns by a vector of names.
The `..` prefix tells data.table to look for the variable in the parent environment, not inside the table itself.
  Example: `dt[, ..my_cols]` selects only the columns listed in `my_cols`.

## BPCells
- `import_matrix_market()` : read MTX matrix from disk without loading into RAM
- `write_matrix_dir()` : save BPCells matrix to disk in BPCells format
- `t()` : transpose a matrix (swap rows and columns). 
  Used here because the MTX file is stored as cells x genes, 
  but Seurat expects genes x cells.
  - `unlink()` : delete a file or directory from disk. 
  Use recursive = TRUE to delete a folder and all its contents.
  Used here to remove an incomplete BPCells output before rewriting.
- `convert_matrix_type()` : convert BPCells matrix to a specific numeric type.
  type = "uint32_t" stores values as unsigned 32-bit integers (whole numbers only),
  required for efficient compression. Raw counts are always integers.

## Concepts
- **Embedding** : a lower-dimensional representation of the data. PCA embedding: 50 PCs. Harmony embedding: batch-corrected PCs. UMAP embedding: 2D visualization. Each is a different representation of the same cells in a progressively reduced space.

## ProjecTILs
- `get.reference.maps()` : download and load all available ProjecTILs reference atlases.
  Returns a nested list organized by species (human/mouse) and cell type (CD8, CD4, DC, MoMac).
  NOTE: downloads ~2GB total on first use; run once, references cached locally.
  Example: ref <- get.reference.maps()$human$CD8

- `Run.ProjecTILs()` : project query T cells onto a reference TIL atlas.
  ref: reference Seurat object (e.g. human CD8 atlas).
  filter.cells = TRUE: removes cells too distant from reference (low quality projection).
  Output: adds functional.cluster column to Seurat object with T cell state labels.
  NOTE: uses RNA assay by default — log-normalization applied internally.

- `list.reference.maps()` : list all available ProjecTILs references with metadata
  (species, cell type, filename, figshare ID).

## sctype
- `gene_sets_prepare()` : load and prepare cell type marker database from sctype.
  Takes a database file path and tissue type as input.
  Example: gene_sets_prepare(db_path, "Lung")

- `sctype_score()` : compute cell type scores for each cell based on marker expression.
  scRNAseqData: scaled expression matrix (genes x cells).
  gs: positive markers list. gs2: negative markers list.
  NOTE: requires standard matrix — not compatible with BPCells directly.

## Seurat
- `subset()` : filter a Seurat object by cells or by metadata conditions
Used here to extract T cell clusters from the full TME object.
  Example: subset(seu, subset = seurat_clusters %in% c(1,2,3))
- `merge()` : combine two or more Seurat objects into one
- `DimPlot()` : visualize cells in a 2D embedding (UMAP, PCA, etc.).
  group.by: color cells by a metadata column.
  label: show cluster labels on the plot.
  raster = FALSE: in DimPlot(), disables automatic point rasterization. Default behavior rasterizes points when >100,000 cells (image becomes blurry). raster = FALSE keeps each point as a vector object; sharper at any zoom level.
  Recommended for publication-quality figures.
  split.by: split the plot by a metadata variable.
- `@reductions <- list()` : reset all dimensional reductions (PCA, UMAP, Harmony) stored in the Seurat object. Required after subsetting to force recalculation on the new cell population from scratch.

- `@graphs <- list()` : reset all neighbor graphs stored in the Seurat object. Required after subsetting for the same reason as @reductions.

## Seurat / Harmony
- `RunHarmony()` : batch correction in PCA space across patients.
  group.by.vars = "sampleID": corrects inter-patient technical variation.
  Output: a new "harmony" embedding that replaces "pca" for downstream steps.
  
- `FindNeighbors()` : builds a nearest-neighbor graph between cells.
  reduction = "harmony": uses Harmony-corrected embedding (not raw PCA).
  dims = 1:40: number of PCs used, chosen based on ElbowPlot inspection.

- `FindClusters()` : groups cells into clusters based on neighbor graph.
  resolution: controls cluster granularity. Lower = fewer broader clusters.
  Start conservative (0.5) for global TME annotation, refine later.

- `RunUMAP()` : projects cells into 2D space for visualization.
  reduction = "harmony": uses Harmony-corrected embedding.
  Does not change the biology — only the visualization.

- `JoinLayers()` : merges multiple layers into one in a Seurat v5 object.
Required after merge() which creates one layer per sample.Must be called before FindAllMarkers() or any layer-dependent operation.
NOTE: Seurat v5 replaced 'slot' with 'layer' — use layer = "data" instead of slot = "data" in GetAssayData().

