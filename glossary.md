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

## data.table
- `fread()` : fast read of large files (csv, gz) into a data.table object
- `uniqueN()` : count number of unique values in a column

## BPCells
- `import_matrix_market()` : read MTX matrix from disk without loading into RAM
- `write_matrix_dir()` : save BPCells matrix to disk in BPCells format
- `t()` : transpose a matrix (swap rows and columns). 
  Used here because the MTX file is stored as cells x genes, 
  but Seurat expects genes x cells.
  - `unlink()` — delete a file or directory from disk. 
  Use recursive = TRUE to delete a folder and all its contents.
  Used here to remove an incomplete BPCells output before rewriting.
- `convert_matrix_type()` — convert BPCells matrix to a specific numeric type.
  type = "uint32_t" stores values as unsigned 32-bit integers (whole numbers only),
  required for efficient compression. Raw counts are always integers.