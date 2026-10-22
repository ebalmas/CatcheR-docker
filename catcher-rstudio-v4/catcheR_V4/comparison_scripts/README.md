# Library comparison — standalone scripts

Deliberately **outside** the catcheR package. The package produces the general
per-run results; these two files take one of its run objects and answer two
experiment-specific questions, so they can be edited per experiment without
touching the package or reinstalling anything.

| file | what it is |
|---|---|
| `catcheR_library_comparison.R` | the functions: `compare_bc_libraries()`, `compare_sublibraries()`, plus local copies of the small helpers they need |
| `catcheR_library_comparison.Rmd` | a parameterised report that sources the script and runs both on one run object |

Needs `ggplot2` and `patchwork`; `hdf5r` only if you pass a Cell Ranger `.h5`
as the cell list. Nothing calls catcheR internals, so the scripts keep working
if the package changes underneath them.

## Use

Run the chemistry's own template first — these read the `.rds` it wrote and
re-stream nothing.

```r
# as a script
source("catcheR_library_comparison.R")
res <- readRDS("analysis_EB_PIP/Output/pipcatch_AB019b.rds")
L <- compare_bc_libraries(res,
       called = "analysis_EB_PIP/scratch/AB019b_gex.scRNA.filtered.barcodes.tsv.gz",
       output_dir = "comparison_out", file = "pipseq_AB019b_bc_libraries")

# or as a report
rmarkdown::render("catcheR_library_comparison.Rmd", params = list(
  run_object   = "analysis_EB_PARSE/Output/parsecatch_AB019c.rds",
  called_cells = "analysis_EB_PARSE/scratch/cell_metadata.csv.gz",
  tag          = "parse_AB019c"))
```

## The two questions

* **`compare_bc_libraries()`** — the pooled BC libraries (`set1`, `set2` …):
  recovery against each library's own size, dropout, evenness (CV, top-decile
  share, Gini, ranked curve), cells per construct, and cells carrying calls
  from more than one library as candidate doublets. Also writes a
  `_per_construct.csv`, one row per construct, zero-filled so dropouts are rows
  of zero rather than missing rows.
* **`compare_sublibraries()`** — the sequencing aliquots (`s1`, `s2` …): reads,
  share reaching the cassette, cells, share resolving to one clone, median
  depth, constructs, clones, and clone sharing. A run with one FASTQ pair gets
  a message, not a one-bar chart.

Two things worth knowing before reading the output:

1. **Pass `called`.** Without a called-cell list the counts include the ambient
   floor, which is not evidence about either library.
2. **Control constructs share a name across libraries.** `SCR1`, `SCR2`,
   `LUCIFERASE`, `B2M` and `ARHGDIA` appear in both CSVs — 10 of the 160 names
   on this reference. Keying on the shRNA name alone credits one detection to
   both libraries and inflates every per-library count. The script keys on
   `(shRNA, set)` and uses the library the cell's dominant call came from.

## Verified on

| run | BC libraries | sublibraries |
|---|---|---|
| AB019b (PIPseq) | set1 51/100, set2 39/60; 4.0 vs 14.6 cells per construct | one sequencing library - section declines |
| AB019c (Parse) | set1 48/100, set2 40/60; 10.9 vs 49.1 cells per construct | s1 2,076 cells / 72.9% one clone / median 188; s2 1,961 / 68.5% / 53 |
