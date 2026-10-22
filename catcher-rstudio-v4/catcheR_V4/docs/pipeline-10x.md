# The 10x Genomics pipeline — `catcheR_tenxcatch()`

**Template:** `inst/rmd/catcheR_tenxcatch_analysis_EB.Rmd`
**Status:** knitted end to end on AB019a (27.8 M read pairs) with this release's code - 15.3 min, 60 files written.

---

## What the reads look like

A 10x barcode-enrichment library carries the cell identity on **read 1** and
the CATCHER cassette on **read 2** (verified against `.tenx_adapter()`, which
sets `cell_read = 1`):

```
read 1   [cell barcode 16] [UMI 12] ...
read 2   [rc(reference) 23] [rc(UCI) 6] [rc(barcode) 8] [vector constant 20] ...
```

Read 1 is at fixed offsets — exactly the slice the original shell pipeline cut
with `cut -c 1-16` and `cut -c 17-28` — which is why 10x needs no barcode
*parsing*, only barcode *correction*.

On read 2 the cassette fields are located relative to the **vector constant**,
not to the reference. This matters: a read with a sequencing error inside the
23 nt reference anchor is still recovered, with the right barcode, because the
constant downstream of the barcode is intact.

## What the function does

One streaming pass over the pair, in chunks, in base R:

1. parse read 1 into cell barcode + UMI, and correct the barcode against the
   whitelist to `max.edit.distance` (1 by default; ambiguous neighbours are
   dropped, not assigned);
2. classify read 2 — cassette, empty vector, or neither;
3. slice UCI and shRNA barcode upstream of the constant and look the barcode up
   in the merged reference;
4. deduplicate on (cell, UMI, barcode, UCI) — 10x's 12 nt UMI is a real
   molecular index, so events are molecules;
5. tabulate per (cell, UCI), apply the calling rule, write every table and plot.

Steps 3–5 are the shared `.cat_*` layer: identical code runs for PIPseq and
Parse. Only steps 1–2 are chemistry-specific.

## Arguments that change the answer

| argument | default | what it does |
|---|---|---|
| `cell.whitelist` | `NULL` | barcode list for error correction. Without it no correction happens and a single sequencing error becomes a barcode belonging to no real cell. |
| `called.cells` | `NULL` | the denominator for the final panels. Without it those panels are skipped rather than computed against the wrong denominator. |
| `percentage` | 15 | minimum share of a cell's UMIs a UCI must hold. |
| `ratio` | 2 | dominant-to-second ratio required once a cell is down to one surviving UCI. This is **not** catcheR's historical default of 5 — see "Choosing the calling thresholds" in the README. |
| `UMI.count` | 1 | minimum UMIs on a single (cell, UCI) pair, applied as a strict `>`. |
| `min.cell.UMI` | 20 | minimum cassette UMIs for a cell to be analysed at all. |
| `max.edit.distance` | 1 | correction radius per barcode. |

The same file can serve as both `cell.whitelist` and `called.cells`; on 10x it
usually should, since Cell Ranger's filtered barcode list is both.

Accepted cell-list formats: `barcodes.tsv.gz`, a plain text file of IDs, or
`filtered_feature_bc_matrix.h5` (needs the `hdf5r` package). A `-1` suffix is
stripped on read.

## What it writes

Into `output_dir`:

* `tenxcatch_<sample>.rds` — the run object, input to the filter step and to
  `catcheR_compare_chemistries()`;
* the 13 catcheR plots under their original names, as `.jpg` and `.pdf`;
* `tenxcatch_*.png` — the panels specific to this generation of the pipeline:
  detection funnel, barcode rank, saturation, UCI per cell, call accuracy,
  threshold sweep, clone matrix, QC distributions and thresholds, and — when
  `called.cells` is supplied — called-cell composition, the UMI-floor panel
  restricted to called cells, and yield versus expected wrong calls;
* the tables: `complete_table_fin`, per-cell, clone summary, read-class summary.

## Downstream

```r
flt <- catcheR_tenx_filtercatch(res, UMI.count = 1, percentage = 15, ratio = 2,
                                min.cell.UMI = 20, output_dir = OUT)
noc <- catcheR_tenx_nocatch(res, filtered = flt)
obj <- catcheR_tenx_seurat(res, filtered = flt, empty = noc,
                           mtx = "matrix.mtx.gz", features = "features.tsv.gz",
                           cells = "barcodes.tsv.gz")
```

These three are thin wrappers: the implementation is shared across chemistries
and lives in `R/catcheR_parse_functions.R`. They exist so a 10x script reads as
10x throughout.

## Relation to `catcheR_10Xcatch()`

The original entry point is still in the package and unchanged. It drives the
bundled shell scripts through `system()` and writes intermediates into the
*input* folder. `catcheR_tenxcatch()` is its native-R replacement and is the
recommended path: it writes only into `output_dir`, needs no shell tools, and
shares its downstream half with the other chemistries, so a cross-chemistry
comparison is a single function call. The legacy template is kept as
`inst/rmd/catcheR_10Xcatch_legacy_shell.Rmd` and is marked unexecuted.

## Where the outputs land

`output_dir` is nested one level: everything goes to `<output_dir>/tenxcatch/`,
so another chemistry writing into the same `Output/` cannot overwrite the plots
that carry catcheR's own file names. Each re-filter adds a tagged subfolder of
its own (`filtercatch_p15_u1_r2_f20/`). Pass `subdir = FALSE` for a flat layout.

## Known limits

* Correction is to edit distance 1 by default and only against the list you
  supply. With no whitelist, barcodes are taken as read.
* The barcode-rank panel is less informative on 10x than elsewhere, because the
  whitelist *is* the called-cell list, so the knee and the cell count agree by
  construction.
