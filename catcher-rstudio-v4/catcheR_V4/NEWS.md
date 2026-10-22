# catcheR 0.4.0 (V4)

Major release. Three chemistries now run through one native-R implementation,
and the runs can be compared to each other.

## New

* **`catcheR_tenxcatch()`** — native R entry point for 10x, replacing the
  shell-driven `catcheR_10Xcatch()` (which is unchanged and still supported).
  Streams the FASTQ in R, corrects the cell barcode against a whitelist,
  deduplicates on the 12 nt UMI, and writes only into `output_dir`. With it,
  `catcheR_tenx_filtercatch()`, `catcheR_tenx_nocatch()`,
  `catcheR_tenx_seurat()` — thin chemistry-named wrappers over the shared
  implementations.
* **`catcheR_compare_chemistries()`** — one table and one three-panel figure
  from run objects you already have. Re-streams nothing. Enforces two rules:
  the denominator is the platform's own called-cell list, and per-cell depth is
  labelled as molecules or reads depending on whether the chemistry's index is
  deduplicable.
* **`called.cells`** on every entry point — unlocks the panels whose
  denominator is "cells the platform called" rather than "barcodes we
  observed". Without it those panels are skipped rather than computed against
  the wrong denominator.
* **Templates**: `catcheR_tenxcatch_analysis_EB.Rmd` and
  `catcheR_chemistry_comparison_EB.Rmd`.
* **Per-pipeline explanations** in `docs/`: `pipeline-10x.md`,
  `pipeline-pipseq.md`, `pipeline-parse.md`, `pipeline-comparison.md`. Plain
  markdown, rendered by GitHub, no build step.

## Fixed

* `catcheR_pipcatch_analysis_EB.Rmd` did not pass `called.cells`, so its
  called-cell panels were never produced.
* Output filenames now follow the class of the object throughout: a 10x or
  Parse object run through a shared implementation no longer produces files
  called `pip_*`. (The example `Output/` folders knitted for this release were
  produced just before this fix, so their `*_nocatch_*.rds` still carries the
  old prefix; a fresh run writes the new one.)
* The called-composition headline read the `pooled` row, which single-library
  chemistries no longer emit, and printed `NA`. It is now a cell-weighted mean
  over whichever rows exist.
* Three fixes in the threshold-sweep panel, all found by reading the rendered
  figure on a read-counted library: the fallback floor picked the *highest*
  value that retained cells - the most restrictive - drawing the panel on 4
  cells instead of 2,018; it now takes the value nearest the one requested. The
  title asserted that "lowering the percentage mostly converts 'no call' into
  'ambiguous'", which is what the Parse library shows and is false where the
  curves are flat; it now states the measured change (on AB019b, percentage 20
  to 10 moves 5 cells; on AB019c, 158). And every panel helper defaulted
  `prefix` to `"parsecatch"`, so a direct call attributed the figure to the
  wrong chemistry; the default is now neutral.
* **The threshold definitions were ambiguous and are now exact.** `UMI.count`
  was documented as "minimum UMIs" while applied as a strict `>`, so
  `UMI.count = 1` actually keeps pairs of **two or more** UMIs; `percentage` was
  "minimum share" and is likewise strict, so a UCI at exactly 15% does not pass
  `percentage = 15`; `min.cell.UMI` is the one inclusive comparison (`>=`). All
  three are now stated with their comparison and a worked consequence, in a
  shared `@section` on the entry points and the re-filter, in
  `docs/data-flow.md`, and in the templates. The documentation also states what
  each step defaults the floor to: both the entry points and the re-filter use
  20, so it keeps applying unless overridden (`0` disables it).
* **The re-filter shows how to choose the floor when you have not chosen it.**
  Leaving `min.cell.UMI` at its default makes `catcheR_*_filtercatch()` say which
  value it used and write `catcheR_call_accuracy_<sample>.csv` plus
  `<chemistry>_yield_vs_error.png` next to its results, so 20 can be confirmed or
  replaced on the dataset in front of you rather than taken on trust. Passing the
  argument suppresses both. The calls are identical either way. The entry points
  now also keep the resolved called-cell list on the object, so that panel uses
  the called-cell denominator instead of falling back to every observed barcode,
  ambient floor included (on AB019c, 2,857 rather than 9,313 at no floor).
* **An input manifest.** Each entry point writes
  `<sample>_<chemistry>_input_files.csv` and stores it on the object as
  `$input_files`: role, resolved path, size and modification time for every
  input it opened. It joins the handoff set copied to `scratch/`. An input
  record, not a checksum - it cannot prove sample identity, but it does catch a
  reference from the wrong library or a FASTQ replaced since the run.
* **`docs/data-flow.md`**, and a flow section at the top of each chemistry
  template: the four ordered steps with what each one reads and produces, why
  the order is forced (`nocatch` subtracts the called set, so it cannot precede
  `filtercatch`), which steps are cheap to repeat, where files land, and which
  three files are handoffs rather than deliverables.
* **`to.scratch` on the re-filter step.** A `.csv` in `Output/` is a
  deliverable; a `.csv` the next step reads is an input, and in catcheR's layout
  those belong in `scratch/`. Two of them - `cells_with_at_least_1_UCI.csv` and
  `new_cell_names.csv` - carry neither the chemistry nor the sample in their
  name, which is fine inside one chemistry's folder and wrong the moment two
  pipelines share a scratch. Give `to.scratch` a path and those files, plus
  `complete_table_fin_<sample>.csv`, are COPIED there as
  `<name>_<chemistry>_<sample>.csv`; the originals stay in `output_dir`. The
  three analysis templates pass it by default (`to_scratch: true`). Verified by
  running all three chemistries into one scratch: nine files, no collisions.
* **One subfolder per chemistry, and one per threshold set.** Every entry point
  now writes to `<output_dir>/<chemistry>/` (`subdir = TRUE`), so running the
  10x and the Parse pipeline into the same `Output/` no longer has the second
  replace the first: the 13 plots that carry catcheR's own file names have no
  chemistry in the filename, by design. Re-filtering is a working step and gets
  its own tagged subfolder, `filtercatch_p15_u1_r2_f20/`, so successive attempts
  at different thresholds coexist instead of overwriting one another; the path
  is returned as `$output_dir`. Everything downstream reads `x$output_dir`, so
  the nesting propagates on its own, and the analysis templates take `OUT` from
  the object instead of rebuilding the path. Pass `subdir = FALSE` for the old
  flat layout.
* Panel helpers create their output directory, so regenerating figures from a
  saved `.rds` no longer fails when the original run's folder is gone.
* **`catcheR_threshold_sweep()` now calls the filter instead of restating its
  rule.** It re-implemented the criteria with `>=` on the percentage - the
  filter uses a strict `>`, as catcheR does - and applied no ratio test, so it
  advertised a yield the filter would not deliver: on AB019c it reported 3,676
  cells where the filter validates 3,500. The table now carries both
  `one_clone` (one UCI above threshold) and `valid` (ValidCell), and `valid`
  matches `catcheR_*_filtercatch()` exactly at the same parameters.
* **`catcheR_call_accuracy()` now scores the subsample on the same rule**:
  strict `>` on both criteria plus the dominant-to-second ratio with the
  100%-share escape clause. Scoring a looser rule made the error curve
  optimistic about the threshold it is used to choose. Error among calls on
  AB019c moves from 1.92% to 1.14% at a 20-UMI floor and from 1.08% to 0.25%
  at 30 - the ratio test now rejects ambiguous draws instead of letting them be
  wrong calls. Figures quoted before this release were computed on the looser
  rule.
* A **guard against silent overwriting**: the 13 plots written under catcheR's
  own file names carry no sample in the filename, by design, so a second sample
  written into the same `output_dir` used to replace the first without a word.
  It now warns, naming the sample already present.
* **`gem.well`** on the Seurat join, for `cellranger aggr` matrices. `count`
  and each per-sample output of `multi` carry a single suffix (`-1`), safe to
  drop; `aggr` pools several wells as `-1`/`-2`/..., where the same 16 nt
  barcode appears once per well and dropping the suffix would merge cells from
  different wells. Since catcheR's cell IDs come from one enrichment library,
  pass the suffix of the well this run corresponds to and the matrix is
  restricted to it. Without it, an aggregated matrix is refused with that
  instruction rather than silently merged. To annotate a whole aggregated
  object, run once per well and merge afterwards, attaching the metadata before
  the merge.
* The Seurat join now accepts a Cell Ranger `filtered_feature_bc_matrix.h5`
  directly as `mtx` (`features` and `cells` are then ignored, and a multimodal
  `.h5` uses its first assay). Needs `hdf5r`, via `Seurat::Read10X_h5()`.
* `subset.valid = TRUE` now reports what it drops, by status, instead of
  silently returning a smaller object. The matrix and catcheR do not cover the
  same cells: the matrix holds every cell the chemistry called, catcheR only
  those whose cassette survived the filter, and `catcheR_status` marks each one.
* `catcheR_tenx_seurat()` was routed through `catcheR_parse_seurat()`, whose
  reader expects split-pipe's `all_genes.csv.gz` / `cell_metadata.csv.gz`, while
  its own documentation told you to pass Cell Ranger's `features.tsv.gz` /
  `barcodes.tsv.gz`. It now routes through `catcheR_pip_seurat()`, which reads a
  MatrixMarket triple with `Seurat::ReadMtx()`.
* The Seurat join did not strip Cell Ranger's GEM-well suffix (`-1`) from the
  barcodes, so on 10x **every** cell failed to match and the function reported
  0 calls on a valid matrix. The suffix is now stripped from both the cell list
  and the count-matrix column names; a collision after stripping is an error,
  not a silent merge. Without the second half of that fix the per-cell summary
  printed correct numbers while the object itself carried a metadata frame of
  all `NA`.
* The barcode-width check in the Seurat join was hardcoded to PIPseq's 28 nt and
  warned on every 10x run. It now expects 16 nt for `tenxcatch`, 28 nt for
  `pipcatch`, and skips the check for Parse well-triplet IDs.
* A diagnostic that cannot run on a given library — the accuracy test needs
  cells deep enough to serve as ground truth, the saturation curve needs a
  deduplicable index — was skipped silently. It now says which panel is missing
  and why. On PIPseq this accounts for three of the panels the other two
  chemistries produce.
* Several plotting defaults were tuned to one chemistry's depth and broke on
  another: axis breaks per discrete level, the zoom window, the per-cell UMI
  floor, the threshold sweep on a read-counted library, a duplicate `pooled`
  facet, and the barcode-rank label running off the panel.

## Verification status of the templates

Knitted end to end with the 0.4.0 code, on the HEDGe libraries:

| template | library | read pairs | time |
|---|---|---|---|
| `catcheR_tenxcatch_analysis_EB.Rmd` | AB019a (10x) | 27.8 M | 15.3 min |
| `catcheR_pipcatch_analysis_EB.Rmd` | AB019b (PIPseq) | 18.2 M | 5.4 min |
| `catcheR_parsecatch_analysis_EB.Rmd` | AB019c (Parse) | 21.5 M | 8.8 min |
| `catcheR_chemistry_comparison_EB.Rmd` | the three above | — | 0.2 min |

`catcheR_10Xcatch_legacy_shell.Rmd` has **not** been executed and is kept for
reference only.

## Result on those three libraries

Same shRNA library, same parameters (`percentage` 15, `ratio` 2), cells counted
against each platform's own called-cell list:

| | 10x AB019a | PIPseq AB019b | Parse AB019c |
|---|---|---|---|
| reads with a cell barcode | 88.0% | 0.07% | 84.1% |
| reads with a barcode in reference | 66.4% | 0.02% | 74.2% |
| called cells resolving to one clone | 6,178 / 8,566 (72.1%) | 1,210 / 5,509 (22.0%) | 2,857 / 4,040 (70.7%) |
| median depth per called cell | 221 UMI | 1 **read** | 90 UMI |

---

# catcheR 0.3.0

PIPseq (`pipcatch`) and Parse (`parsecatch`) ports, chemistry-independent
threshold diagnostics (`call_accuracy`, `threshold_sweep`), and the shared
`.cat_*` layer underneath both.
