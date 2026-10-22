# The Fluent PIPseq pipeline — `catcheR_pipcatch()`

**Template:** `inst/rmd/catcheR_pipcatch_analysis_EB.Rmd`
**Step-by-step derivation:** `inst/rmd/catcheR_pipcatch_report.Rmd`
**Status:** knitted end to end on AB019b (18.2 M read pairs) with this release's code - 5.4 min, 54 files written.

---

## What the reads look like

PIPseq puts the cassette on **read 2** and the cell identity on **read 1** —
the mirror of 10x:

```
read 1   [0-3 nt phasing] bc1(8) ATG bc2(6) GAG bc3(6) TCGAG bc4(8) [index 3] ACG (dT)n
read 2   [rc(reference) 23] [rc(UCI) 6] [rc(barcode) 8] [vector constant 20] ...
```

Two features of read 1 drive the whole port:

* **The barcode is built from four ligation tiers** separated by constant
  linkers (`ATG`, `GAG`, `TCGAG`), preceded by 0–3 phasing bases. The fields
  are therefore **not at fixed offsets**, so the 10x slice (`cut -c 1-16`)
  cannot be reused — the barcode has to be parsed with a regular expression.
  The emitted cell ID is the 28 nt concatenation `bc1 bc2 bc3 bc4`.
* **The molecular index is 3 nt followed by a constant `ACG`**, not a 12 nt
  UMI. There is effectively nothing to deduplicate on, so every per-cell count
  in a PIPseq run is a **read count**, not a molecule count. This propagates:
  `events` in the comparison table are reads, the depth panel says so, and the
  automatic threshold estimators have far less to work with.

On the AB019b library the phasing ladder observed was `''` / `A` / `CA` / `TCA`
rather than the documented `None` / `T` / `GT` / `TGA`; tiers 2 and 3 resolved
to 96-sequence whitelists as expected. If your kit differs, check which
chemistry `pipseeker --chemistry` was run with (v3/v4/v5 parse differently).

## Joining to the expression matrix

PIPseeker's filtered `barcodes.tsv.gz` contains the **same 28 nt concatenated
four-tier barcode** that `catcheR_pipcatch()` emits — verified on AB019b, where
all 5,509 entries were 28 nt. The join is a plain string match; no translation
through `metrics/barcodes/barcode_whitelist.txt` is needed. The function warns
if it ever sees a different width.

## Arguments that change the answer

| argument | default | what it does |
|---|---|---|
| `called.cells` | `NULL` | PIPseeker's cell list. Supplying it unlocks the panels whose denominator is "cells the platform called". |
| `percentage` | 15 | minimum share of a cell's reads a UCI must hold. |
| `ratio` | 2 | dominant-to-second ratio once a cell is down to one surviving UCI. |
| `UMI.count` | 0 recommended here | minimum reads on a (cell, UCI) pair, strict `>`. At one read per cell, any positive value discards most of the data. |
| `min.cell.UMI` | 2 for this chemistry | minimum cassette reads per cell. The 20 used on Parse and 10x excludes every PIPseq cell, so the default falls back automatically. |

## What it writes

Into `output_dir`: `pipcatch_<sample>.rds`, the 13 catcheR plots as `.jpg` and
`.pdf`, the `pipcatch_*.png` analysis panels, and the tables. Identical
inventory to the other chemistries — only the prefix changes.

## Downstream

```r
flt <- catcheR_pip_filtercatch(res, UMI.count = 0, percentage = 15, ratio = 2,
                               output_dir = OUT)
noc <- catcheR_pip_nocatch(res, filtered = flt)
obj <- catcheR_pip_seurat(res, filtered = flt, empty = noc,
                          mtx = "matrix.mtx.gz", features = "features.tsv.gz",
                          cells = "barcodes.tsv.gz")
```

## What AB019b showed, and why it matters for your bench work

The library is dominated by a 20 nt primer-dimer product. Because that product
carries no capture oligo, its read 1 is adapter read-through rather than a
barcode read — so the loss shows up at the *cell barcode* stage, before the
cassette is even looked for: only 12,126 of 18,161,354 pairs (0.07%) carried a
parseable four-tier barcode, and 0.02% ended with a reference barcode.

Of PIPseeker's 5,509 called cells, 1,210 (22.0%) resolved to exactly one clone,
against ~71% for 10x and Parse on the same shRNA library. The median called
cell had **1 read**. This is a library-preparation problem, not a chemistry
limit or a thresholding problem: a size selection on the enrichment PCR to drop
the 20 nt product is the fix. No choice of `percentage`, `ratio` or
`UMI.count` recovers information that was never sequenced.

## Where the outputs land

`output_dir` is nested one level: everything goes to `<output_dir>/pipcatch/`,
so another chemistry writing into the same `Output/` cannot overwrite the plots
that carry catcheR's own file names. Each re-filter adds a tagged subfolder of
its own (`filtercatch_p15_u1_r2_f20/`). Pass `subdir = FALSE` for a flat layout.

## Known limits

* **No UMI deduplication**, as above. Every count is reads.
* **Barcode tiers are validated structurally, not against PIPseeker's
  whitelists.** A single error in `bc1` or `bc4` produces a barcode that passes
  the structural check and belongs to no real cell. Correcting each tier
  against the published whitelist would be the improvement; on this library the
  depth problem dominates it.
* `catcheR_pip_nocatch()` cannot build a control population at this depth —
  AB019b yielded 23 empty-reference reads in 17 cells. The function reports
  this rather than pretending otherwise.
