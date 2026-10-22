# Comparing chemistries — `catcheR_compare_chemistries()`

**Template:** `inst/rmd/catcheR_chemistry_comparison_EB.Rmd`
**Status:** knitted on the three runs below with this release's code - 0.2 min, since it re-streams nothing.

---

## What it is for

Three entry points — `catcheR_tenxcatch()`, `catcheR_pipcatch()`,
`catcheR_parsecatch()` — return objects of the same shape, because everything
downstream of read parsing is shared code. That makes a like-for-like
comparison a table lookup rather than a new analysis: this function reads the
`.rds` each run wrote and produces one table and one three-panel figure.

Nothing is re-streamed. Run the three pipelines once each, then the comparison
knits in seconds.

## Usage

```r
C <- catcheR_compare_chemistries(
  x = list("10x"    = "analysis_EB_10X/Output/tenxcatch_AB019a.rds",
           "PIPseq" = "analysis_EB_PIP/Output/pipcatch_AB019b.rds",
           "Parse"  = "analysis_EB_PARSE/Output/parsecatch_AB019c.rds"),
  called = c("10x"    = "analysis_EB_10X/scratch/tenx_called_barcodes.txt",
             "PIPseq" = "analysis_EB_PIP/scratch/AB019b_gex.scRNA.filtered.barcodes.tsv.gz",
             "Parse"  = "analysis_EB_PARSE/scratch/cell_metadata.csv.gz"),
  output_dir = "comparison_Output")
```

The list names are the figure labels. Two chemistries, or four runs of the same
one, work the same way — nothing assumes there are exactly three.

## The two rules the function enforces

**1. The denominator is the platform's own called-cell list.** Counting against
the barcodes catcheR merely observed includes the ambient floor and inflates
coverage — on the AB019c Parse library, 50,675 observed barcodes against 4,040
called cells. If `called` is not supplied for a run, the per-cell columns come
back `NA` rather than being computed against the wrong denominator.

**2. Depth is not comparable across chemistries without a label.** A chemistry
whose molecular index is long enough to deduplicate reports molecules; one
whose index is too short reports reads. The function detects this from the
object (a per-molecule read column is present or it is not), reports it as
`umi_dedup`, and names the read-counted runs in the depth panel's subtitle.
`events` means molecules where `umi_dedup` is `TRUE` and reads where it is
`FALSE` — do not quote the two as one number.

## Columns

| column | meaning |
|---|---|
| `read_pairs` | pairs streamed |
| `cellbc_reads`, `pct_cellbc` | reads with a parseable (and, where a whitelist was given, corrected) cell barcode |
| `cassette_reads` | reads where the cassette was found |
| `assigned_reads`, `pct_assigned` | reads ending with a barcode present in the reference |
| `events` | deduplicated molecules, or reads where `umi_dedup` is `FALSE` |
| `called_cells` | cells the platform's pipeline called |
| `called_with_cassette` | of those, with at least one cassette event |
| `called_one_clone`, `pct_called_one_clone` | of those, resolving to exactly one clone — **the number for a methods section** |
| `median_depth_called` | median cassette molecules (or reads) per called cell |
| `constructs`, `clones` | distinct shRNAs and distinct barcode–UCI combinations detected |

## How to read the figure

**Panel A, where the reads go.** Each stage is a strict subset of the one
above; the axis is logarithmic because libraries can differ by three orders of
magnitude. A library that collapses at *cell barcode parsed* has an amplicon
problem upstream of anything catcheR does — a primer-dimer product carries no
capture oligo, so its cell read is adapter read-through.

**Panel B, clone calls per called cell.** The comparable denominator.

**Panel C, depth per called cell.** Read the subtitle before comparing the
violins; see rule 2 above.

## Result on the three AB019 libraries

Same shRNA library (160 constructs), same parameters (`percentage` 15,
`ratio` 2):

| | 10x AB019a | PIPseq AB019b | Parse AB019c |
|---|---|---|---|
| read pairs | 27.79 M | 18.16 M | 21.45 M |
| with a cell barcode | 88.0% | 0.07% | 84.1% |
| with a barcode in reference | 66.4% | 0.02% | 74.2% |
| cells called by the platform | 8,566 | 5,509 | 4,040 |
| **resolving to one clone** | **6,178 (72.1%)** | **1,210 (22.0%)** | **2,857 (70.7%)** |
| median depth per called cell | 221 UMI | 1 **read** | 90 UMI |
| constructs / clones | 102 / 4,068 | 92 / 178 | 100 / 4,474 |

10x and Parse are equivalent on the number that matters. PIPseq is three orders
of magnitude behind on read yield, and panel A localises the loss precisely: at
the cell-barcode stage, which is the primer-dimer problem described in
`pipeline-pipseq.md` — not a chemistry limit and not a thresholding choice.
