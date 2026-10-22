# The data flow: what reads what, and where it lands

Two questions this answers: *which step must run before which*, and *which file
is a deliverable versus an input to the next step*. Both matter more than usual
here because a barcode-enrichment analysis and an expression analysis meet
halfway, and because three chemistries can share one experiment folder.

---

## The four steps

| # | step | reads | produces | needed by |
|---|---|---|---|---|
| 1 | `catcheR_<chem>catch()` | FASTQ pair(s), barcode-to-gene CSVs, chemistry extras (`barcode_data.csv` for Parse, a called-cell list for 10x) | the run object (`.rds`), `complete_table_fin`, 13 catcheR plots, the QC panels | everything below |
| 2 | `catcheR_<pfx>_filtercatch()` | the object from 1 | `ValidCell`, `cells_with_at_least_1_UCI.csv`, `new_cell_names.csv`, the threshold figures | steps 3, 4 |
| 3 | `catcheR_<pfx>_nocatch()` | the object from 1 **and** the filter from 2 | the empty-vector (unperturbed) cells | step 4, optional |
| 4 | `catcheR_<pfx>_seurat()` | the filter from 2, the empty set from 3, the expression matrix | a Seurat object with `catcheR_status` per cell | your analysis |

`<chem>` is `tenx`, `pip` or `parse`; `<pfx>` the same.

**Why the order is forced.** An "empty" cell in step 3 is one that is *not* in
the called set, so without step 2 there is nothing to subtract from. Step 4
needs step 2 for the calls and step 3 only if you want the controls labelled
rather than left as `not_detected`.

**What is cheap to repeat.** Steps 2–4 act on the in-memory object, so trying
new thresholds costs seconds. Only step 1 reads the FASTQ.

---

## Where files land

```
experiment/
  scratch/                          <- inputs, plus the handoff copies below
  Output/
    tenxcatch/                      <- step 1, one folder per chemistry
      filtercatch_p15_u1_r2_f20/    <- step 2, one folder per threshold set
    pipcatch/
    parsecatch/
```

Step 1 nests under the chemistry name because the 13 plots that carry catcheR's
own file names (`barcode_distribution`, `UMIxUCI`, `true_UCIxcell`, …) have no
chemistry in the filename - deliberately, so they can be cited as in the
published figures. Without the nesting, a Parse run would silently replace a
10x run's figures. Step 2 adds a folder named for its thresholds, so a second
attempt sits beside the first rather than on top of it.

Everything downstream defaults `output_dir` to `x$output_dir`, so the nesting
propagates on its own; the templates take `OUT` from the object rather than
rebuilding the path. `subdir = FALSE` restores the old flat layout.

Two samples of the **same** chemistry in one `output_dir` still collide, since
the 13 names carry no sample either. The code warns, naming the sample already
present; use one `output_dir` per sample.

---

## Deliverable or handoff?

A `.csv` in `Output/` is something you read. A `.csv` the next step reads is an
*input*, and in catcheR's layout inputs live in `scratch/`. Three files are of
the second kind:

| file | consumed by |
|---|---|
| `complete_table_fin_<sample>.csv` | the legacy shell filtering step |
| `cells_with_at_least_1_UCI.csv` | `*_nocatch()` - the empty-cell selection subtracts it |
| `new_cell_names.csv` | your own renaming, downstream of the Seurat join |
| `<sample>_<chemistry>_input_files.csv` | you, when you need to know which file played which role |

The first carries the sample but not the chemistry; the other two carry neither.
Inside one chemistry's folder that is harmless. The moment two pipelines share a
`scratch/` it is wrong - and that is exactly what a multi-chemistry experiment
does.

So the copy is explicit and renamed:

```r
flt <- catcheR_tenx_filtercatch(res, ..., to.scratch = "experiment/scratch")
```

```
scratch/
  cells_with_at_least_1_UCI_tenxcatch_AB019a.csv
  complete_table_fin_tenxcatch_AB019a.csv
  new_cell_names_tenxcatch_AB019a.csv
```

The originals stay in `output_dir`; nothing is moved. The analysis templates
pass `to.scratch` by default (`to_scratch: true` in the YAML header) - set it to
`false` if you would rather keep `scratch/` for inputs only.

In the **native R** path the handoff between steps 1→2→3→4 is the object, not a
file, so a filename collision never breaks catcheR itself. It breaks *your*
downstream use, which is why the renaming exists. In the **legacy shell** path
the files are the handoff, and there the collision is structural.

---

## Knowing where a result came from

Step 1 writes `<sample>_<chemistry>_input_files.csv` beside its outputs and
stores the same table on the object as `$input_files`: one row per input, with
the role it played, the path actually opened, its size and its modification
time. It travels to `scratch/` with the other handoff files.

```
role                      exists  bytes
fastq read 1              TRUE    640634817
fastq read 2              TRUE    943163876
barcode-to-gene reference TRUE    1707
barcode-to-gene reference TRUE    1116
called-cell list          TRUE    32953
```

This is an input *record*, not a checksum and not evidence of biological
identity: it will not catch a correctly-named file containing the wrong sample.
What it does catch is the ordinary accident - a reference from the wrong library,
a FASTQ that was replaced since the run, an input that silently did not exist.

---

## The thresholds, exactly

Two exclude strictly, one includes. This is catcheR's own behaviour and is the
easiest thing in the package to get wrong:

| parameter | comparison | so |
|---|---|---|
| `UMI.count` | `UMIxUCI > UMI.count` | the highest count still **excluded**. `1` keeps pairs of 2+ UMIs; `0` keeps singletons |
| `percentage` | `UMIpercentagexUCI > percentage` | a UCI at exactly 15% does **not** pass `percentage = 15` |
| `min.cell.UMI` | `UMIxcellID >= min.cell.UMI` | **inclusive**: a cell with exactly 20 passes `min.cell.UMI = 20` |

`ratio` applies only once a cell is down to one surviving UCI: it must hold at
least `ratio` times the runner-up's UMIs, unless it holds 100% of the cell, in
which case there is no runner-up.

**The default is a starting point, and the function says so.** When
`min.cell.UMI` is not supplied, `catcheR_*_filtercatch()` reports which value it
used and writes two things beside its results: `catcheR_call_accuracy_<sample>.csv`,
which subsamples the deep cells and scores every outcome as correct, wrong or no
call, and `<chemistry>_yield_vs_error.png`, which puts the yield and the expected
wrong-call rate on the same axis against the floor. On AB019c that is 13.4%
expected wrong calls at a floor of 5, 8.8% at 10 and 1.1% at 20 - which is what
makes 20 the default rather than a guess. Supplying the argument explicitly
suppresses both, on the assumption that you have already chosen. Either way the
calls are identical: the panel is a measurement, not a filter.

**The floor keeps applying by default.** Both step 1 and step 2 default
`min.cell.UMI` to 20, so a re-filter applies the same floor unless you pass
another value; `0` disables it. The two steps share no state - the defaults
simply match, which also means that lowering the floor at extraction does not
lower it at re-filter time. The templates pass it explicitly at both steps so
the value in the YAML header is the one in force.

**On a read-counted chemistry there is no floor to set.** PIPseq's molecular
index is too short to deduplicate, so every count is reads, per-cell totals top
out around ten, and a floor tuned for a UMI-counted library excludes every cell.
`catcheR_threshold_sweep()` detects this and substitutes the nearest usable
floor, saying so in the panel subtitle.
