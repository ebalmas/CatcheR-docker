# The Parse Evercode pipeline — `catcheR_parsecatch()`

**Template:** `inst/rmd/catcheR_parsecatch_analysis_EB.Rmd`
**Status:** knitted end to end on AB019c (21.5 M read pairs, two sublibraries) with this release's code - 8.8 min.

---

## What the reads look like

Parse is the mirror of PIPseq: cassette on **read 1**, cell identity on
**read 2**.

```
read 1   [1 nt] [rc(reference) 23] [rc(UCI) 6] [rc(barcode) 8] [vector constant 20] ...
read 2   [UMI 10] bc3(8) [linker A 12] bc2(8) [linker B 12] bc1(8) (poly-T)
```

Three things follow:

* **The cell barcode is three ligation rounds**, each drawn from a published
  well whitelist, separated by two fixed 12 nt linkers. Each round is corrected
  against its whitelist and converted to a **well index**, so the emitted cell
  ID is split-pipe's own `bc_wells` form: `WW_WW_WW__<sublibrary>`.
* **Read 2 opens with a real 10 nt UMI**, so events are deduplicated by
  molecule. Per-cell counts are molecules, and the saturation analysis is
  meaningful.
* **One leading base precedes the reference on read 1.** The shared classifier
  takes a `lead` offset for exactly this.

## The barcode definition file

`catcheR_parsecatch()` needs split-pipe's `barcode_data.csv`. It is not a plain
CSV: `#` comment lines carry the kit and chemistry and also separate the
rounds, and the single `bci,sequence,...` header sits after the first of them,
so `read.csv()` fails on it. The round each row belongs to is read from the
`bc_round` column.

One detail that silently breaks a hand-rolled whitelist: **round 1 holds two
sequences per well** — one oligo-dT and one random hexamer. Both must map to
the same well index. Deriving the whitelist from abundance instead of from this
file captures only one of the two classes.

## Sublibraries

Pass `fastq.read1` / `fastq.read2` as vectors and name them with
`sublibraries`. Each is streamed separately, the sublibrary tag is appended to
the cell ID, and the panels are faceted per sublibrary plus a pooled panel.
Which FASTQ corresponds to which sublibrary should be confirmed from
well-triplet overlap against `cell_metadata.csv.gz`, not assumed from file
naming order.

## Arguments that change the answer

| argument | default | what it does |
|---|---|---|
| `barcode.data` | `"barcode_data.csv"` | the whitelist source, as above. |
| `called.cells` | `NULL` | split-pipe's `cell_metadata.csv.gz`. The honest denominator. |
| `min.cell.UMI` | 20 | minimum cassette UMIs for a cell to be analysed. The most consequential parameter on Parse data. |
| `percentage` | 15 | minimum share of a cell's UMIs a UCI must hold. |
| `ratio` | 2 | dominant-to-second ratio once a cell is down to one surviving UCI. |
| `UMI.count` | 1 | minimum UMIs on a (cell, UCI) pair, strict `>`. |
| `max.edit.distance` | 1 | correction radius per round. |

An arithmetic effect worth knowing: a cell holding fewer than `100 / percentage`
UMIs has every single-UMI UCI clearing the percentage cut, so several pass and
the cell is rejected as ambiguous. At `percentage = 15` that is any cell under
about 7 UMIs; `UMI.count = 1` removes most of it.

## What AB019c showed

Of split-pipe's 4,040 called cells: 2,857 (70.7%) resolved to exactly one
clone, 203 (5.0%) to two, 28 (0.7%) to three or more, 949 (23.5%) had no
dominant UCI, 3 had no cassette UMI at all. Median depth 90 molecules per
called cell.

**More sequencing would not help this library.** An exact rarefaction over
817,674 molecules carrying 14.5 M reads showed marginal return falling from 298
to 73 new calls per million reads; barcodes with ≥20 UMIs already reach 4,239,
i.e. above the 4,040 cells split-pipe called, so the extra ones are ambient
rather than new cells. The ceiling is cells loaded and ambient cassette load,
not depth. Both are bench-side.

**Where to set the UMI floor.** Yield alone argues for no floor at all, which
is why the package ships a subsampling accuracy test instead
(`catcheR_call_accuracy()`): it scores each outcome as correct, wrong, or
absent, and a wrong call is worse than no call because it adds a cross-clone
assignment. On AB019c the expected wrong-call share flattens at a floor of 20
(2,531 cells, 1.9% expected wrong). The curve is **not monotonic** at the far
left — 5.7% at 2 UMIs, 10.0% at 3, peaking at 13.4% at 5 — because at 2–3 UMIs
almost no cell gets any call, so the wrong-call share is small for the wrong
reason.

## Where the outputs land

`output_dir` is nested one level: everything goes to `<output_dir>/parsecatch/`,
so another chemistry writing into the same `Output/` cannot overwrite the plots
that carry catcheR's own file names. Each re-filter adds a tagged subfolder of
its own (`filtercatch_p15_u1_r2_f20/`). Pass `subdir = FALSE` for a flat layout.

## Known limits

* Rounds are corrected against the supplied whitelist to edit distance 1;
  ambiguous neighbours are dropped rather than assigned.
* The ambient floor is large on enrichment libraries. Always quote the rate
  over split-pipe's called cells, never over observed barcodes — on this
  library the two differ by more than tenfold.
