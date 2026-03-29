# Process poly(A) data split into multiple parts through the Guppy pipeline.

Iterates over a set of poly(A) data part files (produced by
[`split_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/split_polya_data.md)),
processing each sequentially via
[`process_polya_complete`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_complete.md),
saving per-part outputs, and merging all results into a single output
list.

## Usage

``` r
process_polya_parts(
  part_files,
  sequencing_summary,
  workspace,
  num_cores,
  basecall_group,
  pass_only,
  qc,
  save_dir,
  prefix,
  cli_log,
  ...
)
```

## Arguments

- part_files:

  Character vector. File paths to poly(A) data part files (as returned
  by
  [`split_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/split_polya_data.md)).

- sequencing_summary:

  Character string or data frame. Full path of the sequencing summary
  file, or an in-memory data frame.

- workspace:

  Character string. Full path of the directory containing basecalled
  multi-Fast5 files.

- num_cores:

  Numeric `[1]`. Number of physical cores.

- basecall_group:

  Character string `["Basecall_1D_000"]`. Fast5 hierarchy level for data
  extraction.

- pass_only:

  Logical `[TRUE]`. If `TRUE`, only `"PASS"` reads are included.

- qc:

  Logical `[TRUE]`. If `TRUE`, terminal artefact positions are labelled
  with `"-WARN"`.

- save_dir:

  Character string. Output directory path.

- prefix:

  Character string (optional). Output file name prefix.

- cli_log:

  Function. Logging closure defined in
  [`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md).

- ...:

  Additional arguments passed to
  [`process_polya_complete`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_complete.md).

## Value

A named list with two data frames (merged across all parts):

- read_classes:

  Per-read classification data.

- nonadenosine_residues:

  Detailed non-A residue positional data.

## Details

This function is called internally by
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
when the input exceeds the `part_size` threshold. It requires the
`cli_log` closure defined inside
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for formatted logging, and therefore should not be called directly by
the user.

For each part, intermediate results are saved as TSV files in a
dedicated subdirectory (`part_<i>_of_<n>`) within `save_dir`. After all
parts are processed, `read_classes` and `nonadenosine_residues` tables
are merged with `rbind`.

## See also

[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
which calls this function,
[`split_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/split_polya_data.md)
for the input splitting step,
[`process_polya_complete`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_complete.md)
for the per-part processing.
