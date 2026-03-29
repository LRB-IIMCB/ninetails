# Process a single (unsplit) poly(A) data file through the Guppy pipeline.

Runs the full ninetails analysis on a small-to-moderate poly(A) data set
(up to `part_size` rows). Performs format validation, feature
extraction, tail segmentation, GAF computation, CNN classification, and
output assembly.

## Usage

``` r
process_polya_complete(
  polya_data,
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

- polya_data:

  Character string or data frame. Full path of the poly(A) length file
  (`.tsv`), or an in-memory data frame.

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
  [`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
  for formatted console and log file output.

- ...:

  Additional arguments (currently unused).

## Value

A named list with two data frames:

- read_classes:

  Per-read classification data.

- nonadenosine_residues:

  Detailed non-A residue positional data.

## Details

This function is called internally by
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
(for inputs within the `part_size` limit) and by
[`process_polya_parts`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_parts.md)
(for each part of a split input). It requires the `cli_log` closure
defined inside
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for formatted logging, and therefore should not be called directly by
the user.

## See also

[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
which calls this function,
[`process_polya_parts`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_parts.md)
for the split-input variant,
[`check_polya_length_filetype`](https://LRB-IIMCB.github.io/ninetails/reference/check_polya_length_filetype.md)
for format detection,
[`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md),
[`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md),
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md),
[`predict_gaf_classes`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md),
[`create_outputs`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)
for the individual pipeline steps.
