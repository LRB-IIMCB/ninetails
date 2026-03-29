# Process polyT reads using ninetails pipeline

This function processes reads that have been classified as
polyT-containing through the cDNA classification pipeline. Updated to
work with tibble input instead of separate files. It applies the
ninetails analysis pipeline to identify non-thymidine residues within
polyT tails. Currently uses the same model as polyA processing until a
polyT-specific model is trained.

## Usage

``` r
process_polyt_reads_cdna(
  polyt_sequences,
  signal_files,
  num_cores = 1,
  qc = TRUE,
  save_dir,
  prefix = "",
  cli_log = message
)
```

## Arguments

- polyt_sequences:

  Data frame/tibble. PolyT-classified sequences with read_id and
  tail_type columns.

- signal_files:

  Character vector. Paths to signal files containing poly(A) tail
  signals extracted from POD5 files (will be filtered for polyT reads).

- num_cores:

  Integer. Number of CPU cores to use for parallel processing.

- qc:

  Logical. Whether to apply quality control filtering.

- save_dir:

  Character. Directory where processing results will be saved.

- prefix:

  Character. Optional prefix for output file names.

- cli_log:

  Function for logging messages and progress.

## Value

List containing polyT processing results:

- read_classes:

  Data frame with read classification results (with tail_type preserved)

- nonadenosine_residues:

  Data frame with predicted modifications (with tail_type preserved)

- processing_stats:

  Summary statistics for polyT processing

- tail_type:

  Character indicating this is "polyT" data

## Examples

``` r
if (FALSE) { # \dontrun{
polyt_results <- process_polyt_reads_cdna(
  polyt_sequences = polyt_tibble,
  signal_files = c("polya_signal_part1.rds"),
  num_cores = 4,
  qc = TRUE,
  save_dir = "path/to/output/",
  prefix = "experiment1",
  cli_log = message
)
} # }
```
