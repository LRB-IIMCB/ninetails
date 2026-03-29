# Merge polyA and polyT processing results for cDNA analysis

This function combines the results from separate polyA and polyT
processing paths into a unified output in standard ninetails format.
Updated to work with the simplified approach using tibbles instead of
files.

## Usage

``` r
merge_cdna_results(
  polya_results = NULL,
  polyt_results = NULL,
  unidentified_reads = NULL,
  save_dir,
  prefix = "",
  cli_log = message
)
```

## Arguments

- polya_results:

  List. Results from polyA processing (from process_polya_reads_cdna).
  Can be NULL if no polyA reads were found.

- polyt_results:

  List. Results from polyT processing (from process_polyt_reads_cdna).
  Can be NULL if no polyT reads were found.

- unidentified_reads:

  Data frame/tibble. Reads that could not be classified as polyA or
  polyT (with tail_type = "unidentified").

- save_dir:

  Character. Directory where merged results will be saved.

- prefix:

  Character. Optional prefix for output file names.

- cli_log:

  Function for logging messages and progress.

## Value

List containing merged cDNA results in standard ninetails format:

- read_classes:

  Data frame with all read classifications including tail_type

- nonadenosine_residues:

  Data frame with all predicted modifications including tail_type

## Examples

``` r
if (FALSE) { # \dontrun{
merged_results <- merge_cdna_results(
  polya_results = polya_output,
  polyt_results = polyt_output,
  unidentified_reads = unidentified_tibble,
  save_dir = "path/to/output/",
  prefix = "experiment1",
  cli_log = message
)
} # }
```
