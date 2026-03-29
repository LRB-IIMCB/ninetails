# Process and split Dorado summary file into smaller parts

Splits a Dorado summary file or data frame (in-memory file) into
multiple smaller files for downstream analysis. This helps to avoid
memory overflow and data loss during processing.

## Usage

``` r
process_dorado_summary(dorado_summary, save_dir, part_size, cli_log)
```

## Arguments

- dorado_summary:

  Character path to Dorado summary file, or a data frame containing
  summary information.

- save_dir:

  Character path to directory where split summary files will be saved.

- part_size:

  Integer. Number of reads per file part when splitting the summary
  file.

- cli_log:

  Function for logging messages and progress.

## Value

Character vector containing paths to the split summary files.

## Examples

``` r
if (FALSE) { # \dontrun{
summary_files <- process_dorado_summary(
  dorado_summary = "path/to/summary.txt",
  save_dir = "path/to/output/",
  part_size = 40000,
  cli_log = message
)
} # }
```
