# Process Dorado poly(A) signal files for non-A prediction and tail chunk extraction

This function takes a set of Dorado poly(A) signal files (RDS format),
extracts features, segments tails, generates Gramian Angular Fields
(GAFs), applies classification models to predict non-A residues, and
saves both prediction results and tail chunk lists as RDS files in
specified output directories. Progress is logged using a cli logging
function. This function is not intended to be used outside the pipeline
wrapper.

## Usage

``` r
process_dorado_signal_files(
  polya_signal_files,
  nonA_temp_dir,
  polya_chunks_dir,
  num_cores,
  cli_log
)
```

## Arguments

- polya_signal_files:

  Character vector of paths to poly(A) signal files (RDS format)
  generated from Dorado outputs.

- nonA_temp_dir:

  Character path to the directory where non-A prediction results will be
  written. Directory will be created if it does not exist.

- polya_chunks_dir:

  Character path to the directory where tail chunk list RDS files will
  be written. Directory will be created if it does not exist.

- num_cores:

  Integer number of CPU cores to use for parallel feature extraction and
  downstream computations.

- cli_log:

  Function used for logging messages.

## Value

An (invisible) list summarizing processing results for each input file.
Each element contains:

- `success` (logical): whether the file was processed successfully

- `file` (character): path to the non-A prediction RDS file if
  successful

- `chunks_file` (character): path to the tail chunk list RDS file if
  successful

- `error` (character): error message if processing failed

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
signal_files <- list.files("signals/", pattern = "polya_signal.*\\.rds$", full.names = TRUE)
results <- process_dorado_signal_files(
  polya_signal_files = signal_files,
  nonA_temp_dir = "nonA_predictions/",
  polya_chunks_dir = "tail_chunks/",
  num_cores = 4,
  cli_log = function(msg, level, ...) message(sprintf("[%s] %s", level, msg))
)
} # }
```
