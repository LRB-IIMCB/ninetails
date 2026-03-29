# Preprocess Dorado inputs for ninetails cDNA analysis

This function prepares inputs for the cDNA pipeline by processing BAM
files, Dorado summary files, and extracting both basecalled sequences
and poly(A) signals from POD5 files. Unlike the DRS pipeline, this
function additionally extracts basecalled sequences from BAM files for
read orientation classification.

## Usage

``` r
preprocess_inputs_cdna(
  bam_file,
  dorado_summary,
  pod5_dir,
  num_cores,
  qc,
  save_dir,
  prefix,
  part_size,
  cli_log
)
```

## Arguments

- bam_file:

  Character string. Path to BAM file containing aligned cDNA reads with
  basecalled sequences.

- dorado_summary:

  Character or data frame. Path to Dorado summary file or data frame.

- pod5_dir:

  Character. Directory containing pod5 files.

- num_cores:

  Integer. Number of CPU cores to use.

- qc:

  Logical. Whether to perform quality control.

- save_dir:

  Character. Directory where output files will be saved.

- prefix:

  Character. Prefix to add to output file names (optional).

- part_size:

  Integer. Number of reads to process in each chunk.

- cli_log:

  Function for logging messages and progress.

## Value

List containing paths to processed files:

- summary_files:

  Paths to split summary files

- bam_files:

  Paths to split BAM files

- sequence_files:

  Paths to extracted sequence files

- polya_signal_files:

  Paths to extracted poly(A) signal files

## Examples

``` r
if (FALSE) { # \dontrun{
processed_files <- preprocess_inputs_cdna(
  bam_file = "path/to/aligned.bam",
  dorado_summary = "path/to/summary.txt",
  pod5_dir = "path/to/pod5/",
  num_cores = 4,
  qc = TRUE,
  save_dir = "path/to/output/",
  prefix = "experiment1",
  part_size = 40000,
  cli_log = message
)
} # }
```
