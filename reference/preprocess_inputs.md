# Preprocess Dorado inputs for ninetails analysis (no BAM processing)

This function prepares Dorado summary files and extracts poly(A) signals
from POD5 files for downstream analysis. It splits large summary files
into manageable parts and extracts corresponding poly(A) signals.

## Usage

``` r
preprocess_inputs(
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

- dorado_summary:

  Character or data frame. Path to Dorado summary file or data frame.

- pod5_dir:

  Character. Directory containing pod5 files.

- num_cores:

  Integer. Number of CPU cores to use.

- qc:

  Logical. Whether to perform quality control filtering. When `TRUE`,
  applies the following stringent filters to ensure high-quality data:

  Alignment direction

  :   Removes unmapped reads (`alignment_direction == "*"`)

  Mapping quality

  :   Removes reads with zero mapping quality (`alignment_mapq == 0`)

  Poly(A) coordinates

  :   Removes reads with invalid start positions
      (`poly_tail_start == 0`)

  Tail length

  :   Removes reads with tails shorter than 10 nucleotides
      (`poly_tail_length < 10`)

  These filters remove low-quality reads, misaligned sequences, and
  likely artifacts from pore clogging or basecalling errors. Set to
  `FALSE` only for preliminary analysis or when using pre-filtered data.

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

- polya_signal_files:

  Paths to extracted poly(A) signal files

## Quality Control Details

The quality control filtering (when `qc = TRUE`) is implemented via
[`filter_dorado_summary`](https://LRB-IIMCB.github.io/ninetails/reference/filter_dorado_summary.md)
and performs the following operations:

- **Unmapped read removal**: Reads with `alignment_direction = "*"`
  indicate failed alignment to the reference and are excluded.

- **Low mapping quality removal**: Reads with `alignment_mapq = 0` have
  ambiguous or poor-quality alignments and are excluded.

- **Invalid coordinate removal**: Reads with `poly_tail_start = 0`
  indicate failed poly(A) detection, often from pore clogging artifacts.

- **Short tail removal**: Reads with `poly_tail_length < 10` nucleotides
  lack sufficient poly(A) sequence for reliable modification detection.

These criteria ensure that only high-confidence reads with valid poly(A)
tails and proper genomic alignment are processed, significantly
improving downstream classification accuracy and reducing false positive
rates.

## Performance Note

Filtering typically removes 25-30% of raw reads from standard Dorado
output, but this percentage may vary based on sequencing quality, pore
condition, and alignment parameters. The function logs the exact number
of reads removed for quality control tracking.

## Examples

``` r
if (FALSE) { # \dontrun{
processed_files <- preprocess_inputs(
  dorado_summary = "path/to/summary.txt",
  pod5_dir = "path/to/pod5/",
  num_cores = 2,
  qc = TRUE,
  save_dir = "path/to/output/",
  prefix = "experiment1",
  part_size = 40000,
  cli_log = message
)
} # }
```
