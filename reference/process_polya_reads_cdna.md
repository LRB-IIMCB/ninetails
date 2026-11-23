# Process poly(A) reads with the cDNA poly(A) model

Classifies the tail signals of reads sequenced in the poly(A)
orientation (tail at the 3' end of the read) with the dedicated cDNA
model and creates the ninetails output tables for them.

## Usage

``` r
process_polya_reads_cdna(
  polya_sequences,
  signal_files,
  num_cores = 1,
  qc = TRUE,
  save_dir,
  prefix = "",
  cli_log = message
)
```

## Arguments

- polya_sequences:

  Data frame/tibble. Poly(A) reads with their tail coordinates:
  `read_id`, `tail_start`, `tail_end`, `polya_length` and optionally
  `contig` and `qc_tag` (see `read_tail_coordinates_cdna`).

- signal_files:

  Character vector. Paths to the `.rds` files with raw flanked tail
  signals (see `extract_tail_signals_cdna`).

- num_cores:

  Integer. Number of CPU cores (currently unused; the classification is
  parallelised by TensorFlow).

- qc:

  Logical. Whether to apply quality control filtering.

- save_dir:

  Character. Output directory. Intermediate files are written to its
  `polya_processing_dir` (or `polyt_processing_dir`) subdirectory.

- prefix:

  Character. Optional prefix for intermediate file names.

- cli_log:

  Function for logging messages and progress.

## Value

List containing poly(A) processing results:

- read_classes:

  Data frame with read classification results

- nonadenosine_residues:

  Data frame with predicted non-adenosines

- processing_stats:

  Summary statistics for poly(A) processing

- tail_type:

  Character indicating this is "polyA" data

## Details

Signals are classified with `predict_nonA_cdna`, one signal file at a
time. The calls of each file are written to disk, so memory usage does
not grow with the size of the sequencing run. The tables are then built
from these files by
[`create_outputs_dorado_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado_cdna.md).
Reads with invalid coordinates, tails shorter than 10 nt or (with
`qc = TRUE`) a Lachesis QC tag other than `PASS` are not classified.

## Examples

``` r
if (FALSE) { # \dontrun{
polya_results <- ninetails::process_polya_reads_cdna(
  polya_sequences = polya_coordinates,
  signal_files = c("polya_signal_part1.rds"),
  num_cores = 4,
  qc = TRUE,
  save_dir = "path/to/output/",
  prefix = "experiment1",
  cli_log = message
)
} # }
```
