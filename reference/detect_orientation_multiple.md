# Classify multiple cDNA read orientations using Dorado-style poly tail detection

This function analyzes Nanopore cDNA sequences using the Dorado-style
approach for poly tail type detection. It uses edit distance-based
matching of SSP and VNP primer sequences at both ends of reads, testing
forward and reverse orientations to determine if reads contain poly(A)
or poly(T) tails. The method uses sliding window matching and requires
both score and separation thresholds to be met. The tail_type column is
added to original sequence files and all sequences are returned as a
single tibble.

## Usage

``` r
detect_orientation_multiple(sequence_files, num_cores = 1, cli_log = message)
```

## Arguments

- sequence_files:

  Character vector. Paths to sequence files (TSV format) containing
  read_id and sequence columns.

- num_cores:

  Integer. Number of CPU cores to use for parallel processing.

- cli_log:

  Function for logging messages and progress.

## Value

Data frame/tibble containing all sequences with added tail_type column
("polyA", "polyT", or "unidentified") based on Dorado-style edit
distance matching of SSP and VNP primers with score and separation
validation. The tail_type column is also written back to the original
sequence files for persistence.

## Examples

``` r
if (FALSE) { # \dontrun{
classified_sequences <- detect_orientation_multiple(
  sequence_files = c("sequences_part1.tsv", "sequences_part2.tsv"),
  num_cores = 4,
  cli_log = message
)
} # }
```
