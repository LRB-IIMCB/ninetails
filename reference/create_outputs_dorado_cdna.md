# Create Ninetails output tables for Dorado cDNA pipeline

This function extends create_outputs_dorado to handle cDNA data with
additional tail_type information. It produces the same output structure
as the DRS pipeline but preserves the tail_type column for read
orientation information.

## Usage

``` r
create_outputs_dorado_cdna(
  dorado_summary_dir,
  nonA_temp_dir,
  polya_chunks_dir,
  num_cores = 1,
  qc = TRUE
)
```

## Arguments

- dorado_summary_dir:

  Character string. Path to a directory containing Dorado summary files
  (.txt, .tsv, or .csv) with per-read poly(A) tail information and
  tail_type.

- nonA_temp_dir:

  Character string. Path to a directory containing non-adenosine
  prediction RDS files, generated from temporary models.

- polya_chunks_dir:

  Character string. Path to a directory containing poly(A) chunk RDS
  files used for position inference of predictions.

- num_cores:

  Integer. Number of cores to use for parallelized file loading and
  processing. Must be a positive integer. Default is 1.

- qc:

  Logical. Whether to apply quality control filtering of terminal
  predictions (removing predictions near the ends of poly(A) tails).
  Default is TRUE.

## Value

A named list with two data frames (identical to DRS output + tail_type):

- read_classes:

  Data frame with per-read classification results, including columns for
  read name, contig, poly(A) length, QC tag, class, comments, and
  tail_type.

- nonadenosine_residues:

  Data frame with per-chunk predictions of non-adenosine residues,
  including read name, contig, predicted base, estimated position within
  the poly(A) tail, poly(A) length, QC tag, and tail_type.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- create_outputs_dorado_cdna(
dorado_summary_dir = "data/dorado_summaries",
nonA_temp_dir = "data/nonA_predictions",
polya_chunks_dir = "data/polya_chunks",
num_cores = 4,
qc = TRUE
)

# Access read classifications with tail_type
head(results$read_classes)

# Access non-adenosine residues with tail_type
head(results$nonadenosine_residues)
} # }
```
