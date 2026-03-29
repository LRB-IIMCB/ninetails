# Extract data from BAM file for cDNA analysis

This function extracts various informations from BAM files for reads
that are present in the corresponding Dorado summary file. It returns a
data frame with read IDs and their corresponding additional info, such
as basecalled sequences, poly(A) lengths, poly(A) coordinates etc.,
which can then be used in downstream analyses.

## Usage

``` r
extract_data_from_bam(
  bam_file,
  summary_file,
  seq_only = TRUE,
  cli_log = message
)
```

## Arguments

- bam_file:

  Character string. Path to BAM file containing basecalled sequences.

- summary_file:

  Character string. Path to corresponding Dorado summary file containing
  read IDs to extract.

- seq_only:

  logical \[TRUE\]. When `TRUE`, only minimal information is extracted,
  including read id, pod5 file name, and basecalled sequence. If
  `FALSE`, a more comprehensive information is extracted, including
  poly(A) tail length, coordinates etc.

- cli_log:

  Function for logging messages and progress.

## Value

Data frame with columns:

- read_id:

  Character. Read identifier

- sequence:

  Character. Basecalled sequence

- sequence_length:

  Integer. Length of the sequence

- mapping_quality:

  Integer. Mapping quality score

## Examples

``` r
if (FALSE) { # \dontrun{
sequences <- extract_data_from_bam(
  bam_file = "aligned_reads.bam",
  summary_file = "dorado_summary.txt",
  seq_only = TRUE,
  cli_log = message
)
} # }
```
