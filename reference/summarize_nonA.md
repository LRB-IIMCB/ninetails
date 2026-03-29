# Produces summary table of non-A occurrences within an analyzed dataset.

Creates a per-transcript summary table with read counts, non-A residue
counts and hits, and poly(A) tail length statistics, grouped by
user-defined factors (e.g. sample, condition).

## Usage

``` r
summarize_nonA(
  merged_nonA_tables,
  summary_factors = c("group"),
  transcript_id_column = c("ensembl_transcript_id_short")
)
```

## Arguments

- merged_nonA_tables:

  Data frame or tibble. Output of
  [`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md).

- summary_factors:

  Character string or vector of strings. Column name(s) used for
  grouping (default: `"group"`).

- transcript_id_column:

  Character string. Column containing the transcript identifier
  (default: `"ensembl_transcript_id_short"`, as added during data
  pre-processing; can be changed by the user).

## Value

A tibble with one row per transcript per group, containing:

- polya_median:

  Numeric. Median poly(A) tail length for the transcript.

- polya_mean:

  Numeric. Mean poly(A) tail length for the transcript.

- counts_total:

  Integer. Total number of reads mapped to the transcript.

- counts_blank:

  Integer. Number of reads with no non-A residues.

- counts_nonA / hits_nonA:

  Integer. Reads with any non-A / total non-A hits.

- counts_C / hits_C:

  Integer. Reads with C / total C hits.

- counts_G / hits_G:

  Integer. Reads with G / total G hits.

- counts_U / hits_U:

  Integer. Reads with U / total U hits.

## Details

The distinction between **counts** and **hits**:

- **counts** — number of reads containing at least one occurrence of a
  given non-A residue type.

- **hits** — total number of occurrences of a given non-A residue type
  across all reads (a single read may contribute multiple hits).

## See also

[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for preparing the input,
[`calculate_fisher`](https://LRB-IIMCB.github.io/ninetails/reference/calculate_fisher.md)
for statistical testing on merged data,
[`annotate_with_biomart`](https://LRB-IIMCB.github.io/ninetails/reference/annotate_with_biomart.md)
for adding gene-level annotation.

## Examples

``` r
if (FALSE) { # \dontrun{

summarized <- ninetails::summarize_nonA(
  merged_nonA_tables = merged_nonA_tables,
  summary_factors = "group",
  transcript_id_column = "ensembl_transcript_id_short")

} # }
```
