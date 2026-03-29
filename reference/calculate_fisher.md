# Perform Fisher's exact test per transcript with BH p-value adjustment

Runs Fisher's exact test for each transcript in ninetails output data
and applies the Benjamini-Hochberg (BH) step-up procedure to control the
false discovery rate. This is a high-level wrapper for
[`fisher.test`](https://rdrr.io/r/stats/fisher.test.html) and
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) with additional data
wrangling features.

## Usage

``` r
calculate_fisher(
  ninetails_data,
  transcript_id_column = "ensembl_transcript_id_short",
  min_reads = 0,
  min_nonA_reads = 0,
  grouping_factor = "sample_name",
  condition1 = NA,
  condition2 = NA,
  alpha = 0.05,
  base = "C",
  ...
)
```

## Arguments

- ninetails_data:

  Data frame. The output of
  [`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
  (merged tabular output containing read classification and non-A
  position data).

- transcript_id_column:

  Character string. Column with transcript ID data. Default:
  `"ensembl_transcript_id_short"`; can be changed by the user.

- min_reads:

  Numeric. Minimum number of reads representing a given transcript to
  include it in the analysis. Default: 0. Keep in mind that including
  many transcripts with low coverage increases the risk of rejecting the
  true null hypothesis (Benjamini-Hochberg procedure).

- min_nonA_reads:

  Numeric. Minimum number of reads containing non-adenosine residues
  (sum of C, G, U) per transcript to include it in the analysis. This
  prevents considering too many observations as non-significant during
  p-value adjustment. Non-A-containing reads are typically a small
  fraction of the total pool, so additional filtering can provide more
  meaningful results. Default: 0.

- grouping_factor:

  Character string. Name of the grouping variable. Default:
  `"sample_name"`.

- condition1:

  Character string. First level of `grouping_factor` to use for
  comparison. Required when `grouping_factor` has more than 2 levels.

- condition2:

  Character string. Second level of `grouping_factor` to use for
  comparison. Required when `grouping_factor` has more than 2 levels.

- alpha:

  Numeric. Significance threshold for FDR. Default: 0.05.

- base:

  Character string. Letter representing the non-A nucleotide for which
  statistics are computed. Accepted values: `"C"`, `"G"`, `"U"`, or
  `"all"` (for all non-A residues combined). Default: `"C"`.

- ...:

  Additional parameters passed to
  [`nonA_fisher`](https://LRB-IIMCB.github.io/ninetails/reference/nonA_fisher.md)
  (under development).

## Value

A tibble with per-transcript test results, sorted by adjusted p-value.
Columns include:

- \<transcript_id_column\>:

  Transcript identifier

- p.value:

  Numeric. Raw p-value from Fisher's exact test

- stats_code:

  Character. Descriptive status code (see `stat_codes_list` for full
  definitions)

- padj:

  Numeric. BH-adjusted p-value

- significance:

  Character. `"FDR<alpha"` if significant, `"NotSig"` otherwise

## Acknowledgements

Inspired by the Nanotail package written and maintained by Pawel
Krawczyk (smaegol):
<https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R>. Many
thanks to the developer of the original source code.

## See also

[`nonA_fisher`](https://LRB-IIMCB.github.io/ninetails/reference/nonA_fisher.md)
for the per-transcript Fisher's exact test,
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for preparing the input,
[`summarize_nonA`](https://LRB-IIMCB.github.io/ninetails/reference/summarize_nonA.md)
for contingency table computation

## Examples

``` r
if (FALSE) { # \dontrun{

test <- ninetails::calculate_fisher(
  ninetails_data = merged_nonA_tables,
  transcript_id_column = "ensembl_transcript_id_short",
  min_reads = 100,
  min_nonA_reads = 10,
  grouping_factor = "sample_name",
  condition1 = "WT",
  condition2 = "KO",
  alpha = 0.05,
  base = "C"
)

} # }
```
