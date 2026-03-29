# Perform Fisher's exact test on a single transcript in ninetails output

Runs Fisher's exact test for testing the null hypothesis of independence
of rows and columns in a 2x2 contingency table representing a given
transcript in ninetails output data. This is a wrapper for
[`fisher.test`](https://rdrr.io/r/stats/fisher.test.html) with
additional data wrangling features.

## Usage

``` r
nonA_fisher(
  ninetails_data,
  grouping_factor,
  base,
  min_reads = 0,
  transcript_id_column = NA
)
```

## Arguments

- ninetails_data:

  Data frame. The output of
  [`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
  (merged tabular output containing read classification and non-A
  position data).

- grouping_factor:

  Character string. The name of the factor variable defining
  groups/conditions (must have exactly 2 levels).

- base:

  Character string. Letter representing the non-A nucleotide for which
  statistics are computed. Accepted values: `"C"`, `"G"`, `"U"`, or
  `"all"` (for all non-A residues combined). Default: `"C"`.

- min_reads:

  Numeric. Minimum number of reads representing a given transcript to
  include it in the analysis. Default: 0.

- transcript_id_column:

  Character string. Name of the column in which transcript identifiers
  are stored. Default: `NA`.

## Value

A tibble with results for the given transcript, containing:

- p.value:

  Numeric. Raw p-value from Fisher's exact test, or `NA` if conditions
  were not met.

- stats_code:

  Character. Status code describing whether the test conditions were met
  (see `stat_codes_list` for definitions).

## Details

The function is suitable only for pairwise comparisons (2x2 contingency
tables), where two conditions (e.g. WT vs KO) are compared at once. The
user can set a cutoff number of reads required for the analysis.

This function is intended to be called internally by
[`calculate_fisher`](https://LRB-IIMCB.github.io/ninetails/reference/calculate_fisher.md)
and is not typically used standalone.

## Acknowledgements

Inspired by the Nanotail package written and maintained by Pawel
Krawczyk (smaegol):
<https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R>. Many
thanks to the developer of the original source code.

## See also

[`calculate_fisher`](https://LRB-IIMCB.github.io/ninetails/reference/calculate_fisher.md)
for the per-transcript wrapper with multiple-testing correction,
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for preparing the input,
[`summarize_nonA`](https://LRB-IIMCB.github.io/ninetails/reference/summarize_nonA.md)
for contingency table computation

## Examples

``` r
if (FALSE) { # \dontrun{

test <- ninetails::nonA_fisher(
  ninetails_data = merged_nonA_tables,
  grouping_factor = "sample_name",
  base = "C",
  min_reads = 100
)

} # }
```
