# Counts non-A residues found in a nonadenosine_residues data frame produced by the ninetails pipeline.

Tabulates the `prediction` column of the `nonadenosine_residues` data
frame, counting occurrences of each non-A nucleotide type (C, G, U).
Counts represent total hits across all reads, *not* per-read summaries
(a single read may contribute multiple hits).

## Usage

``` r
count_residues(residue_data, grouping_factor = NA)
```

## Arguments

- residue_data:

  Data frame or tibble containing non-A residue predictions produced by
  the ninetails pipeline.

- grouping_factor:

  Character string (default `NA`). Name of a column in `residue_data` to
  use as a grouping variable (e.g. `"sample_name"`).

## Value

A tibble with columns for the grouping variable (if provided),
`prediction` (nucleotide type), and `n` (the count).

## See also

[`count_class`](https://LRB-IIMCB.github.io/ninetails/reference/count_class.md)
for counting read-level classes,
[`read_residue_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
and
[`read_residue_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_multiple.md)
for loading residue data,
[`summarize_nonA`](https://LRB-IIMCB.github.io/ninetails/reference/summarize_nonA.md)
for transcript-level summaries.

## Examples

``` r
if (FALSE) { # \dontrun{

residue_counted <- ninetails::count_residues(
  residue_data = out[[2]],
  grouping_factor = NA)

} # }
```
