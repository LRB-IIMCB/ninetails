# Counts read classes found in a read_classes data frame produced by the ninetails pipeline.

Tabulates the prediction information contained in the `read_classes`
data frame. Counts can be computed at two levels of granularity
(detailed or crude) and optionally stratified by a grouping variable.

## Usage

``` r
count_class(class_data, grouping_factor = NA, detailed = TRUE)
```

## Arguments

- class_data:

  Data frame or tibble containing read_classes predictions produced by
  the ninetails pipeline.

- grouping_factor:

  Character string (default `NA`). Name of a column in `class_data` to
  use as a grouping variable (e.g. `"sample_name"`).

- detailed:

  Logical `[TRUE]`. If `TRUE`, counts are provided based on the
  `comments` column (fine-grained). If `FALSE`, counts are provided
  based on the `class` column (crude: decorated / blank / unclassified
  only).

## Value

A tibble with columns for the grouping variable (if provided), the
classification label (`comments` or `class`), and `n` (the count).

## Details

When `detailed = TRUE`, counts are based on the `comments` column, which
carries fine-grained labels such as `"YAY"` (tail with non-A residues
detected). When `detailed = FALSE`, counts are based on the `class`
column with three broad categories: `"decorated"`, `"blank"`, and
`"unclassified"`.

## See also

[`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
and
[`read_class_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_multiple.md)
for loading class data,
[`count_residues`](https://LRB-IIMCB.github.io/ninetails/reference/count_residues.md)
for the analogous residue-level counts,
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for combining class and residue data.

## Examples

``` r
if (FALSE) { # \dontrun{

class_counted <- ninetails::count_class(
  class_data = out[[1]],
  grouping_factor = NA,
  detailed = TRUE)

} # }
```
