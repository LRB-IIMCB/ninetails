# Counts reads by number of non-A occurrence instances.

Categorises reads into three abundance bins based on how many separate
non-A residue occurrences (instances) were detected per read: `"single"`
(1), `"two"` (2), or `"more"` (\>= 3). Returns the count of reads in
each bin, optionally stratified by a grouping variable.

## Usage

``` r
count_nonA_abundance(residue_data, grouping_factor = NA)
```

## Arguments

- residue_data:

  Data frame or tibble containing non-A residue predictions from the
  ninetails pipeline.

- grouping_factor:

  Character string (default `NA`). Name of a column to use as a grouping
  variable (e.g. `"sample_name"`).

## Value

A tibble with columns for the grouping variable (if provided),
`instances` (one of `"single"`, `"two"`, `"more"`), and `count` (the
number of reads in each bin).

## See also

[`count_residues`](https://LRB-IIMCB.github.io/ninetails/reference/count_residues.md)
for counting total residue hits,
[`count_class`](https://LRB-IIMCB.github.io/ninetails/reference/count_class.md)
for counting read-level classes,
[`read_residue_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
for loading residue data.

## Examples

``` r
if (FALSE) { # \dontrun{

nonA_abundance <- ninetails::count_nonA_abundance(
  residue_data = residue_data,
  grouping_factor = "sample_name")

} # }
```
