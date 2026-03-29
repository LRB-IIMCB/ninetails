# Reads multiple ninetails nonadenosine_residues outputs at once.

Batch-loads any number of `nonadenosine_residues` prediction files with
a single invocation, attaching sample-level metadata from the provided
`samples_table`. Each file is loaded via
[`read_residue_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
and the results are combined into a single long-format tibble.

## Usage

``` r
read_residue_multiple(samples_table, ...)
```

## Arguments

- samples_table:

  Data frame or tibble containing sample metadata and file paths. Must
  have at least two columns:

  residue_path

  :   Character. Path to the `nonadenosine_residues` prediction file for
      each sample.

  sample_name

  :   Character/factor. Unique sample identifier.

  Additional metadata columns are preserved and propagated to the
  output.

- ...:

  Additional parameters passed to
  [`read_residue_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
  (currently accepts `sample_name`).

## Value

A tibble containing non-A residue data for all specified samples, with
metadata from `samples_table` stored as additional columns.

## Acknowledgements

Function based on `read_polya_multiple` from the NanoTail package by P.
Krawczyk (smaegol): <https://github.com/LRB-IIMCB/nanotail/>.

## See also

[`read_residue_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
for the per-file loader,
[`read_class_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_multiple.md)
for the analogous class data loader,
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for combining class and residue outputs.

## Examples

``` r
if (FALSE) { # \dontrun{

samples_table <- data.frame(
  residue_path = c(path_1, path_2, path_3, path_4),
  sample_name = c("wt_1", "mut_1", "wt_2", "mut_2"),
  group = c("wt", "mut", "wt", "mut"))

residues_data <- ninetails::read_residue_multiple(samples_table)

} # }
```
