# Reads multiple ninetails read_classes outputs at once.

Batch-loads any number of `read_classes` prediction files with a single
invocation, attaching sample-level metadata from the provided
`samples_table`. Each file is loaded via
[`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
and the results are combined into a single long-format tibble.

## Usage

``` r
read_class_multiple(samples_table, ...)
```

## Arguments

- samples_table:

  Data frame or tibble containing sample metadata and file paths. Must
  have at least two columns:

  class_path

  :   Character. Path to the `read_classes` prediction file for each
      sample.

  sample_name

  :   Character/factor. Unique sample identifier.

  Additional metadata columns (e.g. `group`, `condition`) are preserved
  and propagated to the output.

- ...:

  Additional parameters passed to
  [`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
  (currently accepts `sample_name`).

## Value

A tibble containing read_classes data for all specified samples, with
metadata from `samples_table` stored as additional columns.

## Acknowledgements

Function based on `read_polya_multiple` from the NanoTail package by P.
Krawczyk (smaegol): <https://github.com/LRB-IIMCB/nanotail/>.

## See also

[`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
for the per-file loader,
[`read_residue_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_multiple.md)
for the analogous residue data loader,
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for combining class and residue outputs.

## Examples

``` r
if (FALSE) { # \dontrun{

samples_table <- data.frame(
  class_path = c(path_1, path_2, path_3, path_4),
  sample_name = c("wt_1", "mut_1", "wt_2", "mut_2"),
  group = c("wt", "mut", "wt", "mut"))

classes_data <- ninetails::read_class_multiple(samples_table)

} # }
```
