# Aggregates nanopolish polya quality control information.

Returns read counts for each `qc_tag` category (e.g. `"PASS"`,
`"SUFFCLIP"`, `"NOREGION"`, etc.), optionally stratified by a
user-defined grouping variable.

## Usage

``` r
nanopolish_qc(class_data, grouping_factor = NA)
```

## Arguments

- class_data:

  Data frame or tibble containing the raw `class_data` output from the
  ninetails pipeline (not the merged output).

- grouping_factor:

  Character string (default `NA`). Variable used for grouping (e.g.
  `"sample_name"`).

## Value

A tibble with columns for the grouping variable (if provided), `qc_tag`,
and `n` (the count).

## Details

**Caution:** do not use the output of
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
as input, because that table is already filtered for quality. Use the
raw `class_data` output from the ninetails pipeline instead.

## Acknowledgements

This is the ninetails implementation of the
`get_nanopolish_processing_info` function originally written by P.
Krawczyk (smaegol) and incorporated in the NanoTail package:
<https://github.com/LRB-IIMCB/nanotail/blob/master/R/polya_stats.R>.
Variable names were adjusted to match the ninetails naming convention.

## See also

[`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
for loading class data,
[`count_class`](https://LRB-IIMCB.github.io/ninetails/reference/count_class.md)
for counting read classification labels.
