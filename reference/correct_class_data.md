# Corrects the classification of reads contained in the class_data table.

Introduces two additional columns (`corr_class` and `corr_comments`) to
the `class_data` table, reclassifying reads based on the positional
quality flags computed by
[`correct_residue_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md).
Reads whose *all* non-A positions were flagged as ambiguous
(`qc_pos == "N"`) are downgraded from `"decorated"` to `"blank"`.

## Usage

``` r
correct_class_data(residue_data_edited, class_data)
```

## Arguments

- residue_data_edited:

  Data frame or tibble. Corrected non-A residue predictions produced by
  [`correct_residue_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md).
  Must contain columns `mode_pos`, `seg_err_quart`, and `qc_pos`.

- class_data:

  Data frame or tibble containing read_classes predictions from the
  ninetails pipeline.

## Value

A tibble identical to `class_data` with two additional columns:

- corr_class:

  Character. Corrected classification: original `class` value, or
  `"blank"` for artefact-only reads.

- corr_comments:

  Character. Corrected comment: original `comments` value, or `"MPU"`
  for reclassified reads.

## Details

The reclassification logic: for each read, if every non-A residue
position has `qc_pos == "N"`, the read is reclassified as `"blank"` with
comment `"MPU"` (to maintain compatibility with the tag system used in
plotting functions). Reads with at least one `qc_pos == "Y"` position
retain their original class and comment.

## Caution

It is recommended to use
[`reclassify_ninetails_data`](https://LRB-IIMCB.github.io/ninetails/reference/reclassify_ninetails_data.md)
for downstream analyses, as it wraps this function together with
[`correct_residue_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md)
and performs the necessary column renaming. If using this function
directly, rename `corr_class` and `corr_comments` to `class` and
`comments` before plotting.

## See also

[`correct_residue_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md)
for the preceding step,
[`reclassify_ninetails_data`](https://LRB-IIMCB.github.io/ninetails/reference/reclassify_ninetails_data.md)
for the high-level wrapper,
[`create_outputs`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)
for the pipeline that produces the input data.

## Examples

``` r
if (FALSE) { # \dontrun{

class_data_corrected <- ninetails::correct_class_data(
  residue_data_edited = residue_data_edited,
  class_data = results[[1]])

} # }
```
