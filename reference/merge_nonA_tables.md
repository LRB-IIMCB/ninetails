# Merges ninetails tabular outputs (read classes and nonadenosine residue data) to produce one concise table.

Combines the `read_classes` and `nonadenosine_residues` data frames into
a single wide-format tibble with one row per read. The `prediction`
column from the residue data is spread into three columns
(`prediction_C`, `prediction_G`, `prediction_U`) via
[`spread_nonA_residues`](https://LRB-IIMCB.github.io/ninetails/reference/spread_nonA_residues.md),
and a CIGAR-like `nonA_residues` summary column is appended.

## Usage

``` r
merge_nonA_tables(class_data, residue_data, pass_only = TRUE)
```

## Arguments

- class_data:

  Data frame or tibble containing read_classes predictions produced by
  the ninetails pipeline. The `qc_tag` column may be character
  (Guppy/nanopolish: `"PASS"`, `"SUFFCLIP"`, etc.) or numeric (Dorado:
  MAPQ score).

- residue_data:

  Data frame or tibble containing non-A residue predictions produced by
  the ninetails pipeline. The `qc_tag` column type must match that of
  `class_data`.

- pass_only:

  Logical `[TRUE]`. Applies only when `qc_tag` is character
  (Guppy/nanopolish pipeline). If `TRUE`, only reads tagged as `"PASS"`
  are included. If `FALSE`, reads tagged as `"PASS"` or `"SUFFCLIP"` are
  included. **Ignored when `qc_tag` is numeric** (Dorado pipeline), in
  which case MAPQ \> 0 filtering is applied instead.

## Value

A tibble with summarised information from both ninetails outputs: all
columns from `class_data` plus `prediction_C`, `prediction_G`,
`prediction_U`, and `nonA_residues`.

## Details

Unclassified reads are excluded before merging. The quality-tag filter
behaviour depends on the pipeline that produced the input data:

**Guppy/nanopolish pipeline** (character `qc_tag`): The `pass_only`
parameter controls filtering. When `TRUE`, only reads tagged as `"PASS"`
are retained. When `FALSE`, reads tagged as `"PASS"` or `"SUFFCLIP"` are
retained.

**Dorado pipeline** (numeric `qc_tag`): The `qc_tag` column contains
MAPQ scores. Filtering is based on mapping quality: only reads with
`qc_tag > 0` are retained. The `pass_only` parameter is ignored for
numeric `qc_tag`.

After the full join, any `NA` values in numeric columns are replaced
with 0.

## Pipeline Compatibility

This function automatically detects whether data originates from the
legacy Guppy/nanopolish pipeline (character `qc_tag`) or the Dorado
pipeline (numeric `qc_tag`) and applies appropriate filtering logic.
Both `class_data` and `residue_data` must originate from the same
pipeline to ensure consistent `qc_tag` types during the merge operation.

## See also

[`spread_nonA_residues`](https://LRB-IIMCB.github.io/ninetails/reference/spread_nonA_residues.md)
for the reshaping step,
[`summarize_nonA`](https://LRB-IIMCB.github.io/ninetails/reference/summarize_nonA.md)
for transcript-level summaries from the merged table,
[`calculate_fisher`](https://LRB-IIMCB.github.io/ninetails/reference/calculate_fisher.md)
for statistical testing on the merged table,
[`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
and
[`read_residue_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
for loading the input data.

## Examples

``` r
if (FALSE) { # \dontrun{

# Guppy/nanopolish data (character qc_tag)
merged_tables <- ninetails::merge_nonA_tables(
  class_data = class_data,
  residue_data = residue_data,
  pass_only = TRUE)

# Dorado data (numeric qc_tag) - pass_only is ignored
merged_tables <- ninetails::merge_nonA_tables(
  class_data = dorado_class_data,
  residue_data = dorado_residue_data,
  pass_only = TRUE)

} # }
```
