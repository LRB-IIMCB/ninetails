# Plot panel characteristics of ninetails output

Generates a five-panel summary visualization describing read categories,
poly(A) tail properties, and non-A residue distributions derived from
the ninetails pipeline output. Those panel charts provide the most
comprehensive characterization of a given pool of reads (representing
particular transcript or set of transcripts, respectively).

## Usage

``` r
plot_panel_characteristics(
  input_residue_data,
  input_class_data = NULL,
  input_merged_nonA_tables_data = NULL,
  type = "default",
  max_length = 300,
  direction_5_prime = TRUE
)
```

## Arguments

- input_residue_data:

  Data frame containing non-A residue predictions.

- input_class_data:

  Optional data frame containing read classification output from
  ninetails. Mutually exclusive with `input_merged_nonA_tables_data`.

- input_merged_nonA_tables_data:

  Optional data frame returned by
  [`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md).
  Mutually exclusive with `input_class_data`.

- type:

  Character string. Either `"default"` or `"moderna"`. The `"moderna"`
  option marks the default UCUAG pentamer position (100 nt).

- max_length:

  Numeric. Maximum poly(A) tail length displayed in distribution panels.

- direction_5_prime:

  Logical. If `TRUE` (default), non-A positions are reported from the 5'
  end of the poly(A) tail. If `FALSE`, positions are recalculated
  relative to the 3' end.

## Value

A patchwork-assembled `ggplot` object containing panels A–E.

## Details

The function assembles panels A–E into a patchwork layout:

- **A** – Read categories (blank, non-A containing, total)

- **B** – Counts of reads containing C, G, or U residues

- **C** – Distribution of poly(A) tail lengths

- **D** – Normalized distribution of non-A positions

- **E** – Raw distribution of non-A positions

## Internal column filtering

This function internally subsets input data to a predefined set of
columns. Any additional columns present in the input data frames are
silently dropped. Therefore, all columns listed below must be present in
the supplied inputs.

## Required columns

### input_residue_data

Must contain at least the following columns:

- readname:

  Character. Unique read identifier.

- prediction:

  Character or factor. Predicted non-A residue (e.g., "C", "G", "U").

- est_nonA_pos:

  Numeric. Estimated position of the non-A residue within the poly(A)
  tail.

- polya_length:

  Numeric. Estimated poly(A) tail length.

These columns are used for:

- binning positions and tail lengths (Panels D and E),

- computing normalized residue frequencies,

- generating residue-level counts and labels.

### input_class_data (if supplied)

Must contain at least:

- readname:

  Character. Required for merging with residue data.

- group:

  Character or factor. Experimental group identifier.

This input is internally merged with `input_residue_data` using
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md).

### input_merged_nonA_tables_data (if supplied)

Must be the output of
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
and contain at least the following columns (used for tail distribution,
summarization and plotting):

- sample:

  Character. Sample identifier.

- group:

  Character or factor. Experimental group identifier.

- readname:

  Character. Unique read identifier.

- prediction:

  Character or factor. Non-A residue prediction.

- est_nonA_pos:

  Numeric. Estimated non-A position.

- polya_length:

  Numeric. Poly(A) tail length.

- class:

  Character. Read classification label.

- comments:

  Character. Additional classification comments.

- transcript:

  Character. Transcript name.

- ensembl_transcript_id_full:

  Character. Full Ensembl transcript ID.

- ensembl_transcript_id_short:

  Character. Short Ensembl transcript ID.

- prediction_C:

  Numeric or logical. Indicator/count of C residues.

- prediction_G:

  Numeric or logical. Indicator/count of G residues.

- prediction_U:

  Numeric or logical. Indicator/count of U residues.

- nonA_residues:

  Character. Encoded non-A residue information; `NA` indicates blank
  tail.

These columns are required because they are internally retained for:

- classification summaries (Panels A and B),

- tail length distributions (Panel C),

- grouping and normalization logic,

- detection of blank vs non-A containing reads.

## Single-group requirement

The input data must represent exactly one experimental group. If
multiple unique values are detected in the `group` column, the function
stops with an error.

Multiple samples within a single group are allowed.

## See also

[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md),
[`summarize_nonA`](https://LRB-IIMCB.github.io/ninetails/reference/summarize_nonA.md),
[`plot_tail_distribution`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_distribution.md)

## Examples

``` r
if (FALSE) { # \dontrun{
residue_data_wt <- residue_data[residue_data$group == "WT", ]

plot_panel_characteristics(
  input_residue_data = residue_data_wt,
  input_class_data = class_data_wt,
  type = "default",
  max_length = 100
)
} # }
```
