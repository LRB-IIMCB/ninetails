# Reclassifies ambiguous non-A residues to mitigate potential errors inherited from nanopolish segmentation.

High-level wrapper that combines
[`correct_residue_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md)
and
[`correct_class_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_class_data.md)
into a single call, producing cleaned class and residue data frames
ready for downstream analysis and visualisation.

## Usage

``` r
reclassify_ninetails_data(
  residue_data,
  class_data,
  grouping_factor = NULL,
  transcript_column,
  ref = NULL
)
```

## Arguments

- residue_data:

  Data frame or tibble containing non-A residue predictions from the
  ninetails pipeline.

- class_data:

  Data frame or tibble containing read_classes predictions from the
  ninetails pipeline.

- grouping_factor:

  Character string or `NULL` (default). A grouping variable (e.g.
  `"sample_name"`).

- transcript_column:

  Character string. Name of the column containing transcript identifiers
  (e.g. `"contig"`, `"ensembl_transcript_id_short"`).

- ref:

  Character string, character vector, or `NULL` (default). Whitelist of
  transcripts with hybrid tails. Built-in options:

  `"athaliana"`

  :   *Arabidopsis thaliana*

  `"hsapiens"`

  :   *Homo sapiens*

  `"mmusculus"`

  :   *Mus musculus*

  `"scerevisiae"`

  :   *Saccharomyces cerevisiae*

  `"celegans"`

  :   *Caenorhabditis elegans*

  `"tbrucei"`

  :   *Trypanosoma brucei*

  A custom character vector may also be provided. Must be consistent
  with the content of `transcript_column`. Using a whitelist is optional
  but allows retrieval of more true positive data.

## Value

A named list with two data frames:

- class_data:

  Data frame. Corrected read classifications with `class` and `comments`
  columns updated. Compatible with plotting functions.

- residue_data:

  Data frame. Filtered non-A residue predictions with ambiguous
  positions removed. Intermediate QC columns are dropped.

## Details

Nanopolish segmentation can misidentify nucleotides from A-rich 3' UTR
regions as part of the poly(A) tail. When tail boundaries are recognised
incorrectly, non-A positions accumulate near the 3' end of the
transcript, significantly affecting analysis results. This function
flags and removes those ambiguous positions and reclassifies affected
reads accordingly.

The procedure:

1.  [`correct_residue_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md)
    annotates each non-A position with a quality flag (`qc_pos`).

2.  [`correct_class_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_class_data.md)
    reclassifies reads whose *all* non-A positions are flagged as
    ambiguous.

3.  Ambiguous positions (`qc_pos == "N"`) are dropped from the residue
    table.

4.  Corrected columns (`corr_class`, `corr_comments`) are renamed to
    `class` and `comments`.

## Caution

Reads containing only ambiguous non-A positions are reclassified as
`"blank"` in the `class` column, and their `comments` are changed from
`"YAY"` to `"MPU"`.

## See also

[`correct_residue_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md)
and
[`correct_class_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_class_data.md)
for the underlying steps,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
and
[`create_outputs`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)
for the pipeline that produces the input data.

## Examples

``` r
if (FALSE) { # \dontrun{

rec_results <- ninetails::reclassify_ninetails_data(
  residue_data = results[[2]],
  class_data = results[[1]],
  transcript_column = "contig")

} # }
```
