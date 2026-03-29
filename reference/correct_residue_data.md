# Marks uncertain positions of non-A residues in ninetails output data.

Annotates each non-A residue position with a quality flag (`qc_pos`)
indicating whether the position is likely genuine (`"Y"`) or a potential
nanopolish segmentation artefact (`"N"`). The assessment is based on
quantile thresholds of poly(A) tail length, modal non-A position per
transcript, and optional species-specific whitelists of transcripts with
A-rich 3' UTRs.

## Usage

``` r
correct_residue_data(
  class_data,
  residue_data,
  grouping_factor = NULL,
  transcript_column,
  ref = NULL
)
```

## Arguments

- class_data:

  Data frame or tibble containing read_classes predictions from the
  ninetails pipeline.

- residue_data:

  Data frame or tibble containing non-A residue predictions from the
  ninetails pipeline.

- grouping_factor:

  Character string or `NULL` (default). A grouping variable (e.g.
  `"sample_name"`, `"group"`).

- transcript_column:

  Character string. Name of the column containing transcript identifiers
  (e.g. `"ensembl_transcript_id_short"`).

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

  A custom character vector of transcript IDs may also be provided. Must
  be consistent with the content of `transcript_column`.

## Value

A tibble based on `residue_data` with the following additional columns:

- mode_pos:

  Integer. Most frequent non-A position reported for the transcript.

- mode_len:

  Integer. Most frequent tail length reported for the transcript.

- seg_err_quart:

  Numeric. 0.05 quantile of tail length for the transcript.

- qc_pos:

  Character. Quality flag: `"Y"` for likely genuine, `"N"` for
  ambiguous.

- pos_err_quart:

  Numeric. 0.05 quantile of non-A position for the transcript and
  prediction type.

- count_nonA:

  Integer. Number of non-A-containing reads for the transcript.

- count:

  Integer. Total number of reads for the transcript.

## Details

Nanopolish segmentation can misidentify nucleotides from A-rich 3' UTR
regions as part of the poly(A) tail. For such transcripts, a peak of
non-A positions accumulates near the transcript body. This function
flags those positions as ambiguous using a combination of:

- The modal position of non-A residues per transcript.

- The 0.05 quantile of poly(A) tail lengths per transcript (segmentation
  error boundary).

- The 0.05 quantile of non-A positions per transcript and prediction
  type.

- Species-specific whitelists of transcripts with hybrid tails (3' UTRs
  with \> 80% A in the last 20 positions).

## See also

[`correct_class_data`](https://LRB-IIMCB.github.io/ninetails/reference/correct_class_data.md)
for the companion function that reclassifies reads based on this output,
[`reclassify_ninetails_data`](https://LRB-IIMCB.github.io/ninetails/reference/reclassify_ninetails_data.md)
for the high-level wrapper,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
and
[`create_outputs`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)
for the pipeline that produces the input data.

## Examples

``` r
if (FALSE) { # \dontrun{

residue_data_edited <- ninetails::correct_residue_data(
  class_data = results[[1]],
  residue_data = results[[2]],
  transcript_column = "contig")

} # }
```
