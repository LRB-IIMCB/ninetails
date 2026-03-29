# Create ninetails output tables (Guppy legacy pipeline)

Integrates nanopolish poly(A) tail data, CNN predictions, and positional
information to generate two main outputs: per-read classifications and
non-adenosine residue predictions with estimated positions along the
poly(A) tail.

## Usage

``` r
create_outputs(
  tail_feature_list,
  tail_chunk_list,
  nanopolish,
  predicted_list,
  num_cores,
  pass_only = TRUE,
  qc = TRUE
)
```

## Arguments

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md).

- tail_chunk_list:

  List object produced by
  [`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md).

- nanopolish:

  Character string or data frame. Full path of the `.tsv` file produced
  by nanopolish polya function, or an in-memory data frame.

- predicted_list:

  List object produced by
  [`predict_gaf_classes`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md).

- num_cores:

  Numeric. Number of physical cores to use in processing. Do not exceed
  1 less than the number of cores at your disposal.

- pass_only:

  Logical. If `TRUE` (default), only reads tagged by nanopolish as
  `"PASS"` are taken into consideration. If `FALSE`, reads tagged as
  `"PASS"` and `"SUFFCLIP"` are both included.

- qc:

  Logical. If `TRUE` (default), quality control of the output
  predictions is performed. Reads with non-A residue positions in
  terminal nucleotides (\< 2 nt from either end of the tail) are labeled
  with `"-WARN"` suffixes, as these are most likely artifacts inherited
  from nanopolish segmentation. It is then up to the user whether to
  include or discard such reads from downstream analysis.

## Value

A named list with two data frames:

- read_classes:

  Data frame with per-read classification results. Columns: `readname`,
  `contig`, `polya_length`, `qc_tag`, `class`, and `comments`.

- nonadenosine_residues:

  Data frame with per-chunk predictions of non-adenosine residues for
  decorated reads only. Columns: `readname`, `contig`, `prediction`,
  `est_nonA_pos`, `polya_length`, `qc_tag`. When `qc = TRUE`, prediction
  values may carry a `"-WARN"` suffix for terminal positions.

## Details

The function implements a complete read accounting system that
classifies all reads from the nanopolish output into biologically
meaningful categories based on both quality control metrics and
modification detection results.

## Read Classification System

Reads are assigned into three main categories with specific comment
codes:

**decorated** - Reads with detected non-adenosine modifications

- YAY: Move transition present, nonA residue detected

**blank** - Reads without detected modifications

- MAU: Move transition absent, nonA residue undetected

- MPU: Move transition present, nonA residue undetected

**unclassified** - Reads that failed quality control

- IRL: Insufficient read length (poly(A) \< 10 nt)

- QCF: Nanopolish QC failed

- NIN: Not included in the analysis (`pass_only = TRUE`)

## Position Estimation

Non-adenosine residue positions are estimated using the formula:

`est_nonA_pos = polya_length - ((polya_length * centr_signal_pos) / signal_length)`

Positions are reported as distance from the 3' end of the tail (position
1 = most 3' nucleotide).

## See also

[`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)
for feature extraction,
[`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md)
for chunk segmentation,
[`predict_gaf_classes`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md)
for CNN classification,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for the complete pipeline wrapper

## Examples

``` r
if (FALSE) { # \dontrun{

create_outputs(
  tail_feature_list = tail_feature_list,
  tail_chunk_list = tail_chunk_list,
  nanopolish = '/path/to/nanopolish_output.tsv',
  predicted_list = predicted_list,
  num_cores = 2,
  pass_only = TRUE,
  qc = TRUE
)

} # }

```
