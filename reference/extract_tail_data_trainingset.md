# Extracts tail features of single RNA read from respective basecalled multi-fast5 file.

This is the training-set variant of the single-read extraction function.
It reads raw signal, basecaller moves and channel metadata from a
multi-Fast5 file, isolates the poly(A) tail region defined by
nanopolish, winsorizes the signal, downsamples by linear interpolation
(to 20% of original length), and computes pseudomoves via
[`filter_signal_by_threshold_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold_trainingset.md).

## Usage

``` r
extract_tail_data_trainingset(
  readname,
  polya_summary,
  workspace,
  basecall_group
)
```

## Arguments

- readname:

  Character string. Name (UUID) of the given read within the analyzed
  dataset.

- polya_summary:

  Data frame. The table containing data extracted from nanopolish and
  the sequencing summary (produced by
  [`extract_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md)).

- workspace:

  Character string. Full path of the directory containing basecalled
  multi-Fast5 files.

- basecall_group:

  Character string `["Basecall_1D_000"]`. Name of the level in the Fast5
  file hierarchy from which the data should be extracted.

## Value

A named list with four elements:

- fast5_filename:

  Character. Name of the source Fast5 file.

- tail_signal:

  Numeric vector. Winsorized and downsampled signal corresponding to the
  poly(A) tail region.

- tail_moves:

  Numeric vector. Downsampled basecaller moves for the tail region.

- tail_pseudomoves:

  Numeric vector. Pseudomove values (-1, 0, 1) indicating potential
  non-A deviations in the tail signal.

Always assign this returned list to a variable; printing the full list
to the console may crash the R session.

## Details

The function differs from its production counterpart
([`extract_tail_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data.md))
in that:

- The signal is downsampled (interpolated) to 20% of its original length
  to speed up downstream training-set preparation.

- Pseudomoves are computed using the training-set-specific threshold
  filter
  ([`filter_signal_by_threshold_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold_trainingset.md)).

- The returned list omits production-only fields and instead provides
  elements suitable for chunk splitting and GAF creation.

## See also

[`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md)
for the parallel wrapper that calls this function,
[`extract_tail_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data.md)
for the production counterpart,
[`winsorize_signal`](https://LRB-IIMCB.github.io/ninetails/reference/winsorize_signal.md)
for the signal clipping step,
[`filter_signal_by_threshold_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold_trainingset.md)
for pseudomove detection.

## Examples

``` r
if (FALSE) { # \dontrun{

extract_tail_data_trainingset(
  readname = 'abc123de-fg45-6789-0987-6543hijk2109',
  polya_summary = polya_summary_table,
  workspace = '/path/to/folder/containing/multifast5s',
  basecall_group = 'Basecall_1D_000')

} # }

```
