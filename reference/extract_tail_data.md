# Extract tail features of a single RNA read from a multi-Fast5 file

Extracts metadata and signal features of a single RNA read from a
multi-Fast5 file basecalled by Guppy. The tail signal, as delimited by
nanopolish polya function, is extracted, winsorized (to remove signal
cliffs), and downsampled to 20% of its original length to facilitate
further analysis. Pseudomoves are computed from the processed signal
using
[`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md).

## Usage

``` r
extract_tail_data(readname, polya_summary, workspace, basecall_group)
```

## Arguments

- readname:

  Character string. Name of the given read (UUID) within the analyzed
  dataset.

- polya_summary:

  Data frame. The table containing data extracted from nanopolish and
  sequencing summary, as produced by
  [`extract_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md).

- workspace:

  Character string. Full path of the directory containing the basecalled
  multi-Fast5 files.

- basecall_group:

  Character string. Name of the level in the Fast5 file hierarchy from
  which data should be extracted (e.g., `"Basecall_1D_000"`).

## Value

A named list containing per-read tail features:

- fast5_filename:

  Character. Name of the source Fast5 file

- tail_signal:

  Numeric vector. Winsorized and downsampled poly(A) tail signal

- tail_moves:

  Numeric vector. Downsampled basecaller moves for the tail region

- tail_pseudomoves:

  Numeric vector. Pseudomove states computed from the tail signal ({-1,
  0, 1})

Always assign the returned list to a variable. Printing the full output
to the console may crash your R session.

## Details

The function performs the following operations:

1.  Reads raw signal from the multi-Fast5 file via rhdf5

2.  Retrieves basecaller moves and stride information

3.  Extracts the poly(A) tail region of the signal based on nanopolish
    coordinates

4.  Applies winsorization to the tail signal
    ([`winsorize_signal`](https://LRB-IIMCB.github.io/ninetails/reference/winsorize_signal.md))

5.  Downsamples both signal and moves to 20% of original length via
    linear interpolation

6.  Computes pseudomoves using
    [`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)

## See also

[`extract_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md)
for preparing the `polya_summary` input,
[`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)
for batch extraction,
[`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)
for pseudomove computation

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::extract_tail_data(
  readname = 'abc123de-fg45-6789-0987-6543hijk2109',
  polya_summary = polya_summary_table,
  workspace = '/path/to/folder/containing/multifast5s',
  basecall_group = 'Basecall_1D_000'
)

} # }
```
