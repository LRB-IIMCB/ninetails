# Draws an entire squiggle for given read from POD5 file.

Creates segmented plot of raw or rescaled ONT RNA signal from Dorado
basecalled POD5 files. A standalone function; depends solely on Dorado
summary output and POD5 files.

## Usage

``` r
plot_squiggle_pod5(
  readname,
  dorado_summary,
  workspace,
  rescale = FALSE,
  residue_data = NULL,
  nonA_flank = 250
)
```

## Arguments

- readname:

  Character string. Name of the given read (read_id) within the analyzed
  dataset.

- dorado_summary:

  Character string or data frame. Either the full path to the Dorado
  summary file (.txt), or a pre-loaded data frame containing at minimum
  the columns: read_id, poly_tail_start, poly_tail_end, filename.

- workspace:

  Character string. Full path of the directory containing POD5 files.

- rescale:

  Logical \[TRUE/FALSE\]. If TRUE, the signal will be rescaled to
  picoamps (pA) per second (s). If FALSE, raw signal per position will
  be plotted. Default is FALSE.

- residue_data:

  Character string, data frame, or NULL. Either the full path to the
  nonadenosine_residues output file, or a pre-loaded data frame.
  Required columns: read_id (or readname), prediction (C/G/U),
  est_nonA_pos, polya_length. If NULL (default), no overlay is drawn.

- nonA_flank:

  Numeric. Number of raw signal positions to highlight on each side of
  the estimated non-A center. Default 250, matching the approximate
  extent of the 100-point interpolated chunk (~500 raw positions).

## Value

ggplot2 object with squiggle plot depicting nanopore read signal.

## Details

The output plot includes an entire squiggle corresponding to the given
ONT read. Vertical lines mark the 5' (red) and 3' (navy blue) termini of
polyA tail according to the Dorado poly(A) estimation. In order to
maintain readability of the graph (and to avoid plotting high cliffs -
e.g. jets of the signal caused by a sudden surge of current in the
sensor) the signal is winsorized.

When `residue_data` is provided, semi-transparent rectangles highlight
estimated non-adenosine residue positions within the poly(A) tail. Each
rectangle spans `+/- nonA_flank` raw signal positions around the
estimated modification center, matching the approximate extent of the
100-point interpolated chunk analyzed by the CNN.

## See also

[`plot_squiggle_fast5`](https://LRB-IIMCB.github.io/ninetails/reference/plot_squiggle_fast5.md)
for the fast5 equivalent,
[`plot_tail_range_pod5`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_range_pod5.md)
for zoomed tail region plot

## Examples

``` r
if (FALSE) { # \dontrun{

# Basic usage (no non-A overlay)
plot <- ninetails::plot_squiggle_pod5(
  readname = "0e8e52dc-3a71-4c33-9a00-e1209ba4d2e9",
  dorado_summary = system.file('extdata', 'test_data', 'pod5_DRS',
                               'aligned_summary.txt',
                               package = 'ninetails'),
  workspace = system.file('extdata', 'test_data', 'pod5_DRS',
                          package = 'ninetails'),
  rescale = FALSE)

# With non-A overlay
plot <- ninetails::plot_squiggle_pod5(
  readname = "0e8e52dc-3a71-4c33-9a00-e1209ba4d2e9",
  dorado_summary = "dorado_summary.txt",
  workspace = "/path/to/pod5/",
  residue_data = "nonadenosine_residues.txt")

print(plot)

} # }
```
