# Extract modification-centered signal fragments from a poly(A) tail

Finds areas in the poly(A) tail signal containing potential
non-adenosine residues and extracts 100-point signal fragments where the
potential modification is always at the center of a given extracted
fragment.

## Usage

``` r
split_tail_centered(readname, tail_feature_list)
```

## Arguments

- readname:

  Character string. Name of the given read (UUID) within the analyzed
  dataset.

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md).
  Must contain per-read entries with `tail_signal`, `tail_moves`, and
  `tail_pseudomoves`.

## Value

A nested list where each element represents one candidate modification
region. Each element is itself a list with:

- chunk_sequence:

  Numeric vector of raw signal values (length 100)

- chunk_start_pos:

  Integer. Start index of the chunk in the original signal

- chunk_end_pos:

  Integer. End index of the chunk in the original signal

- chunk_moves:

  Numeric vector. Basecaller moves for the chunk region

Element names follow the pattern `<readname>_<chunk_index>`. Positions
are indexed from the 3' end.

## Details

Candidate modification regions are identified based on two assumptions:
the presence of significant raw signal distortion (recorded as a
pseudomove by the thresholding algorithm) and the transition of state
(`move == 1`) recorded by Guppy. If only `move == 0` values are present
within a given signal chunk, then that chunk is dropped from the
analysis (the distortion is most likely caused by a sequencing artifact,
not a non-A residue itself).

If the data indicating the presence of modifications are near the signal
ends (3' or 5'), missing upstream or downstream data are imputed based
on the most frequent values in the entire signal.

The extraction procedure is as follows:

1.  Run-length encoding (RLE) of the pseudomove vector

2.  Filtering runs of pseudomoves with length \>= 5

3.  Centering a 100-point window on the midpoint of each qualifying run

4.  Imputing NAs at boundaries with draws from the 5 most frequent
    signal values

5.  Removing chunks where basecaller moves are all zero (likely
    artifacts)

## See also

[`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)
for preparing the input,
[`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md)
for batch segmentation,
[`split_tail_centered_dorado`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_dorado.md)
for the Dorado-specific version

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::split_tail_centered(
  readname = "1234-anexample-r3adn4m3",
  tail_feature_list = tail_feature_list
)

} # }
```
