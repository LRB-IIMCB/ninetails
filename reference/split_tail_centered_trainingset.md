# Extracts decoration-centered fragments of poly(A) tail signal along with positional coordinates.

Splits a single read's poly(A) tail signal into fixed-length (100
data-point) chunks, each centered on a potential non-A decoration
detected from the pseudomove vector. Chunks whose start position falls
before index 1 are left-padded with values randomly sampled from the
three most frequent signal values.

## Usage

``` r
split_tail_centered_trainingset(readname, tail_feature_list)
```

## Arguments

- readname:

  Character string. Name (UUID) of the given read.

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md).

## Value

A named nested list where each element corresponds to one chunk and
contains:

- chunk_sequence:

  Numeric vector (length 100). The signal fragment centered on the
  potential modification.

- chunk_start_pos:

  Integer. Start index of the chunk within the full tail signal (may be
  negative for left-padded chunks).

- chunk_end_pos:

  Integer. End index of the chunk within the full tail signal.

- pseudomoves:

  Numeric vector (length 100). Recomputed pseudomoves for the chunk.
  Coordinates are from the 3' end.

## Details

This training-set variant includes an additional `pseudomoves` element
in each chunk sublist, making the output suitable for supervised
training/validation data preparation. Pseudomoves for each extracted
chunk are recomputed by calling
[`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)
on the chunk sequence.

The centering procedure:

1.  Runs RLE on the pseudomove vector.

2.  Selects runs of length \>= 4 with non-zero values (empirical
    modification threshold).

3.  Centers a 100-element window on the midpoint of each selected run.

Chunk names follow the convention `<readname>_<index>`, where `index` is
the sequential position of the modification within the read (numbered
from the 3' end).

## See also

[`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
for the parallel wrapper that calls this function,
[`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)
for pseudomove recomputation on individual chunks,
[`split_tail_centered`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered.md)
for the production counterpart.

## Examples

``` r
if (FALSE) { # \dontrun{

split_tail_centered_trainingset(
  readname = "1234-anexample-r3adn4m3",
  tail_feature_list = tail_feature_list)

} # }
```
