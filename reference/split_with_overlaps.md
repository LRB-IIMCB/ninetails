# Splits signal to overlapping fragments of equal length.

Divides a single read's poly(A) tail signal into fixed-length,
overlapping segments for data augmentation during training-set
preparation. If the signal is not evenly divisible by the segment
length, trailing `NA` values are imputed with values randomly sampled
from the 10 most frequent signal values.

## Usage

``` r
split_with_overlaps(readname, tail_feature_list, segment, overlap)
```

## Arguments

- readname:

  Character string. Name (UUID) of the given ONT signal.

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md).

- segment:

  Numeric `[1]`. Length (in data points) of each chunk to be created.

- overlap:

  Numeric `[1]`. Length (in data points) of the overlap between
  consecutive chunks.

## Value

A list of numeric vectors, each of length `segment`, representing the
overlapping signal fragments. Trailing `NA` values are replaced by
imputed values.

## Details

This function is used exclusively in the A-only training pipeline
([`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md))
to produce overlapping windows from reads that do *not* contain
modification signals. The overlap acts as data augmentation, multiplying
the number of training examples per read.

## See also

[`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)
which calls this function,
[`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md)
for the preceding step,
[`split_tail_centered_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_trainingset.md)
for the modification-centered splitting alternative.

## Examples

``` r
if (FALSE) { # \dontrun{

split_with_overlaps(
  readname = '12fdcb3-ewfd543-34552-1234ddta345',
  tail_feature_list = tail_feature_list,
  segment = 300,
  overlap = 100)

} # }
```
