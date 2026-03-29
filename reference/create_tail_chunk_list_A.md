# Creates list of tail chunks containing only A nucleotides.

Extracts overlapping fragments of poly(A) tail signals that do *not*
contain non-A nucleotides and organises them in a nested list keyed by
read ID.

## Usage

``` r
create_tail_chunk_list_A(tail_feature_list, num_cores)
```

## Arguments

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md).

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

## Value

A named nested list organised by read IDs, where each element is a list
of numeric vectors (each of length 100) representing the overlapping
signal fragments.

## Details

This training-set variant is designed for A-only reference data
preparation. Unlike
[`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
(which centres chunks on modifications), this function uses
[`split_with_overlaps`](https://LRB-IIMCB.github.io/ninetails/reference/split_with_overlaps.md)
with `segment = 100` and `overlap = 50`, effectively performing data
augmentation by producing overlapping windows across the entire tail
signal.

## See also

[`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md)
for the preceding pipeline step,
[`split_with_overlaps`](https://LRB-IIMCB.github.io/ninetails/reference/split_with_overlaps.md)
for the per-read segmentation logic,
[`create_gaf_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list_A.md)
for the next pipeline step,
[`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
for the non-A variant.

## Examples

``` r
if (FALSE) { # \dontrun{

create_tail_chunk_list_A(
  tail_feature_list = tail_feature_list,
  num_cores = 2)

} # }
```
