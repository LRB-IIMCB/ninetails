# Extracts decoration-centered fragments of poly(A) tails for all reads and appends positional data to a nested list.

Parallel wrapper around
[`split_tail_centered_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_trainingset.md).
For every read in the feature list, it extracts 100-element signal
chunks centered on potential non-A modifications and organises them in a
nested list keyed by read ID.

## Usage

``` r
create_tail_chunk_list_trainingset(tail_feature_list, num_cores)
```

## Arguments

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md).

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

## Value

A named nested list organised by read IDs, where each element is a list
of chunk sublists as returned by
[`split_tail_centered_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_trainingset.md).

## Details

This training-set variant is intended for preparing training and
validation datasets. Each chunk sublist contains four fields:
`chunk_sequence`, `chunk_start_pos`, `chunk_end_pos`, and `pseudomoves`
(see
[`split_tail_centered_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_trainingset.md)).

## See also

[`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md)
for the preceding pipeline step,
[`split_tail_centered_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_trainingset.md)
for the per-read extraction logic,
[`filter_nonA_chunks_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md)
for the next pipeline step,
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
for GAF conversion downstream.

## Examples

``` r
if (FALSE) { # \dontrun{

create_tail_chunk_list_trainingset(
  tail_feature_list = tail_feature_list,
  num_cores = 2)

} # }
```
