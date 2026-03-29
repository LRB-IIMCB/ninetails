# Create list of poly(A) tail chunks centered on significant signal deviations

Extracts fragments of poly(A) tails of ONT RNA reads potentially
containing non-A nucleotides, along with their coordinates, and appends
the data to a nested list organized by read IDs.

## Usage

``` r
create_tail_chunk_list(tail_feature_list, num_cores)
```

## Arguments

- tail_feature_list:

  List object produced by
  [`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md).

- num_cores:

  Numeric. Number of physical cores to use in processing. Do not exceed
  1 less than the number of cores at your disposal.

## Value

A nested list containing the segmented tail data (chunks and
coordinates) organized by read IDs. Each read entry contains one or more
fragments, where each fragment is a list with:

- chunk_sequence:

  Numeric vector. Raw signal values (length 100)

- chunk_start_pos:

  Integer. Starting index of the chunk in the original signal

- chunk_end_pos:

  Integer. Ending index of the chunk in the original signal

## Details

Parallelization is handled with foreach and doSNOW. A progress bar is
displayed during processing. After extraction, empty list elements are
pruned and chunk moves are dropped to minimize memory footprint.

## See also

[`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)
for preparing the input,
[`split_tail_centered`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered.md)
for the per-read segmentation logic,
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
for downstream GAF transformation

## Examples

``` r
if (FALSE) { # \dontrun{

tcl <- ninetails::create_tail_chunk_list(
  tail_feature_list = tfl,
  num_cores = 3
)

} # }
```
