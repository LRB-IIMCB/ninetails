# Create list of Gramian Angular Field matrices from tail chunks

Computes Gramian Angular Field (GAF) representations for all signal
chunks in a tail chunk list. Each chunk is transformed into a combined
GASF + GADF array via
[`combine_gafs`](https://LRB-IIMCB.github.io/ninetails/reference/combine_gafs.md).

## Usage

``` r
create_gaf_list(tail_chunk_list, num_cores)
```

## Arguments

- tail_chunk_list:

  List object produced by
  [`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md).

- num_cores:

  Numeric. Number of physical cores to use in processing. Do not exceed
  1 less than the number of cores at your disposal.

## Value

A named list of GAF arrays (100, 100, 2) organized by read `ID_index`.
Always assign the returned list to a variable. Printing the full output
to the console may crash your R session.

## Details

Parallelization is handled with foreach and doSNOW, and a progress bar
is displayed during processing.

## See also

[`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md)
for preparing the input,
[`combine_gafs`](https://LRB-IIMCB.github.io/ninetails/reference/combine_gafs.md)
for the per-chunk GAF transformation,
[`predict_gaf_classes`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md)
for downstream CNN classification

## Examples

``` r
if (FALSE) { # \dontrun{

gl <- ninetails::create_gaf_list(
  tail_chunk_list = tcl,
  num_cores = 2
)

} # }
```
