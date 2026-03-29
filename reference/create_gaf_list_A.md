# Produces list of GAFs containing exclusively A-nucleotides for neural network training.

Converts A-only signal chunks into Gramian Angular Field (GAF) matrices
in parallel. Each chunk is transformed via
[`combine_gafs`](https://LRB-IIMCB.github.io/ninetails/reference/combine_gafs.md)
and the resulting matrices are collected in a flat named list suitable
for direct input to the CNN training pipeline.

## Usage

``` r
create_gaf_list_A(tail_chunk_list, num_cores)
```

## Arguments

- tail_chunk_list:

  List object produced by
  [`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md).

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

## Value

A named list of GAF matrices, where each element corresponds to one
signal chunk and is named `<readname>_<index>`. Always assign this
returned list to a variable; printing the full list to the console may
crash the R session.

## Details

This training-set variant is designed to work on signal fragments that
contain only A residues (produced by
[`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)).
The naming convention for output elements is `<readname>_<chunk_index>`.

## See also

[`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)
for the preceding pipeline step,
[`combine_gafs`](https://LRB-IIMCB.github.io/ninetails/reference/combine_gafs.md)
for the GAF transformation,
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
for the production counterpart,
[`prepare_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/prepare_trainingset.md)
for the top-level wrapper.

## Examples

``` r
if (FALSE) { # \dontrun{

create_gaf_list_A(
  tail_chunk_list = tail_chunk_list,
  num_cores = 10)

} # }
```
