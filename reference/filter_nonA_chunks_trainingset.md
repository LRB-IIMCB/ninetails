# Filters read chunks containing non-adenosine nucleotides of interest for neural network training-set preparation.

Designed for use with synthetic spike-in data containing a single type
of non-A residue (G, C, or U in the context of 3'-homopolymer A). The
function retains only chunks whose pseudomove vectors contain a
consecutive run of the specified `value` (peak or valley) of length \>=
4 and whose start position is non-negative.

## Usage

``` r
filter_nonA_chunks_trainingset(tail_chunk_list, value, num_cores)
```

## Arguments

- tail_chunk_list:

  List object produced by
  [`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md).

- value:

  Numeric `[1]`. Controls whether valleys (C/U) or peaks (G) are
  retained:

  -1

  :   Retain chunks containing valleys (C or U nucleotides).

  1

  :   Retain chunks containing peaks (G nucleotide).

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

## Value

A named list of filtered chunk sublists, where names correspond to read
IDs. Empty reads (no chunks passing the filter) are removed.

## Details

The filtering takes advantage of the characteristic pseudomove patterns
produced by the z-score signal filter:

- G nucleotides produce a **peak** (pseudomove = +1).

- C and U nucleotides produce a **valley** (pseudomove = -1).

The function does *not* distinguish between C and U (both produce
valleys); classification is handled downstream by the CNN. Therefore it
is essential that training datasets for C and U differ in the transcript
bodies they map to and/or are not sequenced in a single run.

**Important:** this function is *not* suitable for preparing the A-only
reference set. Use
[`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md)
for that purpose.

Before proceeding, visual inspection of at least some filtered signals
is advisable. It may be necessary to manually adjust hardcoded
parameters (contact the developer/maintainer for details).

## See also

[`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
for the preceding pipeline step,
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
for GAF conversion downstream,
[`prepare_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/prepare_trainingset.md)
for the top-level wrapper.

## Examples

``` r
if (FALSE) { # \dontrun{

# filtering G residue:
filter_nonA_chunks_trainingset(
  tail_chunk_list = list_object_with_tail_chunks,
  value = 1, num_cores = 2)

# filtering C/U (user must know which type of residue is dealt with):
filter_nonA_chunks_trainingset(
  tail_chunk_list = list_object_with_tail_chunks,
  value = -1, num_cores = 2)

} # }
```
