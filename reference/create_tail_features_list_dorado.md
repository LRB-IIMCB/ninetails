# Creates a nested list of Dorado tail features (raw signal + pseudomoves).

This function processes raw poly(A) tail signal traces from Dorado and
computes pseudomoves using a threshold-based filter. Each read is stored
in a nested list structure, keyed by its read ID, with both the raw
signal and the computed pseudomove vector.

## Usage

``` r
create_tail_features_list_dorado(signal_list, num_cores)
```

## Arguments

- signal_list:

  list of numeric vectors. Each element must represent a raw poly(A)
  tail signal trace, with list names corresponding to read IDs.

- num_cores:

  numeric \[1\]. Number of physical cores to use in processing. Do not
  exceed 1 less than the total available cores on your system.

## Value

A nested list of tail features, organized by read IDs. Each read entry
contains:

- `tail_signal`: numeric vector of the raw poly(A) tail signal.

- `tail_pseudomoves`: numeric vector of pseudomove states computed by
  [`ninetails::filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md).

## Details

Unlike Guppy-based pipelines, Dorado does not provide moves directly in
BAM files, so pseudomoves are computed empirically from the raw signal
traces.

Parallel execution is handled with foreach and doSNOW, and a progress
bar is displayed during processing.

## Examples

``` r
if (FALSE) { # \dontrun{
features <- ninetails::create_tail_features_list_dorado(
  signal_list = signals,
  num_cores = 4
)
} # }
```
