# Substitute short zero-gaps surrounded by nonzero pseudomoves

Smooths the pseudomove vector by filling in isolated zero-gaps (length
\< 2) flanked by identical nonzero values. This avoids redundancy
introduced by the z-score thresholding algorithm when the signal is
jagged, so a single modification region is reported instead of multiple
fragmented segments.

## Usage

``` r
substitute_gaps(pseudomoves)
```

## Arguments

- pseudomoves:

  Numeric vector. Pseudomove vector produced by the z-score filtering
  algorithm
  ([`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)),
  corresponding to the tail region of the read of interest as delimited
  by nanopolish polya function.

## Value

Numeric vector. Adjusted pseudomoves with short gaps filled.

## See also

[`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)
where this function is called as a post-processing step

## Examples

``` r
if (FALSE) { # \dontrun{

substitute_gaps(pseudomoves = pseudomoves_vector)

} # }
```
