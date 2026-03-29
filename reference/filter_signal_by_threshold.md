# Detect outliers (peaks and valleys) in ONT signal using z-scores

Detects areas in the poly(A) tail signal where values significantly
deviate from those typical of an adenosine homopolymer. The function
produces a vector of pseudomoves with values in the range {-1, 0, 1},
where -1 corresponds to signals significantly below the local mean
(potential C/U), 1 corresponds to signals significantly above the local
mean (potential G), and 0 corresponds to typical adenosine homopolymer
values.

## Usage

``` r
filter_signal_by_threshold(signal)
```

## Arguments

- signal:

  Numeric vector. An ONT read fragment corresponding to the poly(A) tail
  region of the read of interest, as delimited by nanopolish polya
  function. Fragments are stored in
  `tail_feature_list[[1]][[readname]][[2]]` produced by
  [`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md).

## Value

Numeric vector of pseudomoves corresponding to the analyzed tail region,
with the same length as the input signal. Values are integers in the
range {-1, 0, 1}:

- -1:

  Signal valleys (potential C/U nucleotides)

- 0:

  Normal signal (likely A nucleotides)

- 1:

  Signal peaks (potential G nucleotides)

## Details

The pseudomoves vector allows more accurate calibration of nucleotide
positions of potential non-adenosine residues than the moves produced by
the Guppy basecaller.

The algorithm implements a sliding-window approach based on adaptive
z-scores:

1.  A calibration phase (100 synthetic data points drawn from the most
    frequent signal values) establishes baseline signal characteristics

2.  Rolling mean and standard deviation are computed over a 100-point
    adaptive window

3.  Outliers are detected when signal deviates \> 3.5 standard
    deviations from the local mean

4.  Direction of deviation determines the pseudomove value

Terminal positions (first 5) are masked to prevent false positives at
the adapter-tail boundary. This is a temporary safeguard until the
resegmentation model is refined.

## References

Based on: Brakel, J.P.G. van (2014). "Robust peak detection algorithm
using z-scores". Stack Overflow. Available at:
<https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362>
(version: 2020-11-08).

## See also

[`substitute_gaps`](https://LRB-IIMCB.github.io/ninetails/reference/substitute_gaps.md)
for gap handling,
[`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)
for the complete feature extraction pipeline

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::filter_signal_by_threshold(
  signal = tail_feature_list[[1]][["readname"]][[4]]
)

} # }
```
