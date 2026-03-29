# Detection of outliers (peaks & valleys) in ONT signal using z-scores.

This is the training-set variant of the z-score-based peak/valley
detection function. It identifies deviations in the poly(A) tail signal
that may correspond to non-adenosine nucleotides (C, G, or U) by
applying a moving-average z-score filter with empirically calibrated
parameters.

## Usage

``` r
filter_signal_by_threshold_trainingset(signal)
```

## Arguments

- signal:

  Numeric vector. An ONT read fragment corresponding to the poly(A) tail
  region as delimited by nanopolish polya (stored in
  `tail_feature_list[[1]]` produced by
  [`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md)).

## Value

Numeric vector of pseudomoves with values in {-1, 0, 1}:

- 1:

  Signal exceeds baseline (peak) — potential G nucleotide.

- -1:

  Signal falls below baseline (valley) — potential C or U nucleotide.

- 0:

  Signal within baseline range — homopolymer A region.

## Details

The algorithm prepends 100 calibration data points (sampled from the 10
most frequent signal values and the first 10 observations) to stabilise
the baseline estimate. It then applies a sliding window of 100 data
points: any observation exceeding 3.5 standard deviations from the
running mean is classified as a peak (+1) or valley (-1); all other
positions are scored 0. After trimming the calibration prefix, terminal
positions are zeroed out (first 5 and last 5) to prevent edge artefacts
from entering downstream chunk extraction. Finally,
[`substitute_gaps`](https://LRB-IIMCB.github.io/ninetails/reference/substitute_gaps.md)
is applied to fill isolated zero-gaps in the pseudomove vector.

## Acknowledgements

The z-score peak detection algorithm is based on: Brakel, J.P.G. van
(2014). *"Robust peak detection algorithm using z-scores"*. Stack
Overflow. Available at:
<https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362>
(version: 2020-11-08).

## See also

[`extract_tail_data_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md)
which calls this function,
[`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)
for the production counterpart,
[`substitute_gaps`](https://LRB-IIMCB.github.io/ninetails/reference/substitute_gaps.md)
for the gap-filling post-processing step.

## Examples

``` r
if (FALSE) { # \dontrun{

filter_signal_by_threshold_trainingset(
  signal = tail_feature_list[[1]][["readname"]][[4]])

} # }
```
