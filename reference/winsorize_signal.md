# Winsorize nanopore signal

Clips extreme values in the raw nanopore signal to the 0.5th and 99.5th
percentiles, removing high cliffs that extend above or below the main
signal range. This slightly affects signal extremes but preserves the
overall signal shape for downstream analysis.

## Usage

``` r
winsorize_signal(signal)
```

## Arguments

- signal:

  Numeric vector. A raw ONT signal vector.

## Value

Integer vector. The winsorized signal with extreme values clipped to the
0.5% and 99.5% quantiles.

## Acknowledgements

Based on A. Signorell's DescTools:
<https://github.com/AndriSignorell/DescTools/blob/master/R/DescTools.r>.

## See also

[`extract_tail_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data.md)
and
[`extract_tails_from_pod5`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tails_from_pod5.md)
where winsorization is applied during tail signal extraction

## Examples

``` r
if (FALSE) { # \dontrun{

winsorized_signal <- winsorize_signal(signal = sample(200:300))

} # }
```
