# Normalise a raw cDNA tail signal

Centres the signal on the tail median, scales it by the sample-to-sample
noise of the tail and squashes it to the (-1, 1) range:
`tanh(((signal - median(tail)) / noise) / 8)`, where
`noise = 1.4826 * mad(diff(tail)) / sqrt(2)`.

## Usage

``` r
normalise_tail_cdna(signal, tail_start, tail_end)
```

## Arguments

- signal:

  Numeric vector. Raw signal (tail with flanks).

- tail_start:

  Numeric. First sample of the tail within `signal` (0-based).

- tail_end:

  Numeric. End of the tail within `signal` (exclusive).

## Value

Numeric vector of the same length as `signal`.

## Details

The statistics come from the tail only, but the whole provided signal
(tail with flanks) is transformed. The result does not depend on the
signal units (raw ADC values or picoamperes give the same output).

## See also

[`create_tail_windows_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_windows_cdna.md)

## Examples

``` r
if (FALSE) { # \dontrun{

normalised <- ninetails::normalise_tail_cdna(signal = signal,
                                             tail_start = 100,
                                             tail_end = 800)

} # }
```
