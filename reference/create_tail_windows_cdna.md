# Split a normalized cDNA tail signal into sliding windows

Cuts 100-sample windows every 25 samples, with window centers running
from 50 samples before the tail to 50 samples behind it. Windows that
would leave the provided signal are skipped (no padding).

## Usage

``` r
create_tail_windows_cdna(signal, tail_start, tail_end)
```

## Arguments

- signal:

  Numeric vector. Raw signal (tail with flanks).

- tail_start:

  Numeric. First sample of the tail within `signal` (0-based).

- tail_end:

  Numeric. End of the tail within `signal` (exclusive).

## Value

A list with two elements:

- centres:

  Numeric vector. Window centres within `signal` (0-based)

- windows:

  Numeric matrix with one normalised window per row and 100 columns

## See also

[`normalise_tail_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/normalise_tail_cdna.md),
[`create_gaf_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_cdna.md)

## Examples

``` r
if (FALSE) { # \dontrun{

tail_windows <- ninetails::create_tail_windows_cdna(signal = signal,
                                                    tail_start = 100,
                                                    tail_end = 800)

} # }
```
