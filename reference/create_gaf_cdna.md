# Create gramian angular fields from cDNA signal windows

Converts normalised signal windows into two-channel images: the gramian
angular summation field `cos(phi_i + phi_j)` with `phi = acos(x)` and
the field `sin(psi_i + psi_j)` with `psi = asin(x)`.

## Usage

``` r
create_gaf_cdna(windows)
```

## Arguments

- windows:

  Numeric matrix. One normalised window per row (as returned by
  [`create_tail_windows_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_windows_cdna.md)).

## Value

A numeric array with dimensions (number of windows, window length,
window length, 2).

## Details

In contrast to
[`create_gaf`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf.md),
the windows are **not** rescaled to their own minimum and maximum: the
normalisation of
[`normalise_tail_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/normalise_tail_cdna.md)
already maps the signal to (-1, 1) and keeps the amplitude relative to
the whole tail, which carries most of the information in cDNA signals.

## See also

[`create_tail_windows_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_windows_cdna.md),
[`predict_nonA_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/predict_nonA_cdna.md)

## Examples

``` r
if (FALSE) { # \dontrun{

gafs <- ninetails::create_gaf_cdna(windows = tail_windows$windows)

} # }
```
