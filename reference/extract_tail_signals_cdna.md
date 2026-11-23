# Extract raw poly(A)/poly(T) tail signals with flanks from POD5 files (cDNA)

Extracts the raw signal of each tail together with 100 samples on each
side (clipped at the read ends). In contrast to
[`extract_tails_from_pod5`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tails_from_pod5.md),
the signal is neither winsorized nor downsampled, because the cDNA
models work on the raw 5 kHz signal.

## Usage

``` r
extract_tail_signals_cdna(tail_coordinates, pod5_dir, num_cores = 1)
```

## Arguments

- tail_coordinates:

  Data frame. Output of
  [`read_tail_coordinates_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/read_tail_coordinates_cdna.md).

- pod5_dir:

  Character string. Path to the directory with POD5 files.

- num_cores:

  Numeric \[1\]. Number of physical cores to use.

## Value

A named list of numeric vectors (raw signal), one per read with valid
coordinates, named with read IDs. The tail starts at sample
`min(tail_start, 100)` (0-based) of each vector.

## Details

If the `filename` column is not available (Lachesis output), every POD5
file in `pod5_dir` is scanned for the requested reads.

## See also

[`read_tail_coordinates_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/read_tail_coordinates_cdna.md),
[`predict_nonA_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/predict_nonA_cdna.md)

## Examples

``` r
if (FALSE) { # \dontrun{

signal_list <- ninetails::extract_tail_signals_cdna(
  tail_coordinates = coordinates,
  pod5_dir = '/path/to/pod5/',
  num_cores = 2)

} # }
```
