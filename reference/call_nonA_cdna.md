# Call non-adenosine residues from window-wise cDNA predictions

Turns the class probabilities of consecutive signal windows into
non-adenosine calls. Windows whose non-A score (`1 - p(A)`) exceeds the
threshold are merged into runs; each sufficiently long run gives one
call.

## Usage

``` r
call_nonA_cdna(
  centres,
  probabilities,
  tail_start,
  tail_end,
  theta_mid,
  theta_edge,
  min_run = 2
)
```

## Arguments

- centres:

  Numeric vector. Window centres (0-based, same coordinate system as
  `tail_start` and `tail_end`).

- probabilities:

  Numeric matrix. Class probabilities with one row per window and four
  columns (A, C, G, U).

- tail_start:

  Numeric. First sample of the tail (0-based).

- tail_end:

  Numeric. End of the tail (exclusive).

- theta_mid:

  Numeric. Threshold for windows inside the tail.

- theta_edge:

  Numeric. Threshold for windows near a tail boundary.

- min_run:

  Numeric \[2\]. Minimal number of consecutive windows above the
  threshold required for a call.

## Value

A data frame with one row per call:

- prediction:

  Character. Residue in the RNA tail (C, G or U)

- centr_signal_pos:

  Numeric. Position of the call in the signal

- score:

  Numeric. Highest non-A score of the run

## Details

Windows centred within 100 samples of a tail boundary are compared with
`theta_edge`, the remaining ones with `theta_mid`. The call is
positioned at the centre of the highest-scoring window of the run and
its class is the one with the highest summed probability over the run.
Calls positioned outside the tail are dropped.

## See also

[`predict_nonA_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/predict_nonA_cdna.md)

## Examples

``` r
if (FALSE) { # \dontrun{

calls <- ninetails::call_nonA_cdna(centres = tail_windows$centres,
                                   probabilities = probabilities,
                                   tail_start = 100,
                                   tail_end = 800,
                                   theta_mid = 0.95,
                                   theta_edge = 0.98)

} # }
```
