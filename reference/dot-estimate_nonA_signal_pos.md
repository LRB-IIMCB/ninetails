# Convert estimated non-A position to raw signal coordinate

Reverses the position estimation formula used in
[`create_outputs_dorado`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado.md)
to map `est_nonA_pos` (nucleotide distance from the 3' end of the
poly(A) tail) back to a position in the raw signal vector.

## Usage

``` r
.estimate_nonA_signal_pos(
  est_nonA_pos,
  poly_tail_length,
  poly_tail_start,
  poly_tail_end
)
```

## Arguments

- est_nonA_pos:

  Numeric vector. Estimated non-A position(s) from 3' end (as reported
  in `nonadenosine_residues` table).

- poly_tail_length:

  Numeric. Total poly(A) tail length (nt).

- poly_tail_start:

  Numeric. Poly(A) start coordinate in raw signal.

- poly_tail_end:

  Numeric. Poly(A) end coordinate in raw signal.

## Value

Numeric vector of raw signal positions (rounded to integer).

## Details

The original estimation:
`est_nonA_pos = poly_tail_length - (poly_tail_length * centr_signal_pos / signal_length)`

where `signal_length = 0.2 * (poly_tail_end - poly_tail_start)`.

The 0.2x interpolation and 5x rescale cancel algebraically, yielding:
`raw_pos = poly_tail_start + (poly_tail_length - est_nonA_pos) * (poly_tail_end - poly_tail_start) / poly_tail_length`
