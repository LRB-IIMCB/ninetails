# Predict non-adenosine residues in cDNA tails

Classifies sliding windows of raw tail signals with the orientation-
specific cDNA model and calls non-adenosine residues.

## Usage

``` r
predict_nonA_cdna(signal_list, tail_coordinates, tail_type, keras_model_path)
```

## Arguments

- signal_list:

  Named list of numeric vectors. Raw tail signals with flanks, as
  returned by
  [`extract_tail_signals_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_cdna.md).

- tail_coordinates:

  Data frame. Output of
  [`read_tail_coordinates_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/read_tail_coordinates_cdna.md)
  (columns `read_id`, `tail_start`, `tail_end` and `polya_length` are
  used).

- tail_type:

  Character string. Orientation of the reads: `"polyA"` or `"polyT"`.

- keras_model_path:

  Character string or missing. Full path of the `.h5` file with a Keras
  model. If missing (default), the pretrained cDNA model for the given
  `tail_type` is loaded.

## Value

A data frame with one row per called non-adenosine residue:

- read_id:

  Character. Read identifier

- prediction:

  Character. Residue in the RNA tail (C, G or U)

- est_nonA_pos:

  Numeric. Estimated position from the 3' end of the tail (nucleotides)

- score:

  Numeric. Highest non-A score of the call

## Details

For every read the tail signal is normalised
([`normalise_tail_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/normalise_tail_cdna.md)),
cut into windows
([`create_tail_windows_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_windows_cdna.md)),
converted to gramian angular fields
([`create_gaf_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_cdna.md))
and classified. Calls are made with
[`call_nonA_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/call_nonA_cdna.md)
using thresholds frozen together with the models.

The position of each call is reported as the distance from the 3' end of
the RNA tail in nucleotides. Poly(A) reads end with the 3' end of the
tail, poly(T) reads begin with it:

- poly(A): `(tail_end - position) / samples_per_nt`

- poly(T): `(position - tail_start) / samples_per_nt`

where `samples_per_nt = (tail_end - tail_start) / polya_length`.

## See also

[`create_outputs_dorado_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado_cdna.md),
[`load_keras_model`](https://LRB-IIMCB.github.io/ninetails/reference/load_keras_model.md)

## Examples

``` r
if (FALSE) { # \dontrun{

nonA_calls <- ninetails::predict_nonA_cdna(
  signal_list = signal_list,
  tail_coordinates = coordinates,
  tail_type = "polyA")

} # }
```
