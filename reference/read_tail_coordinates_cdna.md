# Read poly(A)/poly(T) tail coordinates for the cDNA pipeline

Converts the tail coordinates produced by Lachesis or Dorado into one
common table used by the cDNA pipeline, so that the downstream functions
do not depend on the tool that delimited the tail.

## Usage

``` r
read_tail_coordinates_cdna(tail_coordinates, source = "lachesis")
```

## Arguments

- tail_coordinates:

  Character string or data frame. Path to the Lachesis output / Dorado
  summary file, or the same table in memory.

- source:

  Character string \["lachesis"\]. Tool that produced the coordinates:
  `"lachesis"` or `"dorado"`.

## Value

A data frame with one row per read:

- read_id:

  Character. Read identifier

- filename:

  Character. POD5 file name (`NA` if not available)

- contig:

  Character. Reference name (`NA` if not available)

- tail_start:

  Numeric. First sample of the tail (0-based)

- tail_end:

  Numeric. End of the tail (exclusive)

- polya_length:

  Numeric. Estimated tail length in nucleotides

- qc_tag:

  Character. Lachesis QC tag or Dorado mapping quality

- tail_type:

  Character. `"polyA"`, `"polyT"`, `"unidentified"` or `NA`

## Details

Coordinates are positions in the raw signal (samples), with the tail
spanning `signal[(tail_start + 1):tail_end]`.

## Lachesis input

The tail is delimited by `adapter_end` and `polya_end` in both
orientations. The orientation is taken from the `strand` column: `"+"`
is a poly(A) read, `"-"` is a poly(T) read.

## Dorado input

The tail is delimited by `poly_tail_start` and `poly_tail_end`. The
summary does not carry the read orientation, so `tail_type` is `NA` and
has to be filled in by the caller (see
[`detect_orientation_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_multiple.md)).
Both `alignment_mapq` and `alignment_mapping_quality`, as well as
`filename` and `input_filename`, are accepted.

## See also

[`extract_tail_signals_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_cdna.md),
[`check_tails_dorado_cDNA`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_cDNA.md)

## Examples

``` r
if (FALSE) { # \dontrun{

coordinates <- ninetails::read_tail_coordinates_cdna(
  tail_coordinates = '/path/to/lachesis_polya.tsv',
  source = "lachesis")

} # }
```
