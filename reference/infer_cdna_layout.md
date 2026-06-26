# Infer cDNA signal layout from pre- and post-tail region sizes

Independent visual cross-check of the orientation algorithm
([`detect_orientation_single`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_single.md)).
Compares the sizes of the pre-tail and post-tail regions of the raw
signal; the larger side is assumed to be the transcript body. Per the
empirical orientation finding, cDNA polyA reads carry the transcript
body before the tail (body-first), while cDNA polyT reads carry it after
(adapter-first).

## Usage

``` r
infer_cdna_layout(
  poly_tail_start,
  poly_tail_end,
  signal_length,
  min_region_samples = 100L,
  strictness_ratio = 1.5
)
```

## Arguments

- poly_tail_start:

  Integer. 1-based start index of the poly tail in signal coordinates.

- poly_tail_end:

  Integer. 1-based end index of the poly tail in signal coordinates.

- signal_length:

  Integer. Total length of the raw signal vector.

- min_region_samples:

  Integer. Minimum acceptable pre- or post-tail region size in samples.
  Defaults to 100.

- strictness_ratio:

  Numeric. Required ratio of larger to smaller region for a
  non-ambiguous call. Defaults to 1.5.

## Value

Character string: `"polyA_layout"` (body \| tail \| adapter),
`"polyT_layout"` (adapter \| tail \| body), or `"ambiguous"`.

## Details

Returns `"ambiguous"` when either region is shorter than a minimum size
guard, or when the larger / smaller ratio is below a strictness
threshold. Both guards are deliberate - the function is intended as a
strict secondary check, so borderline reads are flagged for manual
review rather than silently committed to one orientation.

Used by
[`launch_cdna_signal_browser`](https://LRB-IIMCB.github.io/ninetails/reference/launch_cdna_signal_browser.md)
to colour signal regions and surface algorithm-vs-layout disagreements.

## See also

[`cdna_layout_agreement_marker`](https://LRB-IIMCB.github.io/ninetails/reference/cdna_layout_agreement_marker.md),
[`launch_cdna_signal_browser`](https://LRB-IIMCB.github.io/ninetails/reference/launch_cdna_signal_browser.md).
