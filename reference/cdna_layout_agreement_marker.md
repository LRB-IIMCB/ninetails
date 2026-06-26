# Agreement marker between the orientation algorithm and the layout inference

Pairs the `tail_type` call from
[`detect_orientation_single`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_single.md)
with the layout call from
[`infer_cdna_layout`](https://LRB-IIMCB.github.io/ninetails/reference/infer_cdna_layout.md)
and returns a single-character marker so disagreements stand out in the
cDNA validation browser.

## Usage

``` r
cdna_layout_agreement_marker(algo_tail_type, inferred_layout)
```

## Arguments

- algo_tail_type:

  Character. The `tail_type` value from the orientation algorithm:
  `"polyA"`, `"polyT"`, or `"unidentified"`.

- inferred_layout:

  Character. The value returned by
  [`infer_cdna_layout`](https://LRB-IIMCB.github.io/ninetails/reference/infer_cdna_layout.md):
  `"polyA_layout"`, `"polyT_layout"`, or `"ambiguous"`.

## Value

Single-character string. One of: check mark (Unicode U+2713) for
agreement, not-equal sign (U+2260) for disagreement, question mark for
layout ambiguity, or em dash (U+2014) for unidentified. Returned as
Unicode escapes in source (`"\u2713"`, `"\u2260"`, `"\u2014"`) for ASCII
portability.

## See also

[`infer_cdna_layout`](https://LRB-IIMCB.github.io/ninetails/reference/infer_cdna_layout.md),
[`launch_cdna_signal_browser`](https://LRB-IIMCB.github.io/ninetails/reference/launch_cdna_signal_browser.md).
