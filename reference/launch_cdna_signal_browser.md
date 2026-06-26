# Launch the cDNA orientation validation browser

Minimum-useful companion to
[`launch_signal_browser`](https://LRB-IIMCB.github.io/ninetails/reference/launch_signal_browser.md)
specifically for verifying cDNA orientation classifications against the
raw signal layout. Shows one signal at a time with the algorithm's
`tail_type` call, an independent size-ratio-inferred layout (see
[`infer_cdna_layout`](https://LRB-IIMCB.github.io/ninetails/reference/infer_cdna_layout.md)),
and an agreement marker so disagreements stand out visually.

## Usage

``` r
launch_cdna_signal_browser(dorado_summary, pod5_dir, orientation_data, ...)
```

## Arguments

- dorado_summary:

  Character string or data frame. Path to a Dorado summary file or a
  data frame with at least `read_id`, `filename`, `poly_tail_start`,
  `poly_tail_end` columns. The Dorado \>=1.4.0 `input_filename` alias is
  also accepted.

- pod5_dir:

  Character string. Path to the directory containing POD5 files for the
  run.

- orientation_data:

  Character string or data frame. Path to a TSV produced by
  [`detect_orientation_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_multiple.md),
  or that function's in-memory return value. Must contain `read_id` and
  `tail_type`; `reference` and `sequence` are used if present.

- ...:

  Additional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

Launches a Shiny application (does not return a value).

## Details

Intended as a validation tool while
[`check_tails_dorado_cDNA`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_cDNA.md)
is still under construction; will be superseded by a fuller cDNA
dashboard once cDNA-native CNN submodels are in place.

## See also

[`launch_signal_browser`](https://LRB-IIMCB.github.io/ninetails/reference/launch_signal_browser.md),
[`detect_orientation_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_multiple.md),
[`infer_cdna_layout`](https://LRB-IIMCB.github.io/ninetails/reference/infer_cdna_layout.md),
[`cdna_layout_agreement_marker`](https://LRB-IIMCB.github.io/ninetails/reference/cdna_layout_agreement_marker.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Using files on disk
ninetails::launch_cdna_signal_browser(
  dorado_summary   = "/path/to/dorado_summary.txt",
  pod5_dir         = "/path/to/pod5/",
  orientation_data = "/path/to/sequence_with_tail_type.tsv"
)

# Using an in-memory orientation tibble (e.g. straight out of the pipeline)
classified <- ninetails::detect_orientation_multiple(
  sequence_files = my_seq_files, num_cores = 4
)
ninetails::launch_cdna_signal_browser(
  dorado_summary   = "/path/to/dorado_summary.txt",
  pod5_dir         = "/path/to/pod5/",
  orientation_data = classified
)
} # }
```
