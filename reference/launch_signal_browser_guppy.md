# Launch the Ninetails Analysis Dashboard (Guppy Legacy)

Opens an interactive Shiny application for exploring ninetails results
from the Guppy legacy pipeline (`check_tails_guppy`). This is the
fast5-based counterpart of
[`launch_signal_browser`](https://LRB-IIMCB.github.io/ninetails/reference/launch_signal_browser.md),
designed for data basecalled with Guppy \\\le\\ 6.0.0 using Nanopolish
poly(A) coordinates.

## Usage

``` r
launch_signal_browser_guppy(
  config = NULL,
  nanopolish_file = NULL,
  sequencing_summary_file = NULL,
  workspace = NULL,
  class_file = NULL,
  residue_file = NULL,
  basecall_group = "Basecall_1D_000",
  ...
)
```

## Arguments

- config:

  Character string (optional). Path to a YAML configuration file
  defining multiple samples. When provided, single-sample arguments are
  ignored.

- nanopolish_file:

  Character string (optional). Path to the Nanopolish polya output file.
  Required for the Signal Viewer tab in single-sample mode.

- sequencing_summary_file:

  Character string (optional). Path to the Guppy sequencing summary
  file. Required for the Signal Viewer tab in single-sample mode.

- workspace:

  Character string (optional). Path to the directory containing
  multi-fast5 files. Required for the Signal Viewer tab in single-sample
  mode.

- class_file:

  Character string (optional). Path to the `read_classes` output file
  from
  [`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md).

- residue_file:

  Character string (optional). Path to a `nonadenosine_residues` output
  file from ninetails.

- basecall_group:

  Character string. Fast5 hierarchy level for basecall data extraction.
  Default: `"Basecall_1D_000"`.

- ...:

  Additional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Value

Launches a Shiny application (does not return a value).

## Details

The dashboard provides the same six tabs as the Dorado version
(Classification, Residues, Poly(A) length, Signal Viewer, Download,
About), with two Guppy-specific additions:

- **Nanopolish QC**:

  Additional plot in the Classification tab showing the distribution of
  Nanopolish QC tags (PASS, ADAPTER, NOREGION, SUFFCLIP, etc.) via
  [`plot_nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nanopolish_qc.md).

- **Fast5 Signal Viewer**:

  Uses
  [`plot_squiggle_fast5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_squiggle_fast5.md)
  and
  [`plot_tail_range_fast5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_range_fast5.md)
  instead of POD5-based functions. Includes a `moves` toggle to
  show/hide basecaller move transitions. No Python dependency required.

The YAML configuration file should have the following structure:


    samples:
      sample_label:
        sample_name: WT_rep1
        group: WT
        class_path: /path/to/read_classes.txt
        residue_path: /path/to/nonadenosine_residues.txt
        polya_path: /path/to/nanopolish_output.tsv         # optional
        seq_summary: /path/to/sequencing_summary.txt      # optional
        workspace: /path/to/fast5/                        # optional

## See also

[`launch_signal_browser`](https://LRB-IIMCB.github.io/ninetails/reference/launch_signal_browser.md)
for the Dorado DRS version,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for the Guppy pipeline.

## Examples

``` r
if (FALSE) { # \dontrun{

# Multi-sample mode
ninetails::launch_signal_browser_guppy(config = "config_guppy.yml")

# Single-sample: all tabs
ninetails::launch_signal_browser_guppy(
  nanopolish_file = "/path/to/nanopolish_output.tsv",
  sequencing_summary_file  = "/path/to/sequencing_summary.txt",
  workspace = "/path/to/fast5/",
  class_file = "/path/to/read_classes.txt",
  residue_file = "/path/to/nonadenosine_residues.txt"
)

} # }
```
