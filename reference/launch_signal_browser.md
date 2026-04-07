# Launch the Ninetails Analysis Dashboard

Opens an interactive Shiny application for exploring and visualizing
ninetails poly(A) tail composition analysis results. The dashboard
provides a comprehensive overview of read classification, non-adenosine
residue composition, poly(A) tail length distributions, and raw nanopore
signal inspection — all in a single interface with interactive filters,
per-sample breakdowns, and a configurable report generator.

## Usage

``` r
launch_signal_browser(
  config = NULL,
  summary_file = NULL,
  pod5_dir = NULL,
  class_file = NULL,
  residue_file = NULL,
  ...
)
```

## Arguments

- config:

  Character string (optional). Path to a YAML configuration file
  defining multiple samples with file paths and experimental groups.
  When provided, single-sample arguments are ignored. See Details for
  the expected format.

- summary_file:

  Character string (optional). Path to a Dorado summary file
  (tab-separated, with columns `read_id`, `filename`, `poly_tail_start`,
  `poly_tail_end`). Required for the Signal Viewer tab in single-sample
  mode.

- pod5_dir:

  Character string (optional). Path to the directory containing POD5
  files from the sequencing run. Required for the Signal Viewer tab in
  single-sample mode.

- class_file:

  Character string (optional). Path to the `read_classes` output file
  from a ninetails pipeline (`check_tails_dorado_DRS`,
  `check_tails_dorado_cDNA`, or `check_tails_guppy`). Enables the
  Classification, Residues, Poly(A) length, and Download tabs.

- residue_file:

  Character string (optional). Path to a `nonadenosine_residues` output
  file from ninetails. Enables residue-specific plots (abundance,
  residue counts, rug density) and non-A overlay on signal plots.

- ...:

  Additional arguments passed to
  [`runApp`](https://rdrr.io/pkg/shiny/man/runApp.html), such as `port`,
  `host`, or `launch.browser`.

## Value

Launches a Shiny application (does not return a value).

## Details

The dashboard is organized into six tabs:

- **Classification**:

  Summary value boxes (samples, transcripts, total/blank/decorated
  reads) and bar charts showing read classification and non-A abundance
  per sample or condition. Three classification views are available:
  summary, detailed (by comment code), and decorated-only.

- **Residues**:

  Distribution of non-adenosine residue types (C, G, U) with per-sample
  rug density plots showing positional distribution along poly(A) tails
  (subsampled to 1,000 points per residue type per sample), and a
  searchable summary table when merged data is available.

- **Poly(A) length**:

  Density distributions of poly(A) tail lengths with condition
  filtering, central tendency overlays, 8 color palettes, and a summary
  statistics table (n, mean, median, SD, SEM) that updates with all
  filters.

- **Signal Viewer**:

  Raw nanopore signal visualization from POD5 files with two sub-tabs:
  Static Viewer (ggplot2 full signal and zoomed poly(A) region) and
  Dynamic Explorer (interactive Plotly with zoom/pan). Includes non-A
  residue overlay highlighting, filterable read lists, and Previous/Next
  navigation.

- **Download**:

  Configurable report generator producing a self-contained HTML file.
  Users select which sections to include (classification, abundance,
  residue frequency, rug density, poly(A) distribution, example signal
  plots) and configure plot settings. Supports per-transcript
  sub-reports for up to 3 selected transcripts. All plots include
  descriptive annotations.

- **About**:

  Package version, citation (Nat Commun 2025), links to GitHub, Wiki,
  pkgdown site, Zenodo DOI, laboratory and developer contact
  information.

The dashboard supports two usage modes:

**Multi-sample mode** (`config`): Loads a YAML configuration file
pointing to multiple samples with class and residue data. Enables
comparative analysis across samples and experimental groups. All
plotting functions use `sample_name` and `group` columns for grouping.

**Single-sample mode**: Provide individual file paths. All analysis tabs
are available when `class_file` and `residue_file` are supplied. The
Signal Viewer tab requires `summary_file` and `pod5_dir`. Default
`sample_name` and `group` columns are added automatically if missing.

The YAML configuration file should have the following structure:

    samples:
      sample_label_1:
        sample_name: KO_rep1
        group: KO
        class_path: /path/to/read_classes.txt
        residue_path: /path/to/nonadenosine_residues.txt
        dorado_summary: /path/to/dorado_summary.txt   # optional
        pod5_dir: /path/to/pod5/                      # optional

Required fields per sample: `sample_name`, `group`, `class_path`,
`residue_path`. The `dorado_summary` and `pod5_dir` fields are optional
and only needed for the Signal Viewer tab.

A template configuration file is included in the package:
`system.file("extdata", "config_template.yml", package = "ninetails")`

The dashboard requires the following packages (listed in Suggests):
`shiny`, `plotly`, `htmltools`, `DT`, `base64enc`, `cowplot`. For
multi-sample mode: `yaml`. For the Signal Viewer: a Python environment
with the `pod5` package, accessible via `reticulate`.

Static assets (`logo.png`, `favicon.ico`, `IIMCB_logo.png`) should be
placed in `inst/app/www/`.

## See also

[`check_tails_dorado_DRS`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_DRS.md)
for the main analysis pipeline,
[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for merging output tables,
[`annotate_with_biomart`](https://LRB-IIMCB.github.io/ninetails/reference/annotate_with_biomart.md)
for adding gene symbols,
[`plot_class_counts`](https://LRB-IIMCB.github.io/ninetails/reference/plot_class_counts.md)
and other plotting functions for static equivalents of dashboard plots.

## Examples

``` r
if (FALSE) { # \dontrun{

# Multi-sample mode (recommended for comparative analysis)
ninetails::launch_signal_browser(config = "/path/to/config.yml")

# Single-sample: all tabs (classification + residues + signal viewer)
ninetails::launch_signal_browser(
  summary_file = "/path/to/dorado_summary.txt",
  pod5_dir     = "/path/to/pod5/",
  class_file   = "/path/to/read_classes.txt",
  residue_file = "/path/to/nonadenosine_residues.txt"
)

# Single-sample: signal viewer only (no analysis tabs)
ninetails::launch_signal_browser(
  summary_file = "/path/to/dorado_summary.txt",
  pod5_dir     = "/path/to/pod5/"
)

# Custom port and host for remote access
ninetails::launch_signal_browser(
  config = "config.yml",
  port = 8080,
  host = "0.0.0.0"
)

} # }
```
