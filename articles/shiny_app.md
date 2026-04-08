# Analysis dashboard (Shiny)

## Overview

Ninetails includes an interactive Shiny dashboard for exploring analysis
results. The dashboard provides six tabs for classification overview,
residue composition, poly(A) length distribution, raw signal
visualization, configurable report download, and package information.

The dashboard supports two usage modes:

- **Single-sample mode** — provide file paths directly; all tabs are
  available when classification and residue data are supplied.
- **Multi-sample mode** — provide a YAML configuration file pointing to
  multiple samples; enables comparative analysis across samples and
  experimental groups.

## Requirements

The dashboard requires several optional packages. Install them before
first use:

``` r
install.packages(c("shiny", "plotly", "htmltools", "DT", "base64enc"))

# For multi-sample mode (YAML config):
install.packages("yaml")

# For rug density plots (position distribution):
install.packages("cowplot")
```

The Signal Viewer tab additionally requires a working Python environment
with the `pod5` package, accessible via `reticulate`. See the
[installation guide](https://github.com/LRB-IIMCB/ninetails/wiki) for
details.

## Quick start

### Single-sample mode

The simplest way to launch the dashboard with results from a single
ninetails run:

``` r
ninetails::launch_signal_browser(
  summary_file = "/path/to/dorado_summary.txt",
  pod5_dir     = "/path/to/pod5/",
  class_file   = "/path/to/read_classes.txt",
  residue_file = "/path/to/nonadenosine_residues.txt"
)
```

All four arguments are optional. When `class_file` and `residue_file`
are provided, the Classification, Residues, and Poly(A) length tabs
become active. When `summary_file` and `pod5_dir` are provided, the
Signal Viewer tab enables raw signal visualization.

If only signal data paths are provided (no `class_file`), only the
Signal Viewer tab is functional:

``` r
ninetails::launch_signal_browser(
  summary_file = "/path/to/dorado_summary.txt",
  pod5_dir     = "/path/to/pod5/"
)
```

### Multi-sample mode

For comparative analysis across multiple samples and experimental
conditions, prepare a YAML configuration file and launch with:

``` r
ninetails::launch_signal_browser(
  config = "/path/to/config.yml"
)
```

## YAML configuration file

The YAML configuration file defines samples, experimental groups, and
file paths. Each sample entry must include `class_path` and
`residue_path` (paths to ninetails output files). The `dorado_summary`
and `pod5_dir` fields are optional and only needed for the Signal Viewer
tab.

### Format

``` yaml
samples:
  KO_rep1:
    sample_name: KO_rep1
    group: KO
    class_path: /path/to/KO_rep1/read_classes.txt
    residue_path: /path/to/KO_rep1/nonadenosine_residues.txt
    dorado_summary: /path/to/KO_rep1/dorado_summary.txt   # optional
    pod5_dir: /path/to/KO_rep1/pod5/                      # optional

  KO_rep2:
    sample_name: KO_rep2
    group: KO
    class_path: /path/to/KO_rep2/read_classes.txt
    residue_path: /path/to/KO_rep2/nonadenosine_residues.txt

  WT_rep1:
    sample_name: WT_rep1
    group: WT
    class_path: /path/to/WT_rep1/read_classes.txt
    residue_path: /path/to/WT_rep1/nonadenosine_residues.txt
    dorado_summary: /path/to/WT_rep1/dorado_summary.txt
    pod5_dir: /path/to/WT_rep1/pod5/
```

### Field reference

| Field            | Required | Description                                        |
|------------------|----------|----------------------------------------------------|
| `sample_name`    | Yes      | Display name used in plots and legends             |
| `group`          | Yes      | Experimental condition for grouping                |
| `class_path`     | Yes      | Path to `read_classes.txt` from ninetails          |
| `residue_path`   | Yes      | Path to `nonadenosine_residues.txt` from ninetails |
| `dorado_summary` | No       | Path to Dorado summary file (for Signal Viewer)    |
| `pod5_dir`       | No       | Path to POD5 file directory (for Signal Viewer)    |

A config template is included in the package:

``` r
system.file("extdata", "config_template.yml", package = "ninetails")
```

## Dashboard tabs

### Classification

The first tab provides an overview of the analysis results.

**Value boxes** at the top display five metrics in a row of equal-width
cards, each with a colored border and matching label: number of samples
analyzed, number of transcripts found, total read count, number of blank
reads, and number of decorated reads.

**Read Classification plot** shows the distribution of read categories
(decorated, blank, unclassified). Three views are available via the
“Classification view” dropdown:

- **Summary (N)** — simplified view with main categories
- **Detailed (R)** — expanded view with all comment codes (YAY, MPU,
  MAU, IRL, BAC, UNM)
- **Decorated (A)** — decorated reads only

**Non-A Abundance plot** shows the frequency of reads containing one,
two, or three or more separate non-adenosine residues per read
([`plot_nonA_abundance()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nonA_abundance.md)).

When the browser window is wide enough (≥1200px), both plots are
displayed side-by-side in a responsive flex layout. On narrower screens
they stack vertically.

**Sidebar controls:**

- *Grouping variable* — choose between `sample_name` or `group`
  (detected dynamically from the data)
- *Filter by transcript* — server-side searchable dropdown for
  transcript-level filtering (uses `symbol` if annotated, otherwise
  `contig`)
- *Show as frequency* — toggle between raw counts and proportions
- *Show descriptions* — toggle explanatory text below each plot,
  including classification code definitions (YAY, MPU, MAU, IRL, BAC,
  UNM)

### Residues

The Residues tab provides three views of non-adenosine residue
composition.

**Residue Counts plot** shows the distribution of C, G, and U residues
([`plot_residue_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_residue_counts.md)).
Toggle between counting by read (each read counted once per residue
type) or by residue (total occurrences).

**Non-A Position Distribution** displays rug density plots showing where
non-adenosine residues are located along the poly(A) tail relative to
tail length
([`plot_rug_density()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_rug_density.md)).
One row of three plots (C, G, U) is generated per sample. To keep the
display readable and rendering fast, each residue type is subsampled to
a maximum of 1,000 points per sample (without replacement). If fewer
than 1,000 points exist for a given residue type in a sample, all
available points are plotted. If a residue type is absent from a sample,
a “No X residues” placeholder is shown.

**Summary Table** (when merged data is available) provides
per-transcript statistics including read counts, non-A counts and hits,
and poly(A) length metrics
([`summarize_nonA()`](https://LRB-IIMCB.github.io/ninetails/reference/summarize_nonA.md)).
The table is searchable and sortable via DataTables.

**Sidebar controls:**

- *Grouping variable* — `sample_name` or `group`
- *Filter by transcript* — transcript-level filter
- *Show as frequency* — counts vs proportions
- *Count by read* — toggle between by-read and by-residue counting
- *Max tail length (rug plot)* — controls rug density plot x-axis
- *Show descriptions* — toggle explanatory text below plots

### Poly(A) length

The Poly(A) length tab shows poly(A) tail length distributions across
samples or groups
([`plot_tail_distribution()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_distribution.md)).

**Sidebar controls:**

- *Grouping variable* — `sample_name` or `group`
- *Filter by transcript* — transcript-level filter
- *Select condition(s)* — multi-select dropdown to plot only specific
  conditions (samples or groups). All conditions are selected by
  default. Remove items to compare subsets; changes take effect
  immediately.
- *Central tendency* — overlay mean, median, mode, or none
- *Max tail length* — upper limit for the x-axis (slider, 0–500)
- *Normalized density* — toggle between raw density and normalized
  density
- *Color palette* — choose from 8 built-in palettes: kanto, johto,
  hoenn, sinnoh, hisui, unova, kalos, alola
- *Reset all filters* — restores every control to its default value
- *Show descriptions* — toggle explanatory text below the plot

### Signal Viewer

The Signal Viewer tab provides interactive visualization of raw nanopore
signals from POD5 files. It requires `dorado_summary` and `pod5_dir`
paths (either from the launcher or from the YAML config).

In multi-sample mode, a dropdown at the top selects which sample’s
signals to browse. In single-sample mode, paths are loaded automatically
from the launcher arguments.

Two sub-tabs are available:

**Static Viewer** — renders two ggplot2 plots:

1.  *Entire Signal* — full-length squiggle with colored segments
    (adapter in blue, poly(A) in orange, transcript in black)
2.  *Poly(A) Region* — zoomed view centered on the poly(A) tail with
    ±250 positions flanking

Both plots show vertical dashed lines marking the 5’ (red) and 3’ (navy)
boundaries of the poly(A) tail. When non-A residue data is available,
semi-transparent rectangles highlight the estimated positions of
detected modifications, colored by residue type (C, G, U) with a letter
label at the top.

**Dynamic Explorer** — renders an interactive Plotly chart of the full
signal with zoom, pan, and hover capabilities. Non-A highlights are
shown as translucent rectangles. Double-click to reset the view.

**Signal Viewer sidebar controls:**

- *Minimum poly(A) length* — filter reads by tail length
- *Decoration status* — filter by comment code (YAY, MPU, MAU, etc.)
- *Non-A residue type* — show only reads with a specific residue (C, G,
  U)
- *Alignment genome* — filter by reference (if available in summary)
- *Mapping quality range* — MAPQ slider filter
- *Select Read ID* — server-side searchable dropdown with Previous/Next
  navigation buttons

### Download

The Download tab provides a configurable report generator. The report is
a self-contained HTML file with all plots embedded as base64-encoded PNG
images, suitable for sharing or archiving without external dependencies.

**Report sections** (each controlled by a checkbox):

- *Classification plot* — global read classification summary
- *Non-A abundance* — frequency of reads with 1, 2, 3+ non-A residues
- *Residue frequency* — distribution of C, G, U residues
- *Rug density plots* — per-sample positional distribution of non-A
  residues with automatic downsampling to 1,000 points and annotations
  indicating whether downsampling was applied
- *Poly(A) length distribution* — density plot with configurable central
  tendency
- *Example signal plots* — for each sample with POD5 access, 5 randomly
  selected reads per category (blank, decorated-C, decorated-G,
  decorated-U) are rendered as full-width tail range plots with non-A
  overlay. If a category has fewer than 5 reads, all available reads are
  shown. Categories absent from a sample display a placeholder message.

**Plot settings** (apply to all poly(A) plots in the report):

- *Central tendency* — mean, median, mode, or none (annotated in the
  report)
- *Max poly(A) length* — upper x-axis limit
- *Color palette* — palette applied to poly(A) density curves

**Per-transcript sections** (optional, up to 3 transcripts):

Select up to three transcripts for detailed sub-reports. Each transcript
section includes its own classification, residue frequency, and poly(A)
length distribution plots, filtered to reads mapping to that transcript.

All plots in the report include descriptive annotations explaining what
is shown, consistent with the “Show descriptions” text available in the
interactive tabs.

### About

The About tab displays the ninetails logo, package version, full
citation (Gumińska et al., *Nat Commun* 2025), links to the GitHub
repository, Wiki documentation, pkgdown website, and Zenodo DOI, as well
as the IIMCB logo, laboratory information (Laboratory of RNA Biology,
ERA Chairs Group), and developer contact details.

## Annotated data

When the input data has been annotated with
[`annotate_with_biomart()`](https://LRB-IIMCB.github.io/ninetails/reference/annotate_with_biomart.md)
(adding `symbol` and Ensembl transcript ID columns), the dashboard
automatically uses `symbol` as the transcript label in all “Filter by
transcript” dropdowns. Without annotation, the raw `contig` column is
used instead.

This applies to both single-sample and multi-sample modes. Annotation
should be performed before launching the dashboard:

``` r
# After running the ninetails pipeline:
class_data_annotated <- ninetails::annotate_with_biomart(
  class_data,
  species = "mmusculus"
)

residue_data_annotated <- ninetails::annotate_with_biomart(
  residue_data,
  species = "mmusculus"
)

# Save annotated data, then point the dashboard to the annotated files
```

## Launcher reference

``` r
?ninetails::launch_signal_browser
```

| Argument       | Type      | Description                                                                                   |
|----------------|-----------|-----------------------------------------------------------------------------------------------|
| `config`       | character | Path to YAML config (multi-sample mode)                                                       |
| `summary_file` | character | Path to Dorado summary file                                                                   |
| `pod5_dir`     | character | Path to POD5 directory                                                                        |
| `class_file`   | character | Path to `read_classes.txt`                                                                    |
| `residue_file` | character | Path to `nonadenosine_residues.txt`                                                           |
| `...`          |           | Additional arguments passed to [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) |

When `config` is provided, single-sample arguments are ignored. All
arguments are optional; the dashboard adapts its active tabs based on
which data is available.

Additional
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) arguments
can be passed through `...`:

``` r
# Custom port and host
ninetails::launch_signal_browser(
  config = "config.yml",
  port = 8080,
  host = "0.0.0.0"
)

# Suppress browser launch
ninetails::launch_signal_browser(
  config = "config.yml",
  launch.browser = FALSE
)
```

## Static assets

The dashboard expects three image files in `inst/app/www/`:

- `logo.png` — ninetails logo (copy from `man/figures/logo.png`)
- `favicon.ico` — browser tab icon (copy from `pkgdown/favicon/`)
- `IIMCB_logo.png` — IIMCB institute logo (displayed in the About tab)

## Deploying to Shiny Server

The dashboard can be deployed to Shiny Server (Open Source or Pro) for
shared access over a network. A deployment wrapper is included in the
package that handles data loading from a YAML configuration file and
delegates rendering to the installed ninetails app.

### Setup

``` r
# 1. Copy the deployment wrapper to the Shiny Server app directory
deploy_dir <- "/srv/shiny-server/ninetails"
dir.create(deploy_dir, recursive = TRUE, showWarnings = FALSE)

file.copy(
  system.file("deployment", "app.R", package = "ninetails"),
  file.path(deploy_dir, "app.R"),
  overwrite = TRUE
)

# 2. Copy your YAML configuration
file.copy("config.yml", file.path(deploy_dir, "config.yml"))

# 3. Symlink static assets from the installed package
file.symlink(
  system.file("app", "www", package = "ninetails"),
  file.path(deploy_dir, "www")
)
```

### Configuring the deployment wrapper

Before launching, edit the top of the deployed `app.R` to set two
variables:

- **`config_path`** — path to the YAML configuration file (default:
  `config.yml` in the same directory as `app.R`)
- **`PYTHON_PATH`** — path to the Python binary that has the `pod5`
  module installed. Required for the Signal Viewer tab. Set to `NULL` if
  signal visualization is not needed.

``` r
# In /srv/shiny-server/ninetails/app.R:
config_path <- file.path(getwd(), "config.yml")
PYTHON_PATH <- "/home/user/miniconda3/envs/r-reticulate/bin/python"
```

Common Python path examples:

| Environment   | Typical path                                         |
|---------------|------------------------------------------------------|
| System Python | `/usr/bin/python3`                                   |
| Conda env     | `/home/user/miniconda3/envs/r-reticulate/bin/python` |
| virtualenv    | `/home/user/.virtualenvs/ninetails/bin/python`       |
| pyenv         | `/home/user/.pyenv/versions/3.11.0/bin/python`       |

You can find the correct path by running `which python` in the
environment where `pod5` is installed.

### Directory structure

    /srv/shiny-server/ninetails/
    ├── app.R          # deployment wrapper (from inst/deployment/)
    ├── config.yml     # your YAML configuration
    └── www/           # symlink to inst/app/www/ in the installed package
        ├── logo.png
        ├── favicon.ico
        └── IIMCB_logo.png

### Requirements

- ninetails and all Suggests dependencies must be installed system-wide
  (accessible to the `shiny` user)
- All data file paths in `config.yml` must be readable by the `shiny`
  user
- For the Signal Viewer tab:
  - `PYTHON_PATH` must be set in the deployment wrapper
  - The specified Python environment must have the `pod5` module
    installed (`pip install pod5`)
  - The `shiny` user must have execute permissions on the Python binary
    and read access to the POD5 files

After placing the files, restart Shiny Server:

``` bash
sudo systemctl restart shiny-server
```

The dashboard will be accessible at
`http://your-server:3838/ninetails/`.

> **Note:** When you update the ninetails package, the deployed app
> picks up changes automatically (it sources from the installed package
> at runtime). To change the dataset, edit `config.yml` and restart
> Shiny Server. To change the Python environment, edit `PYTHON_PATH` in
> `app.R` and restart.

## Troubleshooting

**Rug density plots are empty or show an error**
[`plot_rug_density()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_rug_density.md)
requires the `cowplot` package. Install it with
`install.packages("cowplot")`.

**Download button produces an error** The report generation requires
`base64enc` for embedding plots as base64 images. Install with
`install.packages("base64enc")`.

**Signal Viewer shows “No POD5 file found”** Ensure the `pod5_dir` path
points to the directory containing `.pod5` files (not a parent
directory). The app searches recursively within that directory.

**Signal Viewer shows “Extraction failed”** The Python `pod5` package
must be installed and accessible via `reticulate`. Check with:

``` r
reticulate::py_module_available("pod5")
```

**Classification or Residue plots are blank** If using single-sample
mode, make sure you provided `class_file` and/or `residue_file`. These
are optional arguments — without them, the corresponding tabs show
placeholder messages.

**“Column `more` doesn’t exist” in Non-A Abundance plot** This occurs
when no reads have 3 or more non-A residues. Update
[`plot_nonA_abundance()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nonA_abundance.md)
to the latest version which handles missing columns after pivot.

**Y-axis shows scientific notation** The dashboard disables scientific
notation globally with `options(scipen = 999)` at startup. If you still
see scientific notation, ensure you are running the latest version of
the app.
