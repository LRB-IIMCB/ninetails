# Plotting ninetails data

This vignette documents all plotting functions provided by **ninetails**
for visualizing classification results, residue composition, poly(A)
tail distributions, and positional patterns. For interactive exploration
of these plots, see
[`vignette("shiny_app")`](https://LRB-IIMCB.github.io/ninetails/articles/shiny_app.md).

All plotting functions return `ggplot2` objects and can be further
customized with standard ggplot2 layers (themes, scales, labels).

## Read classification

[`plot_class_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_class_counts.md)
produces bar charts of read classification results. The `type` argument
controls the level of detail.

``` r
plt <- ninetails::plot_class_counts(
  class_data       = class_data,
  grouping_factor  = "sample_name",
  frequency        = TRUE,
  type             = "N"
)
print(plt)
```

### Parameters

| Parameter         | Type       | Default    | Description                                                |
|-------------------|------------|------------|------------------------------------------------------------|
| `class_data`      | data.frame | *required* | Read classification table from ninetails                   |
| `grouping_factor` | character  | `NA`       | Column name to group by (e.g., `"sample_name"`, `"group"`) |
| `frequency`       | logical    | `FALSE`    | If `TRUE`, show proportions instead of raw counts          |
| `type`            | character  | `"N"`      | Level of detail (see below)                                |

### Classification views

| Type  | Description    | Categories shown                                 |
|-------|----------------|--------------------------------------------------|
| `"N"` | Summary        | decorated, blank, unclassified                   |
| `"R"` | Detailed       | All comment codes (YAY, MAU, MPU, IRL, UNM, BAC) |
| `"A"` | Decorated only | Only decorated reads                             |

![Bar chart showing read classes at the decorated/blank/unclassified
level](../reference/figures/read_classes_2.png)

Bar chart showing read classes at the decorated/blank/unclassified level

### Detailed classification

All reads — decorated, blank, and unclassified — are included and
colour-coded by their comment code.

``` r
plt <- ninetails::plot_class_counts(
  class_data       = class_data,
  grouping_factor  = "sample_name",
  frequency        = TRUE,
  type             = "R"
)
print(plt)
```

![Bar chart showing all read classes broken down by comment
code](../reference/figures/read_classes_1.png)

Bar chart showing all read classes broken down by comment code

------------------------------------------------------------------------

## Residue counts

[`plot_residue_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_residue_counts.md)
shows the distribution of non-adenosine residue types (C, G, U) across
samples or conditions. Two counting modes are available:

- **By residue** (default): Counts every individual non-A residue
  occurrence. A read with two C residues contributes 2 to the C count.
- **By read** (`by_read = TRUE`): Counts each read once per residue
  type. A read with two C residues contributes 1 to the C count.

``` r
plt <- ninetails::plot_residue_counts(
  residue_data     = residue_data,
  grouping_factor  = "sample_name",
  frequency        = TRUE,
  by_read          = FALSE
)
print(plt)
```

### Parameters

| Parameter         | Type       | Default    | Description                                    |
|-------------------|------------|------------|------------------------------------------------|
| `residue_data`    | data.frame | *required* | Non-A residue table from ninetails             |
| `grouping_factor` | character  | `NA`       | Grouping column                                |
| `frequency`       | logical    | `FALSE`    | Show proportions instead of counts             |
| `by_read`         | logical    | `FALSE`    | Count by read instead of by individual residue |

![Bar chart showing counts of C, G, and U residues per
sample](../reference/figures/residues_1.png)

Bar chart showing counts of C, G, and U residues per sample

------------------------------------------------------------------------

## Non-A abundance

[`plot_nonA_abundance()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nonA_abundance.md)
shows the frequency of reads containing one, two, or three or more
separate non-adenosine residues per read. Frequencies are computed
relative to the total number of decorated reads in each sample.

``` r
plt <- ninetails::plot_nonA_abundance(
  residue_data     = residue_data,
  grouping_factor  = "sample_name"
)
print(plt)
```

### Parameters

| Parameter         | Type       | Default    | Description                        |
|-------------------|------------|------------|------------------------------------|
| `residue_data`    | data.frame | *required* | Non-A residue table from ninetails |
| `grouping_factor` | character  | `NA`       | Grouping column                    |

The legend shows three categories: **single** (one non-A per read),
**two** (exactly two), and **more** (three or more).

------------------------------------------------------------------------

## Poly(A) tail length distribution

[`plot_tail_distribution()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_distribution.md)
plots density distributions of poly(A) tail lengths across samples or
conditions. Central tendency lines can be overlaid.

``` r
plt <- ninetails::plot_tail_distribution(
  input_data       = class_data,
  grouping_factor  = "sample_name",
  max_length       = 200,
  value_to_show    = "median",
  ndensity         = TRUE
)
print(plt)
```

### Parameters

| Parameter         | Type       | Default    | Description                                                           |
|-------------------|------------|------------|-----------------------------------------------------------------------|
| `input_data`      | data.frame | *required* | Data frame with `polya_length` column (class data or merged table)    |
| `grouping_factor` | character  | `NA`       | Grouping column                                                       |
| `max_length`      | numeric    | `200`      | Upper limit of the x-axis                                             |
| `value_to_show`   | character  | `NA`       | Central tendency line: `"mean"`, `"median"`, `"mode"`, or `NA` (none) |
| `ndensity`        | logical    | `TRUE`     | If `TRUE`, normalize density so the peak equals 1 across groups       |

Custom color palettes can be applied:

``` r
# Apply a custom palette
plt + ggplot2::scale_color_manual(
  values = c("#D7191C", "#2C7BB6", "#FDAE61", "#ABD9E9")
)
```

------------------------------------------------------------------------

## Non-A position distribution (rug density)

[`plot_rug_density()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_rug_density.md)
visualizes the positional distribution of non-adenosine residues along
poly(A) tails. Each point represents one detected modification, plotted
against its estimated position from the 3’ end (x-axis) and the total
tail length of the read (y-axis). Marginal density curves on both axes
show the overall distribution shape. Call the function once per
nucleotide type.

``` r
plt_C <- ninetails::plot_rug_density(
  residue_data = residue_data,
  base         = "C",
  max_length   = 100
)
print(plt_C)
```

### Parameters

| Parameter      | Type       | Default    | Description                             |
|----------------|------------|------------|-----------------------------------------|
| `residue_data` | data.frame | *required* | Non-A residue table                     |
| `base`         | character  | *required* | Nucleotide type: `"C"`, `"G"`, or `"U"` |
| `max_length`   | numeric    | *required* | Maximum tail length for axis limits     |

> **Note:** This function requires the `cowplot` package. Install with
> `install.packages("cowplot")`.

For large datasets, consider subsampling before plotting to keep the
scatter readable:

``` r
# Subsample to 1000 points per base type
set.seed(42)
rd_sub <- residue_data[residue_data$prediction == "C", ]
if (nrow(rd_sub) > 1000) {
  keep <- sample.int(nrow(rd_sub), 1000, replace = FALSE)
  rd_sub <- rd_sub[keep, ]
}
plt <- ninetails::plot_rug_density(rd_sub, base = "C", max_length = 200)
```

![Scatterplot with rug and density overlay showing positions of C
residues](../reference/figures/rug_plot.png)

Scatterplot with rug and density overlay showing positions of C residues

------------------------------------------------------------------------

## Panel characteristics

[`plot_panel_characteristics()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_panel_characteristics.md)
produces a multi-panel summary of the full tail composition for a
dataset or subset (e.g., a single sample, group, or transcript). The
resulting plot can be aligned from either the 5’ or 3’ end of the tail.

``` r
plt <- ninetails::plot_panel_characteristics(
  input_residue_data            = residue_data,
  input_class_data              = class_data,
  input_merged_nonA_tables_data = NULL,
  type                          = "default",
  max_length                    = 100,
  direction_5_prime             = TRUE
)
print(plt)
```

### Parameters

| Parameter                       | Type       | Default     | Description                                    |
|---------------------------------|------------|-------------|------------------------------------------------|
| `input_residue_data`            | data.frame | *required*  | Non-A residue table                            |
| `input_class_data`              | data.frame | *required*  | Read classification table                      |
| `input_merged_nonA_tables_data` | data.frame | `NULL`      | Merged table (optional)                        |
| `type`                          | character  | `"default"` | Panel layout type                              |
| `max_length`                    | numeric    | `100`       | Maximum tail length                            |
| `direction_5_prime`             | logical    | `TRUE`      | Align from 5’ end (`TRUE`) or 3’ end (`FALSE`) |

The panel contains the following subplots:

- **A** — Read counts
- **B** — Frequency of non-adenosine residues
- **C** — Poly(A) tail length distribution
- **D** — Normalised distribution of non-adenosines (binned by length)
- **E** — Raw distribution of non-adenosines

![Multi-panel plot](../reference/figures/plot_panel.png)

Multi-panel plot

------------------------------------------------------------------------

## Nanopolish QC tags

> **Note:** This section applies only to the Guppy legacy pipeline
> ([`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)).

The
[`nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/nanopolish_qc.md)
and
[`plot_nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nanopolish_qc.md)
functions allow visualization of QC tags assigned by the Nanopolish
polya function.

``` r
# Summarise QC tags per sample
qc_summary <- ninetails::nanopolish_qc(
  class_data,
  grouping_factor = "sample_name"
)

# Plot by read count
plt_count <- ninetails::plot_nanopolish_qc(qc_summary, frequency = FALSE)

# Plot by frequency
plt_freq <- ninetails::plot_nanopolish_qc(qc_summary, frequency = TRUE)

print(plt_count)
print(plt_freq)
```

### Parameters (nanopolish_qc)

| Parameter         | Type       | Default    | Description               |
|-------------------|------------|------------|---------------------------|
| `class_data`      | data.frame | *required* | Read classification table |
| `grouping_factor` | character  | `NA`       | Grouping column           |

### Parameters (plot_nanopolish_qc)

| Parameter         | Type       | Default    | Description                                                                                       |
|-------------------|------------|------------|---------------------------------------------------------------------------------------------------|
| `processing_info` | data.frame | *required* | Output from [`nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/nanopolish_qc.md) |
| `frequency`       | logical    | `FALSE`    | Show proportions instead of counts                                                                |

![Bar charts showing Nanopolish QC tag
distribution](../reference/figures/nanopolish_QC.png)

Bar charts showing Nanopolish QC tag distribution

------------------------------------------------------------------------

## Interactive dashboard

For interactive exploration of all plots described above, ninetails
provides a Shiny dashboard. See
[`vignette("shiny_app")`](https://LRB-IIMCB.github.io/ninetails/articles/shiny_app.md)
for details.

``` r
ninetails::launch_signal_browser(
  class_file   = "/path/to/read_classes.txt",
  residue_file = "/path/to/nonadenosine_residues.txt"
)
```

------------------------------------------------------------------------

## Summary of plotting functions

| Function                                                                                                                                                                    | Description                             | Required data                                                                                |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------|----------------------------------------------------------------------------------------------|
| [`plot_class_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_class_counts.md)                                                                               | Read classification bar chart           | `class_data`                                                                                 |
| [`plot_residue_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_residue_counts.md)                                                                           | Residue type distribution (C/G/U)       | `residue_data`                                                                               |
| [`plot_nonA_abundance()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nonA_abundance.md)                                                                           | Reads with 1, 2, 3+ non-A residues      | `residue_data`                                                                               |
| [`plot_tail_distribution()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_distribution.md)                                                                     | Poly(A) tail length density             | `class_data` or merged table                                                                 |
| [`plot_rug_density()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_rug_density.md)                                                                                 | Non-A position scatterplot with density | `residue_data`                                                                               |
| [`plot_panel_characteristics()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_panel_characteristics.md)                                                             | Multi-panel composition summary         | `class_data` + `residue_data`                                                                |
| [`plot_nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nanopolish_qc.md)                                                                             | Nanopolish QC tag distribution          | [`nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/nanopolish_qc.md) output |
| [`plot_squiggle_fast5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_squiggle_fast5.md)                                                                           | Full read signal (fast5)                | Fast5 files                                                                                  |
| [`plot_squiggle_pod5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_squiggle_pod5.md)                                                                             | Full read signal (POD5)                 | POD5 files                                                                                   |
| [`plot_tail_range_fast5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_range_fast5.md)                                                                       | Poly(A) tail signal (fast5)             | Fast5 files                                                                                  |
| [`plot_tail_range_pod5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_range_pod5.md)                                                                         | Poly(A) tail signal (POD5)              | POD5 files                                                                                   |
| [`plot_tail_chunk()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_chunk.md)                                                                                   | Signal segment                          | Intermediate data                                                                            |
| [`plot_gaf()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_gaf.md) / [`plot_multiple_gaf()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_multiple_gaf.md) | Gramian Angular Fields                  | Intermediate data                                                                            |
