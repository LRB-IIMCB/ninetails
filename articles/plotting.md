# Plotting ninetails data

## Read classification

**Ninetails** provides insights into the non-adenosine landscape of
sequenced samples. The
[`plot_class_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_class_counts.md)
function produces graphical representations of read classification
results. The `type` argument controls the level of detail:

- `type = "R"` — detailed view: all comment codes shown (YAY, MAU, MPU,
  IRL, UNM, BAC)
- `type = "N"` — coarser view: three classes only (decorated, blank,
  unclassified)
- `type = "A"` — decorated reads only

### Detailed classification

All reads — decorated, blank, and those that could not be classified —
are included and colour-coded by their comment code.

``` r
plt1 <- ninetails::plot_class_counts(class_data = class_data,
                                     grouping_factor = "sample_name",
                                     frequency = TRUE,
                                     type = "R")
print(plt1)
```

![Bar chart showing all read classes broken down by comment code (YAY,
MAU, MPU, IRL, UNM/BAC)](../reference/figures/read_classes_1.png)

Bar chart showing all read classes broken down by comment code (YAY,
MAU, MPU, IRL, UNM/BAC)

### Less detailed classification

With `type = "N"`, only the three top-level classes are shown, making it
easier to compare the proportion of decorated reads across samples.

``` r
plt1 <- ninetails::plot_class_counts(class_data = class_data,
                                     grouping_factor = "sample_name",
                                     frequency = TRUE,
                                     type = "N")
print(plt1)
```

![Bar chart showing read classes at the decorated/blank/unclassified
level](../reference/figures/read_classes_2.png)

Bar chart showing read classes at the decorated/blank/unclassified level

## Residue classification

The
[`plot_residue_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_residue_counts.md)
function plots the frequency or counts of each predicted non-adenosine
residue type (C, G, U).

``` r
plt1 <- ninetails::plot_residue_counts(residue_data = residue_data,
                                       grouping_factor = "sample_name")
print(plt1)
```

![Bar chart showing counts or frequencies of C, G, and U residues per
sample](../reference/figures/residues_1.png)

Bar chart showing counts or frequencies of C, G, and U residues per
sample

## Plot panel characteristics

The
[`plot_panel_characteristics()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_panel_characteristics.md)
function produces a multi-panel summary of the full tail composition for
a dataset or subset (e.g. a single sample, group, or transcript). The
resulting plot can be aligned from either the 5’ or 3’ end of the tail.

``` r
plt <- ninetails::plot_panel_characteristics(
  input_residue_data = residue_data,
  input_class_data = class_data,
  input_merged_nonA_tables_data = NULL,
  type = "default",
  max_length = 100,
  direction_5_prime = TRUE
)
print(plt)
```

The panel contains the following subplots:

- **A** — read counts
- **B** — frequency of non-adenosine residues
- **C** — poly(A) tail length distribution
- **D** — normalised distribution of non-adenosines (binned and
  normalised by length)
- **E** — raw distribution of non-adenosines

![Multi-panel plot showing read counts, residue frequencies, tail length
distribution, and positional distributions of
non-adenosines](../reference/figures/plot_panel.png)

Multi-panel plot showing read counts, residue frequencies, tail length
distribution, and positional distributions of non-adenosines

## Position scatterplot

The
[`plot_rug_density()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_rug_density.md)
function visualises non-adenosine positions relative to tail lengths,
combining a scatterplot, density plot, and rug in one figure. Call the
function once per nucleotide type.

``` r
plt <- ninetails::plot_rug_density(residue_data = residue_data,
                                   base = "C",
                                   max_length = 100)
print(plt)
```

![Scatterplot with rug and density overlay showing positions of C
residues along poly(A) tails](../reference/figures/rug_plot.png)

Scatterplot with rug and density overlay showing positions of C residues
along poly(A) tails

## Non-A abundance

The
[`plot_nonA_abundance()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nonA_abundance.md)
function visualises the overall abundance of non-adenosine residues
across samples.

``` r
plt <- ninetails::plot_nonA_abundance(residue_data = residue_data,
                                      grouping_factor = "sample_name")
print(plt)
```

## Tail length distribution

The
[`plot_tail_distribution()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_distribution.md)
function plots the distribution of poly(A) tail lengths across samples.

``` r
plt <- ninetails::plot_tail_distribution(input_data = class_data,
                                         grouping_factor = "sample_name")
print(plt)
```

## Nanopolish QC tags

> **Note**
>
> This section applies only to the Guppy legacy pipeline
> ([`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)).

The
[`nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/nanopolish_qc.md)
and
[`plot_nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nanopolish_qc.md)
functions allow visualisation of QC tags assigned by the Nanopolish
polya function, either by read count or by frequency.

``` r
# Summarise QC tags per sample
qc_summary <- ninetails::nanopolish_qc(class_data,
                                       grouping_factor = "sample_name")

# Plot by read count
plt1 <- ninetails::plot_nanopolish_qc(qc_summary, frequency = FALSE)

# Plot by frequency
plt2 <- ninetails::plot_nanopolish_qc(qc_summary, frequency = TRUE)

print(plt1)
print(plt2)
```

![Bar charts showing Nanopolish QC tag distribution by count and by
frequency](../reference/figures/nanopolish_QC.png)

Bar charts showing Nanopolish QC tag distribution by count and by
frequency
