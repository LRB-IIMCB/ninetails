# Scatterplot of nonA residue positions within poly(A) tail

This function allows to produce the scatterplot of raw non-A residue
predictions (y-axis) along the user-defined tail span (x-axis). The
resulting plot also contains the rugs and ridges along the axes, which
are a visual aid designed to better show the distributions of given
residue.

## Usage

``` r
plot_rug_density(residue_data, base, max_length)
```

## Arguments

- residue_data:

  A dataframe or tibble containig non-A residue predictions made by
  ninetails pipeline.

- base:

  character. One of the following \["C"/"G","U"\]. This parameter
  defines which nonadenosine is to be plotted. It is obligatory to
  prevent the overplotting.

- max_length:

  numeric \[1\]. This parameter controls maximum length of the poly(A)
  tail to be taken into consideration for plotting.

## Value

a ggplot2 object

## Examples

``` r
if (FALSE) { # \dontrun{
ninetails::plot_rug_density(residue_data=residue_data, base="C", max_length=100)
} # }
```
