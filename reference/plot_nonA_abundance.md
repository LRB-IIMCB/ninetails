# Plot abundances of reads with given amount of non-A residues per read

This function plots frequencies of reads containing one, two or more
separate instances (occurrences) of non-As reported by ninetails. The
frequency is computed with respect to the total amount of decorated
reads in the analyzed dataset.

## Usage

``` r
plot_nonA_abundance(residue_data, grouping_factor = NA)
```

## Arguments

- residue_data:

  A dataframe or tibble containig non-A residue predictions made by
  ninetails pipeline

- grouping_factor:

  character string. A grouping variable (e.g. "sample_name")

## Value

ggplot2 object

## Examples

``` r
if (FALSE) { # \dontrun{
plt <- ninetails::plot_nonA_abundance(residue_data = residue_data,
                                      grouping_factor = "sample_name")
plt
} # }
```
