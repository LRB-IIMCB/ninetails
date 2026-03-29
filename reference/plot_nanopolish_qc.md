# Plots qc data (qc_tag) inherited from nanopolish polya function.

Provides insight into the proportion of each quality tag category per
sample or experiment condition defined by the user.

## Usage

``` r
plot_nanopolish_qc(processing_info, frequency = TRUE)
```

## Arguments

- processing_info:

  the output of nanopolish_qc function

- frequency:

  logical \[TRUE/FALSE\].

## Value

a ggplot object

## Details

This is the ninetails' implementation of the `name` function originally
written by P. Krawczyk (smeagol) and incorporated within the NanoTail
package.

For original source code, see:
https://github.com/LRB-IIMCB/nanotail/blob/master/R/polya_plots.R

The variable names were adjusted according to the naming convention
within ninetails to avoid confusion.

Many thanks to the author of original source code for help and advice.
