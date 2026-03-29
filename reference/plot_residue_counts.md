# Plot counts of nonadenosine residues found in ninetails output data

This function may report the non-A occurrences in two flavours: by read
or by residue. To illustrate this, let's have a look at the following
example:

## Usage

``` r
plot_residue_counts(
  residue_data,
  grouping_factor = NA,
  by_read = FALSE,
  frequency = TRUE
)
```

## Arguments

- residue_data:

  A dataframe or tibble containig non-A residue predictions made by
  ninetails pipeline

- grouping_factor:

  grouping variable (e.g. "sample_name")

- by_read:

  logical \[TRUE/FALSE\]. If TRUE, the count/frequency data per reads
  with given residues will be plotted (i.e. how many reads contain given
  residue. If FALSE, the number of residues will be plotted). Set to
  FALSE by default.

- frequency:

  logical \[TRUE/FALSE\]. If TRUE, the frequency will be plotted. If
  FALSE, raw counts will be shown. This parameter is set to TRUE by
  default.

## Value

a ggplot2 object

## Details

read_1: contains 3xU residues

read_2: contains 2xC residue, 1xG residue, 1xU residue

function launched in by read mode (by_read==TRUE): read_1 - taken into
account once, because it non-As of the same (single) type read_2 - taken
into account thrice, because it has non-As of three types

function launched in by residue mode (by_read==FALSE): read_1 - taken
into account thrice, because it has 3 non-A residues reported read_2 -
taken into account fourfold, because it has 4 nonAs reported

The user can switch between those modes depending on the desired
informations to be reported.

## Examples

``` r
if (FALSE) { # \dontrun{
ninetails::plot_residue_counts(residue_data=nonadenosine_residues_dataframe,
                               grouping_factor="sample_name",
                               frequency=TRUE)
} # }
```
