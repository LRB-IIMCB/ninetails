# Plotting read classes data per category assigned to the analyzed reads.

This function requires as input the dataframe with read_classes provided
by ninetails pipeline. Function works in 3 flavours, by plotting either:

- detailed classification (based on column "comments")

- crude classification (based on column "class")

- reads decorated with non-As exclusively

## Usage

``` r
plot_class_counts(
  class_data,
  grouping_factor = NA,
  frequency = TRUE,
  type = "R"
)
```

## Arguments

- class_data:

  A dataframe or tibble containing read_classes predictions made by
  ninetails pipeline

- grouping_factor:

  character string. A grouping variable (e.g. "sample_name")

- frequency:

  logical \[TRUE/FALSE\]. If TRUE, the frequency will be plotted. If
  FALSE, raw counts will be shown. This parameter is set to TRUE by
  default.

- type:

  character string \["R"/"N"/"A"\]. This variable controls the level of
  detail of the resulting plot:

  - "R" - detailed classification (based on column "comments")

  - "N" - crude classification (based on column "class")

  - "A" - reads decorated with non-As exclusively

  By default, the "R" option is set.

## Value

ggplot object with read class prediction

## Details

Function based on the Nanotail equivalent:
https://github.com/LRB-IIMCB/nanotail Many thanks to smeagol (Pawel
Krawczyk) for advice & support!

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::plot_class_counts(class_data=read_classes_dataframe,
                             grouping_factor="sample_name",
                             frequency=TRUE,
                             type="R")
} # }

```
