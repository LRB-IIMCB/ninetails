# Correct read class labels for backward compatibility

Resolves naming convention changes between pre-release and current
versions of ninetails. Previously, read classes were labeled
`"modified"`, `"unmodified"`, and `"unclassified"`. These have been
changed to `"decorated"`, `"blank"`, and `"unclassified"`, since the
presence of non-A within a poly(A) tail does not necessarily result from
post-transcriptional modification; it may be caused by semi-templated
tails, polymerase slippage, or other mechanisms.

## Usage

``` r
correct_labels(df)
```

## Arguments

- df:

  Data frame. Ninetails read classification results (`class_data`) or
  the output of
  [`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md).
  Must contain a `class` column.

## Value

Data frame with corrected class labels in the `class` column.

## Examples

``` r
if (FALSE) { # \dontrun{

class_data <- ninetails::correct_labels(class_data)

} # }
```
