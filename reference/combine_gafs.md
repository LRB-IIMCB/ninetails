# Combine GASF and GADF into a two-channel array

Creates a two-dimensional array containing both the Gramian Angular
Summation Field (GASF) and the Gramian Angular Difference Field (GADF)
produced from a single ONT tail chunk. Using both representations
increases the sensitivity of the CNN classification by overcoming the
limitations of each method individually.

## Usage

``` r
combine_gafs(tail_chunk)
```

## Arguments

- tail_chunk:

  Numeric vector. A 100-element signal chunk representing a fragment of
  the poly(A) tail.

## Value

An array of dimensions (100, 100, 2) where the first channel contains
the GASF and the second channel contains the GADF, both produced by
[`create_gaf`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf.md).

## See also

[`create_gaf`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf.md)
for individual GAF computation,
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
for batch processing,
[`predict_gaf_classes`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md)
for downstream CNN classification

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::combine_gafs(tail_chunk = tail_chunk)

} # }
```
