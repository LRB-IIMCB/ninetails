# Convert ONT signal to Gramian Angular Field

Represents time series data (ONT squiggle) in a polar coordinate system
instead of the typical Cartesian coordinates. This is a pure-R
equivalent of the Python pyts.image implementation.

## Usage

``` r
create_gaf(tail_chunk, method = "s")
```

## Arguments

- tail_chunk:

  Numeric vector. A 100-element signal chunk representing a fragment of
  the poly(A) tail within the analyzed dataset.

- method:

  Character string specifying the type of Gramian Angular Field. `"s"`
  (default) produces a summation field (GASF), `"d"` produces a
  difference field (GADF).

## Value

An array of dimensions (100, 100, 1) with values representing the
Gramian Angular Field transformation of the input signal chunk. Values
are in the range \[0, 1\].

## Details

Two methods of transformation are available: Gramian Angular Summation
Field (GASF) and Gramian Angular Difference Field (GADF).

The transformation proceeds as follows:

1.  Rescale signal values to the interval \[-1, 1\]

2.  Compute the inverse trigonometric function (`acos` for GASF, `asin`
    for GADF) of each value

3.  Form a pairwise matrix by summing (GASF) or differencing (GADF) the
    angle vectors

4.  Apply the corresponding trigonometric function (`cos` or `sin`) to
    the matrix

5.  Reshape to a 100x100x1 array

6.  Rescale result to the interval \[0, 1\]

## See also

[`combine_gafs`](https://LRB-IIMCB.github.io/ninetails/reference/combine_gafs.md)
for combining GASF and GADF into a single array,
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
for batch GAF computation

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::create_gaf(
  tail_chunk = tail_chunk,
  method = "s"
)

} # }
```
