# Calculate the statistical mode of a numeric vector

Computes the statistical mode (most frequent value) using one of two
methods. The `"density"` method (default) returns the mode from the
normalized kernel density estimate and always produces a single value.
The `"value"` method returns the actual most frequent value(s) and may
return multiple values for multimodal data.

## Usage

``` r
get_mode(x, method = "density", na.rm = FALSE)
```

## Arguments

- x:

  Numeric vector. Values to compute the mode for.

- method:

  Character string. Which method to use for computing the statistical
  mode. Two options: `"density"` (default) for the kernel density
  estimate mode, or `"value"` for the actual most frequent value(s).

- na.rm:

  Logical. Whether to remove `NA` values before computation. Default:
  `FALSE`.

## Value

Numeric. The statistical mode of the input vector. Returns a single
value when `method = "density"`, or one or more values when
`method = "value"` (for multimodal data).

## Acknowledgements

Written based on the following Stack Overflow thread:
<https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode>.
Special thanks to Chris and hugovdberg.

## Examples

``` r
if (FALSE) { # \dontrun{

test1 <- c(rep(2, 5), rep(3, 4), rep(1, 4), rep(8, 2))
test2 <- c(rep(2, 5), rep(3, 4), rep(1, 4), rep(8, 5))

result <- ninetails::get_mode(x = test1, method = "value")
print(result)

} # }
```
