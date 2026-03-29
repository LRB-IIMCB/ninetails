# Count trailing occurrences of a character in a string

Counts how many times a specific character appears consecutively at the
end of a string, working backwards from the last character. Used
internally by the cDNA pipeline for poly(A)/poly(T) tail length
estimation from basecalled sequences.

## Usage

``` r
count_trailing_chars(string, char)
```

## Arguments

- string:

  Character string. The sequence to analyze.

- char:

  Character string of length 1. The character to count at the end of the
  string.

## Value

Integer. Count of trailing character occurrences.

## See also

[`detect_orientation_single`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_single.md)
where trailing character counts are used for read orientation
classification

## Examples

``` r
if (FALSE) { # \dontrun{

count_trailing_chars("ACGTTTTT", "T")
# Returns 5

count_trailing_chars("ACGT", "A")
# Returns 0

} # }
```
