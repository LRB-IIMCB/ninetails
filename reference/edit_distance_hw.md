# Calculate edit distance with sliding window (HW mode)

Computes the minimum edit distance between a query sequence and a target
sequence using a sliding window approach (HW: Hamming-like Window mode).
The query is aligned to all possible substrings of the target of equal
length, and the minimum distance across all alignments is returned.

## Usage

``` r
edit_distance_hw(query, target)
```

## Arguments

- query:

  Character string. The primer or query sequence to search for.

- target:

  Character string. The target sequence window to search within.

## Value

Integer. The minimum Damerau–Levenshtein edit distance found across all
possible sliding window alignments.

## Details

Edit distance is computed using the Damerau–Levenshtein distance, as
implemented by [`adist`](https://rdrr.io/r/utils/adist.html), which
accounts for insertions, deletions, substitutions, and transpositions.

The sliding-window strategy for edit distance calculation is adapted
from approaches commonly used in natural language processing; see the
reference below for an illustrative example.

## References

Silge, J. (2019). \*Natural Language Processing in R: Edit Distance\*.
R-bloggers.
<https://www.r-bloggers.com/2019/04/natural-language-processing-in-r-edit-distance/>

## See also

[`detect_orientation_single`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_single.md)
where this function is used for fuzzy primer matching in read
orientation classification.

## Examples

``` r
if (FALSE) { # \dontrun{

edit_distance_hw("ATCG", "XXATCGXX")
# Returns 0 (perfect match inside target)

edit_distance_hw("ATCG", "ATTT")
# Returns minimum distance across sliding windows

} # }
```
