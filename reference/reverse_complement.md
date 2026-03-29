# Generate reverse complement of a DNA sequence

Creates the reverse complement of a DNA sequence string, handling
standard nucleotide bases (A, T, G, C) and ambiguous bases (N). Used
internally by the cDNA pipeline for primer matching on both strands.

## Usage

``` r
reverse_complement(sequence)
```

## Arguments

- sequence:

  Character string. A DNA sequence containing only the characters A, T,
  G, C, and N.

## Value

Character string. The reverse complement of the input sequence.

## See also

[`detect_orientation_single`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_single.md)
where reverse complement sequences are used for primer alignment

## Examples

``` r
if (FALSE) { # \dontrun{

reverse_complement("ATCG")
# Returns "CGAT"

reverse_complement("AAATTT")
# Returns "AAATTT"

} # }
```
