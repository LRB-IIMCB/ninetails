# Detect poly tail type for a single sequence using Dorado-style algorithm

This function implements the Dorado-style poly tail detection algorithm
for a single cDNA sequence. It uses edit distance matching of SSP and
VNP primer sequences at both ends of the read to determine orientation,
testing both forward (polyA) and reverse (polyT) strand configurations.

## Usage

``` r
detect_orientation_single(sequence)
```

## Arguments

- sequence:

  Character string containing the DNA sequence to classify

## Value

Character string: "A" for polyA orientation, "T" for polyT orientation,
or "unknown" for unclassified sequences

## Examples

``` r
if (FALSE) { # \dontrun{
detect_orientation_single("TTTCTGTTGGTGCTGATATTGCTTT...")  # Returns "A"
detect_orientation_single("ACTTGCCTGTCGCTCTATCTTCAG...")   # Returns "T"
} # }
```
