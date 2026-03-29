# Reshapes nonadenosine_residues data frame to wide format.

Pivots the `prediction` column of the `nonadenosine_residues` data frame
into three separate columns (`prediction_C`, `prediction_G`,
`prediction_U`), each containing the count of respective non-A residues
per read. An additional `nonA_residues` column provides a CIGAR-like
summary string listing all non-A positions from 5' to 3', separated by
`":"`.

## Usage

``` r
spread_nonA_residues(residue_data)
```

## Arguments

- residue_data:

  Data frame or tibble containing non-A residue predictions produced by
  the ninetails pipeline (filename typically ends with
  `_nonadenosine_residues.tsv`). Required columns: `readname`, `group`,
  `prediction`, `est_nonA_pos`. The `qc_tag` column may be character
  (Guppy/nanopolish) or numeric (Dorado MAPQ).

## Value

A tibble in wide format with one row per read. In addition to the
original columns (minus `est_nonA_pos` and `prediction`), the following
are added:

- prediction_C:

  Integer. Number of C residues detected in the read.

- prediction_G:

  Integer. Number of G residues detected in the read.

- prediction_U:

  Integer. Number of U residues detected in the read.

- nonA_residues:

  Character. CIGAR-like string summarising all non-A positions (e.g.
  `"C12:G25:U40"`).

## Details

This function supports data from both the legacy Guppy/nanopolish
pipeline and the Dorado pipeline. The `qc_tag` column may be either
character (Guppy: `"PASS"`, `"SUFFCLIP"`, etc.) or numeric (Dorado: MAPQ
score). Both types are handled transparently and preserved in the
output.

## See also

[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
which calls this function internally,
[`read_residue_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
for loading residue data.

## Examples

``` r
if (FALSE) { # \dontrun{

spread_table <- ninetails::spread_nonA_residues(
  residue_data = residue_data)

} # }
```
