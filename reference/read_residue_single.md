# Reads ninetails nonadenosine_residues data from file.

Imports a single `nonadenosine_residues` output file produced by
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
(or other ninetails pipelines) into R. Optionally attaches a sample
identifier and parses GENCODE-style contig names into Ensembl transcript
IDs.

## Usage

``` r
read_residue_single(residue_path, sample_name = NA)
```

## Arguments

- residue_path:

  Character string. Path to the ninetails `nonadenosine_residues` output
  file (`.tsv`).

- sample_name:

  Character string (optional, default `NA`). If specified, added as a
  factor column `sample_name`. Overwrites any existing `sample_name`
  column with a warning.

## Value

A tibble containing non-A residue predictions with the following
additional columns (when GENCODE contigs are present):

- transcript:

  Character. Gene symbol parsed from GENCODE contig.

- ensembl_transcript_id_full:

  Character. Ensembl transcript ID with version number.

- ensembl_transcript_id_short:

  Character. Ensembl transcript ID without version number.

## Details

If the `contig` column contains GENCODE-formatted identifiers, three
additional columns are created: `transcript`,
`ensembl_transcript_id_full`, and `ensembl_transcript_id_short` (see
[`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
for the same logic applied to class data).

## Acknowledgements

Function based on `read_polya_single` from the NanoTail package by P.
Krawczyk (smaegol): <https://github.com/LRB-IIMCB/nanotail/>.

## See also

[`read_residue_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_multiple.md)
for batch loading,
[`read_class_single`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
for the analogous class data loader,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for the pipeline that produces the input file.

## Examples

``` r
if (FALSE) { # \dontrun{

residue_path <- "/directory/with/ninetails/nonadenosine_residues_output.tsv"
residue_data <- ninetails::read_residue_single(residue_path)

} # }
```
