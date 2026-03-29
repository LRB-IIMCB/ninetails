# Reads ninetails read_classes data frame from file.

Imports a single `read_classes` output file produced by
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
(or other ninetails pipelines) into R. Optionally attaches a sample
identifier and parses GENCODE-style contig names into Ensembl transcript
IDs.

## Usage

``` r
read_class_single(class_path, sample_name = NA)
```

## Arguments

- class_path:

  Character string. Path to the ninetails `read_classes` output file
  (`.tsv`).

- sample_name:

  Character string (optional, default `NA`). If specified, added as a
  factor column `sample_name`. Overwrites any existing `sample_name`
  column with a warning.

## Value

A tibble containing read_classes predictions with the following
additional columns (when GENCODE contigs are present):

- transcript:

  Character. Gene symbol parsed from GENCODE contig.

- ensembl_transcript_id_full:

  Character. Ensembl transcript ID with version number.

- ensembl_transcript_id_short:

  Character. Ensembl transcript ID without version number.

## Details

After loading, the function applies
[`correct_labels`](https://LRB-IIMCB.github.io/ninetails/reference/correct_labels.md)
to ensure backward compatibility with the label changes introduced in
v0.9 (decorated/blank/unclassified). If the `contig` column contains
GENCODE-formatted identifiers, three additional columns are created:
`transcript` (gene symbol), `ensembl_transcript_id_full` (with version),
and `ensembl_transcript_id_short` (without version).

## Acknowledgements

Function based on `read_polya_single` from the NanoTail package by P.
Krawczyk (smaegol): <https://github.com/LRB-IIMCB/nanotail/>.

## See also

[`read_class_multiple`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_multiple.md)
for batch loading,
[`correct_labels`](https://LRB-IIMCB.github.io/ninetails/reference/correct_labels.md)
for the backward-compatibility step,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for the pipeline that produces the input file.

## Examples

``` r
if (FALSE) { # \dontrun{

class_path <- "/directory/with/ninetails/read_class_output.tsv"
class_data <- ninetails::read_class_single(class_path)

} # }
```
