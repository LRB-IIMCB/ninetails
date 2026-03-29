# Annotate ninetails output data with biomaRt

Retrieves gene-level annotation from Ensembl for transcripts in
ninetails output data. This is a convenience wrapper for
[`getBM`](https://huber-group-embl.github.io/biomaRt/reference/getBM.html)
with built-in organism presets and support for custom mart objects.

## Usage

``` r
annotate_with_biomart(
  input_data,
  attributes_to_get = c("ensembl_transcript_id", "external_gene_name", "description",
    "transcript_biotype"),
  filters = "ensembl_transcript_id",
  organism = NULL,
  mart_to_use = NULL
)
```

## Arguments

- input_data:

  Data frame. Tabular output of the ninetails pipeline. Must contain a
  column named `ensembl_transcript_id_short`.

- attributes_to_get:

  Character vector. Annotation attributes to retrieve from biomaRt.
  Default:
  `c("ensembl_transcript_id", "external_gene_name", "description", "transcript_biotype")`.

- filters:

  Character string. Column of the input data frame to match with the
  target mart. Default: `"ensembl_transcript_id"`.

- organism:

  Character string or `NULL`. Organism shorthand for built-in mart
  presets. Currently available:

  `"athaliana"`

  :   *Arabidopsis thaliana* (Ensembl Plants)

  `"hsapiens"`

  :   *Homo sapiens* (Ensembl)

  `"mmusculus"`

  :   *Mus musculus* (Ensembl)

  `"scerevisiae"`

  :   *Saccharomyces cerevisiae* (Ensembl Fungi)

  Mutually exclusive with `mart_to_use`. Default: `NULL`.

- mart_to_use:

  Mart object or `NULL`. A mart object created with
  [`useMart`](https://huber-group-embl.github.io/biomaRt/reference/useMart.html)
  or
  [`useEnsembl`](https://huber-group-embl.github.io/biomaRt/reference/useEnsembl.html).
  Mutually exclusive with `organism`. Default: `NULL`.

## Value

A data frame with the original ninetails output data joined with the
retrieved annotation attributes via left join on transcript IDs.

## Details

Requires the biomaRt package version \>= 2.40.

Exactly one of `organism` or `mart_to_use` must be provided. The two
arguments are mutually exclusive; if both are declared the function
throws an error.

## Acknowledgements

Based on PK (smaegol) NanoTail `annotate_with_biomart()`:
<https://github.com/LRB-IIMCB/nanotail/>. Many thanks to the NanoTail
developer for help and kind advice.

## See also

[`merge_nonA_tables`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
for preparing the input,
[`getBM`](https://huber-group-embl.github.io/biomaRt/reference/getBM.html)
for the underlying biomaRt query,
[`useMart`](https://huber-group-embl.github.io/biomaRt/reference/useMart.html)
for creating custom mart objects

## Examples

``` r
if (FALSE) { # \dontrun{

# With built-in organism preset
annot <- ninetails::annotate_with_biomart(
  input_data = merged_nonA_tables,
  organism = "mmusculus"
)

# With custom mart object
mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)
annot <- ninetails::annotate_with_biomart(
  input_data = merged_nonA_tables,
  mart_to_use = mart
)

} # }
```
