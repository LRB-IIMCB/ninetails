# Extracts features of poly(A) tails containing only A nucleotides for training-set preparation.

Training-set variant of the feature extraction wrapper designed
exclusively for A-only signals. It processes all reads in parallel via
[`extract_tail_data_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md),
then applies an inverse filtering criterion: only reads whose pseudomove
vectors do *not* contain consecutive non-zero runs of length \>= 4 are
retained. This ensures the resulting dataset represents pure homopolymer
A tails without modification artefacts.

## Usage

``` r
create_tail_feature_list_A(
  nanopolish,
  sequencing_summary,
  workspace,
  num_cores,
  basecall_group,
  pass_only = TRUE
)
```

## Arguments

- nanopolish:

  Character string. Full path of the `.tsv` file produced by
  `nanopolish polya`.

- sequencing_summary:

  Character string. Full path of the `.txt` file with the sequencing
  summary.

- workspace:

  Character string. Full path of the directory containing basecalled
  multi-Fast5 files.

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

- basecall_group:

  Character string `["Basecall_1D_000"]`. Name of the level in the Fast5
  file hierarchy from which the data should be extracted.

- pass_only:

  Logical `[TRUE]`. If `TRUE`, only reads tagged by nanopolish as
  `"PASS"` are retained. Otherwise, reads tagged as `"PASS"` or
  `"SUFFCLIP"` are included.

## Value

A named list with one element:

- tail_feature_list:

  Named list of per-read feature lists (as returned by
  [`extract_tail_data_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md))
  containing only reads with pure A tails.

Always assign this returned list to a variable; printing the full list
to the console may crash the R session.

## Details

The inverse filtering uses a sliding-window approach
([`stats::embed`](https://rdrr.io/r/stats/embed.html)) to detect runs of
\>= 4 consecutive non-zero pseudomoves. Reads that pass this filter
(i.e. have no such runs) are collected as the A-only reference set.
Unlike
[`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md),
the returned list does not include `zeromoved_readnames` or
`nonpseudomoved_readnames` categories.

## See also

[`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md)
for the non-A variant,
[`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)
for the next pipeline step,
[`extract_tail_data_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md)
for per-read extraction,
[`prepare_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/prepare_trainingset.md)
for the top-level wrapper.

## Examples

``` r
if (FALSE) { # \dontrun{

create_tail_feature_list_A(
  nanopolish = '/path/to/file',
  sequencing_summary = '/path/to/file',
  workspace = '/path/to/guppy/workspace',
  num_cores = 10,
  basecall_group = 'Basecall_1D_000')

} # }
```
