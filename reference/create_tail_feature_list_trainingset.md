# Extracts features of poly(A) tails of ONT RNA reads required for finding non-A nucleotides within the given tails.

This is the training-set variant of the feature extraction wrapper. It
processes all reads in parallel, calling
[`extract_tail_data_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md)
for each read, then filters out reads with zero-moved tails or
pseudomove chains too short to represent potential modifications.

## Usage

``` r
create_tail_feature_list_trainingset(
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

A named list with three elements:

- tail_feature_list:

  Named list of per-read feature lists (as returned by
  [`extract_tail_data_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md)).

- zeromoved_readnames:

  Character vector. Read IDs discarded because their tail moves summed
  to zero.

- nonpseudomoved_readnames:

  Character vector. Read IDs discarded because their pseudomove chains
  were too short (\< 4).

Always assign this returned list to a variable; printing the full list
to the console may crash the R session.

## Details

The function differs from its production counterpart in that it retains
reads whose pseudomove chains satisfy a length \>= 4 criterion, which is
required for subsequent modification-centered chunk splitting. Two
categories of discarded reads are tracked (zero-moved and
non-pseudomoved) and returned alongside the valid feature list.

## See also

[`extract_tail_data_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md)
for the per-read extraction step,
[`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
for the next pipeline step,
[`extract_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md)
for input data preparation.

## Examples

``` r
if (FALSE) { # \dontrun{

create_tail_feature_list_trainingset(
  nanopolish = '/path/to/file',
  sequencing_summary = '/path/to/file',
  workspace = '/path/to/guppy/workspace',
  num_cores = 10,
  basecall_group = 'Basecall_1D_000')

} # }
```
