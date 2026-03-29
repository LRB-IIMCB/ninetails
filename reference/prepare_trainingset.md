# Filters out signals of a given nucleotide type for neural network training-set preparation.

Top-level convenience wrapper that orchestrates the complete
training-set production pipeline for a single nucleotide category.
Depending on the selected `nucleotide`, it chains the appropriate
feature extraction, chunk splitting, filtering, and GAF creation
functions into a single call.

## Usage

``` r
prepare_trainingset(
  nucleotide,
  nanopolish,
  sequencing_summary,
  workspace,
  num_cores = 1,
  basecall_group = "Basecall_1D_000",
  pass_only = TRUE
)
```

## Arguments

- nucleotide:

  Character. One of `"A"`, `"C"`, `"G"`, or `"U"`. Defines the type of
  signal filtering applied to produce training data for the desired
  nucleotide context.

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

A named list of GAF matrices organised by `<read_ID>_<index>`. Always
assign this returned list to a variable; printing the full list to the
console may crash the R session.

## Details

The internal pipeline differs by nucleotide:

- `"A"`:

  Uses the A-only branch:
  [`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md)
  \\\rightarrow\\
  [`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)
  \\\rightarrow\\
  [`create_gaf_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list_A.md).

- `"C"`:

  Uses the non-A branch with `value = -1` (valley filtering):
  [`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md)
  \\\rightarrow\\
  [`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
  \\\rightarrow\\
  [`filter_nonA_chunks_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md)
  \\\rightarrow\\
  [`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md).

- `"G"`:

  Uses the non-A branch with `value = 1` (peak filtering).

- `"U"`:

  Uses the non-A branch with `value = -1` (valley filtering), same
  filter direction as C.

## See also

[`create_tail_feature_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md),
[`create_tail_feature_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md),
[`create_tail_chunk_list_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md),
[`create_tail_chunk_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md),
[`filter_nonA_chunks_trainingset`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md),
[`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md),
[`create_gaf_list_A`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list_A.md).

## Examples

``` r
if (FALSE) { # \dontrun{

prepare_trainingset(
  nucleotide = "A",
  nanopolish = '/path/to/file',
  sequencing_summary = '/path/to/file',
  workspace = '/path/to/guppy/workspace',
  num_cores = 10,
  basecall_group = 'Basecall_1D_000',
  pass_only = TRUE)

} # }
```
