# Create list of poly(A) tail features from multi-Fast5 files

Extracts tail features of RNA reads from multi-Fast5 files basecalled by
Guppy and poly(A) tail characteristics (coordinates) produced by
nanopolish polya function. Processing is parallelized across reads using
foreach and doSNOW. A progress bar is displayed during extraction.

## Usage

``` r
create_tail_feature_list(
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

  Character string or data frame. Full path of the `.tsv` file produced
  by nanopolish polya function, or an in-memory data frame.

- sequencing_summary:

  Character string or data frame. Full path of the `.txt` file with
  sequencing summary, or an in-memory data frame.

- workspace:

  Character string. Full path of the directory containing the basecalled
  multi-Fast5 files.

- num_cores:

  Numeric. Number of physical cores to use in processing. Do not exceed
  1 less than the number of cores at your disposal.

- basecall_group:

  Character string. Name of the level in the Fast5 file hierarchy from
  which data should be extracted (e.g., `"Basecall_1D_000"`).

- pass_only:

  Logical. If `TRUE` (default), only reads tagged by nanopolish as
  `"PASS"` are taken into consideration. If `FALSE`, reads tagged as
  `"PASS"` and `"SUFFCLIP"` are both included.

## Value

A named list with three elements:

- tail_feature_list:

  Named list of per-read tail features. Each element contains
  `fast5_filename`, `tail_signal`, `tail_moves`, and `tail_pseudomoves`
  (see
  [`extract_tail_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data.md)).

- zeromoved_readnames:

  Character vector. Read IDs discarded because all basecaller moves in
  their tail region were zero.

- nonpseudomoved_readnames:

  Character vector. Read IDs discarded because their pseudomove chain
  was too short (\< 5 consecutive positions) to indicate a potential
  modification.

Always assign the returned list to a variable. Printing the full output
to the console may crash your R session.

## Details

After extraction, reads with zero-moved tails and reads that do not
satisfy the pseudomove condition (minimum run length of 5) are filtered
out and their identifiers are stored separately for downstream
classification.

## See also

[`extract_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md)
for reading nanopolish and sequencing summary data,
[`extract_tail_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data.md)
for single-read extraction,
[`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md)
for downstream chunk segmentation

## Examples

``` r
if (FALSE) { # \dontrun{

tfl <- ninetails::create_tail_feature_list(
  nanopolish = system.file('extdata',
                           'test_data',
                           'nanopolish_output.tsv',
                           package = 'ninetails'),
  sequencing_summary = system.file('extdata',
                                   'test_data',
                                   'sequencing_summary.txt',
                                   package = 'ninetails'),
  workspace = system.file('extdata',
                          'test_data',
                          'basecalled_fast5',
                          package = 'ninetails'),
  num_cores = 2,
  basecall_group = 'Basecall_1D_000',
  pass_only = TRUE
)

} # }
```
