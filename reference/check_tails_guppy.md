# Wrapper function for complete DRS processing by ninetails package (legacy mode).

Orchestrates the full ninetails analysis pipeline for direct RNA
sequencing (DRS) data basecalled with Guppy and characterised by
nanopolish (or tailfindr). The function validates inputs, extracts
poly(A) tail features from multi-Fast5 files, segments tails into
chunks, computes Gramian Angular Fields (GAFs), classifies reads with a
pretrained CNN, and writes output files.

## Usage

``` r
check_tails_guppy(
  polya_data,
  sequencing_summary,
  workspace,
  num_cores = 1,
  basecall_group = "Basecall_1D_000",
  pass_only = TRUE,
  qc = TRUE,
  save_dir,
  prefix = "",
  part_size = 1e+06
)
```

## Arguments

- polya_data:

  Character string or data frame. Full path of the `.tsv` file produced
  by nanopolish polya or tailfindr (DRS only), or an in-memory data
  frame. Tailfindr input is converted automatically.

- sequencing_summary:

  Character string or data frame. Full path of the `.txt` sequencing
  summary file, or an in-memory data frame.

- workspace:

  Character string. Full path of the directory containing basecalled
  multi-Fast5 files.

- num_cores:

  Numeric `[1]`. Number of physical cores to use. Do not exceed 1 less
  than the number of cores at your disposal.

- basecall_group:

  Character string `["Basecall_1D_000"]`. Name of the level in the Fast5
  file hierarchy from which data should be extracted.

- pass_only:

  Logical `[TRUE]`. If `TRUE`, only reads tagged by nanopolish as
  `"PASS"` are included. If `FALSE`, reads tagged `"PASS"` or
  `"SUFFCLIP"` are included.

- qc:

  Logical `[TRUE]`. If `TRUE`, terminal non-A residue positions likely
  arising from nanopolish segmentation errors are labelled with a
  `"-WARN"` suffix. It is up to the user whether to include or discard
  such reads.

- save_dir:

  Character string. Full path of the directory where output files should
  be stored.

- prefix:

  Character string (optional, default `""`). If provided, inserted into
  output file names between the timestamp and the file-type suffix.

- part_size:

  Numeric `[1000000]`. Maximum number of rows processed at once. Must be
  \>= 1000. If the input exceeds this value, it is split and processed
  sequentially.

## Value

A named list with two data frames:

- read_classes:

  Data frame. Per-read classification including columns `class` and
  `comments`.

- nonadenosine_residues:

  Data frame. Detailed positional information for all detected non-A
  residues.

A log file and TSV output files (`read_classes`,
`nonadenosine_residues`) are also written to `save_dir`.

## Details

**Legacy mode.** Due to Oxford Nanopore Technologies' transition from
Guppy to Dorado, from Fast5 to POD5, and from R9 to R10 chemistry, this
pipeline is maintained for backward compatibility only and will not be
further optimised. For current data, use
[`check_tails_dorado_DRS`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_DRS.md)
instead.

The function accepts either nanopolish or tailfindr output (DRS only).
Tailfindr input is automatically converted to nanopolish-like format via
[`check_polya_length_filetype`](https://LRB-IIMCB.github.io/ninetails/reference/check_polya_length_filetype.md).

For large input files (exceeding `part_size` rows), the poly(A) table is
split into parts and processed sequentially via
[`process_polya_parts`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_parts.md)
to avoid memory overflow.

The pipeline steps are:

1.  Input validation and Fast5 format checking.

2.  Feature extraction
    ([`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)).

3.  Tail segmentation
    ([`create_tail_chunk_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md)).

4.  GAF computation
    ([`create_gaf_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)).

5.  CNN classification
    ([`predict_gaf_classes`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md)).

6.  Output creation
    ([`create_outputs`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)).

The `comments` column in the output encodes detailed read
classification:

- IRL:

  Insufficient read length.

- QCF:

  Nanopolish QC failed.

- MAU:

  Move transition absent, non-A residue undetected.

- MPU:

  Move transition present, non-A residue undetected.

- NIN:

  Not included in analysis (`pass_only = TRUE`).

- YAY:

  Move transition present, non-A residue detected.

Filtering criteria: basecaller move value = 1, `qc_tag %in% c("PASS")`
(or `c("PASS", "SUFFCLIP")`), and nanopolish-estimated tail length \>=
10 nt.

## See also

[`check_tails_dorado_DRS`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_DRS.md)
for the current (Dorado-based) pipeline,
[`process_polya_complete`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_complete.md)
and
[`process_polya_parts`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_parts.md)
for the internal processing functions,
[`check_polya_length_filetype`](https://LRB-IIMCB.github.io/ninetails/reference/check_polya_length_filetype.md)
for input format detection,
[`create_outputs`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)
for the output assembly step.

## Examples

``` r
if (FALSE) { # \dontrun{

results <- ninetails::check_tails_guppy(
  polya_data = system.file('extdata', 'test_data',
                           'nanopolish_output.tsv',
                           package = 'ninetails'),
  sequencing_summary = system.file('extdata', 'test_data',
                                   'sequencing_summary.txt',
                                   package = 'ninetails'),
  workspace = system.file('extdata', 'test_data',
                          'basecalled_fast5',
                          package = 'ninetails'),
  num_cores = 2,
  basecall_group = 'Basecall_1D_000',
  pass_only = TRUE,
  qc = TRUE,
  save_dir = '~/Downloads',
  prefix = "prefix",
  part_size = 2000)

} # }
```
