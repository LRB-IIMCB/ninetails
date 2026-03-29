# Check and convert poly(A) length file format

Detects whether an input poly(A) length file originates from nanopolish
or tailfindr and, if necessary, converts it to the nanopolish-like
format expected by the ninetails legacy (Guppy/Fast5) pipeline.
Tailfindr cDNA output is explicitly rejected because the legacy pipeline
supports only DRS (direct RNA sequencing) data.

## Usage

``` r
check_polya_length_filetype(input)
```

## Arguments

- input:

  Character string or data frame. Either the path to a poly(A) length
  file (nanopolish `.tsv` or tailfindr `.csv`) or an in-memory data
  frame. Only DRS data are accepted for tailfindr input.

## Value

A named list with two elements:

- data:

  Data frame / tibble. Standardised poly(A) length data in
  nanopolish-like format, ready for downstream ninetails processing.

- file_type:

  Character. One of `"nanopolish"` or `"tailfindr_drs"`, indicating the
  detected input format.

## Details

File-type detection is based on column names:

- nanopolish:

  Identified by the presence of a `qc_tag` column. Returned as-is.

- tailfindr DRS:

  Identified by the presence of `read_id`, `tail_start`, `tail_end`,
  `samples_per_nt`, and `tail_length`. Converted via
  [`convert_tailfindr_output`](https://LRB-IIMCB.github.io/ninetails/reference/convert_tailfindr_output.md).

- tailfindr cDNA:

  Identified by the presence of a `tail_is_valid` column. Raises an
  error because cDNA data are not compatible with this pipeline.

If none of the above patterns match, an error is raised.

This function is part of the legacy pipeline (Guppy basecaller,
multi-Fast5 files). It may be retired if the Fast5 format is deprecated.

## See also

[`convert_tailfindr_output`](https://LRB-IIMCB.github.io/ninetails/reference/convert_tailfindr_output.md)
for the tailfindr conversion logic,
[`extract_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md)
for subsequent processing of the standardised output,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for the legacy pipeline entry point.

## Examples

``` r
if (FALSE) { # \dontrun{

# run on nanopolish output
test <- ninetails::check_polya_length_filetype(
  input = system.file('extdata', 'test_data',
                      'nanopolish_output.tsv',
                      package = 'ninetails'))

# run on tailfindr output
test <- ninetails::check_polya_length_filetype(
  input = system.file('extdata', 'test_data',
                      'tailfindr_output.csv',
                      package = 'ninetails'))

} # }
```
