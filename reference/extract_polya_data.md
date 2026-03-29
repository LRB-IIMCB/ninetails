# Extract poly(A) data from nanopolish output and sequencing summary

Extracts features of poly(A) tails of selected RNA reads from the output
table provided by nanopolish polya function and the sequencing summary
provided by the sequencer. Filenames are taken from the sequencing
summary file. Only reads with tail lengths estimated as \>= 10 nt by
nanopolish polya function are taken into account.

## Usage

``` r
extract_polya_data(nanopolish, sequencing_summary, pass_only = TRUE)
```

## Arguments

- nanopolish:

  Character string or data frame. Either the full path of the `.tsv`
  file produced by nanopolish polya function or an in-memory data frame
  containing nanopolish data.

- sequencing_summary:

  Character string or data frame. Either the full path of the `.txt`
  file with sequencing summary or an in-memory data frame containing
  sequencing summary data.

- pass_only:

  Logical. If `TRUE` (default), only reads tagged by nanopolish as
  `"PASS"` are taken into consideration. If `FALSE`, reads tagged as
  `"PASS"` and `"SUFFCLIP"` are both included in the analysis.

## Value

A data frame containing read information organized by the read ID.
Columns include:

- readname:

  Character. Read identifier

- polya_start:

  Integer. Start position of the poly(A) tail in the raw signal

- transcript_start:

  Integer. Start position of the transcript in the raw signal

- polya_length:

  Numeric. Estimated poly(A) tail length in nucleotides

- qc_tag:

  Character. Nanopolish quality control tag

- filename:

  Character. Name of the source Fast5 file

Always assign the returned data frame to a variable. Printing the full
output to the console may crash your R session.

## Details

The function performs the following operations:

1.  Reads and validates nanopolish and sequencing summary inputs
    (accepts both file paths and in-memory data frames)

2.  Filters reads by QC tag (`pass_only` parameter)

3.  Joins nanopolish poly(A) data with sequencing summary by read name

4.  Filters reads with poly(A) tail length \>= 10 nt

5.  Removes duplicate entries from secondary alignments

## See also

[`extract_tail_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data.md)
for extracting tail features from individual reads,
[`create_tail_feature_list`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)
for batch feature extraction

## Examples

``` r
if (FALSE) { # \dontrun{

ninetails::extract_polya_data(
  nanopolish = '/path/to/nanopolish/polya/output.tsv',
  sequencing_summary = '/path/to/sequencing_summary.txt',
  pass_only = TRUE
)

} # }
```
