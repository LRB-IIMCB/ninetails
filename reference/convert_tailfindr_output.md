# Converts tailfindr results to format compatible with ninetails

Reformats the output of the tailfindr pipeline so it can be passed to
the ninetails legacy (Guppy) pipeline in place of nanopolish polya
output. The function renames key columns, creates derived columns, and
introduces dummy quality tags to approximate the nanopolish output
schema.

## Usage

``` r
convert_tailfindr_output(tailfindr_output)
```

## Arguments

- tailfindr_output:

  Character string or data frame. Either the full path of the `.csv`
  file produced by tailfindr, or an in-memory data frame containing
  tailfindr result data. Required columns: `read_id`, `tail_start`,
  `tail_end`, `tail_length`.

## Value

A tibble containing tailfindr results reformatted to resemble the output
of `nanopolish polya`, with columns:

- readname:

  Character. Read identifier (renamed from `read_id`).

- polya_start:

  Integer. Start coordinate of the poly(A) tail (renamed from
  `tail_start`).

- tail_end:

  Integer. End coordinate of the poly(A) tail (retained from tailfindr).

- polya_length:

  Numeric. Estimated tail length (renamed from `tail_length`).

- transcript_start:

  Integer. Transcript start coordinate (`tail_end + 1`).

- contig:

  Character. Dummy mapping target (`"tailfindr_out"`).

- qc_tag:

  Character. Quality tag (`"PASS"` or `"NOREGION"`).

Always assign this returned tibble to a variable; printing the full
tibble to the console may crash the R session.

## Details

The following column mappings are applied:

- `read_id` \\\rightarrow\\ `readname`:

  Read identifier renamed to match the ninetails naming convention.

- `tail_start` \\\rightarrow\\ `polya_start`:

  Start coordinate of the poly(A) tail.

- `tail_length` \\\rightarrow\\ `polya_length`:

  Estimated poly(A) tail length.

The following columns are created:

- `transcript_start`:

  Set to `tail_end + 1`.

- `contig`:

  Dummy value `"tailfindr_out"` (nanopolish returns mapping information,
  which is absent in tailfindr output).

- `qc_tag`:

  Simplified quality tag. Because the R9.4.1 pore detection region spans
  5 nucleotides, tails shorter than 10 nt (2 full adjacent 5-mers with
  no overlap) are unreliable. Therefore `qc_tag` is set to `"PASS"` for
  tails \>= 10 nt and `"NOREGION"` otherwise.

## Warning

Ninetails is optimised to work with nanopolish polya output. The
HMM-based approach of nanopolish provides more robust tail boundary
predictions than the slope estimator used by tailfindr. In particular,
ninetails relies on the quality tags produced by nanopolish (`"PASS"` /
`"SUFFCLIP"`). When using tailfindr output, the signal quality cannot be
inferred, which may lead to poor-quality signals being passed to the CNN
and consequently misclassified. **Use tailfindr input at your own
risk.**

## See also

[`check_polya_length_filetype`](https://LRB-IIMCB.github.io/ninetails/reference/check_polya_length_filetype.md)
which uses this function internally,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for the legacy pipeline that accepts the converted output as its
`nanopolish` argument,
[`extract_polya_data`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md)
for how nanopolish output is normally processed.

## Examples

``` r
if (FALSE) { # \dontrun{

df <- ninetails::convert_tailfindr_output(
  tailfindr_output = '/path/to/tailfindr_out.csv')

# The output can be passed to the ninetails pipeline as the
# nanopolish argument:
results <- ninetails::check_tails(
  nanopolish = df,
  sequencing_summary = '/path/to/sequencing_summary.txt',
  workspace = '/path/to/workspace',
  num_cores = 2,
  basecall_group = 'Basecall_1D_000',
  pass_only = TRUE,
  save_dir = '~/Downloads')

} # }
```
