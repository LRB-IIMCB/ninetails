# Filter Dorado summary for reads fulfilling ninetails quality criteria

Takes a Dorado summary file (from ONT Dorado basecaller) or a data frame
containing equivalent summary information, and filters out reads that do
not meet ninetails quality and alignment criteria.

## Usage

``` r
filter_dorado_summary(dorado_summary)
```

## Arguments

- dorado_summary:

  Character string or data frame. Path to a tab-delimited Dorado summary
  file, or a data frame containing equivalent summary information. Must
  contain at least the columns `read_id` and `poly_tail_length`.

## Value

A filtered tibble or data frame containing only reads that meet the
filtering criteria.

## Details

Reads must meet all four criteria to pass the filter:

1.  **Mapped**: `alignment_direction` is `"+"` or `"-"`, not `"*"`

2.  **High mapping quality**: `alignment_mapq` \> 0

3.  **Valid coordinates**: `poly_tail_start` != 0 (in DRS, the adapter
    passes through the pore first, so a start position of 0 indicates a
    likely artifact)

4.  **Sufficient length**: `poly_tail_length` \>= 10 nt (shorter tails
    cannot be reliably processed by the CNN)

## See also

[`preprocess_inputs`](https://LRB-IIMCB.github.io/ninetails/reference/preprocess_inputs.md)
and
[`preprocess_inputs_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/preprocess_inputs_cdna.md)
where this function is called during pipeline initialization,
[`process_dorado_summary`](https://LRB-IIMCB.github.io/ninetails/reference/process_dorado_summary.md)
for splitting large summary files

## Examples

``` r
if (FALSE) { # \dontrun{

# From file
filtered <- filter_dorado_summary("dorado_summary.txt")

# From data frame
df <- data.frame(
  read_id = c("read1", "read2"),
  alignment_direction = c("+", "*"),
  alignment_mapq = c(60, 0),
  poly_tail_start = c(100, 0),
  poly_tail_length = c(20, 5)
)
filtered <- filter_dorado_summary(df)

} # }
```
