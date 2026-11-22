# Split large poly(A) data file into smaller parts.

Divides a poly(A) length table (nanopolish or tailfindr output) into
smaller files of at most `part_size` rows each, saving them to a
`polya_data_parts` subdirectory within `save_dir`. This prevents memory
overflow when processing large datasets.

## Usage

``` r
split_polya_data(polya_data, part_size = 1e+05, save_dir)
```

## Arguments

- polya_data:

  Character string or data frame. Full path of the poly(A) length file
  (`.tsv`), or an in-memory data frame.

- part_size:

  Numeric `[100000]`. Maximum number of rows per output part file.

- save_dir:

  Character string. Full path of the directory where the
  `polya_data_parts` subdirectory will be created.

## Value

Character vector of file paths to the created part files, named
`polya_data_part_<i>_of_<n>.tsv`.

## Details

This function is called internally by
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
when the input exceeds the `part_size` threshold. It is not intended to
be called directly by the user.

## See also

[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
which calls this function,
[`process_polya_parts`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_parts.md)
which processes the resulting parts.
