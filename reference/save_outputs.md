# Save pipeline outputs to files.

Writes the ninetails pipeline output data frames (`read_classes` and
`nonadenosine_residues`) to tab-separated files in `save_dir`. File
names include a timestamp and an optional user-defined prefix.

## Usage

``` r
save_outputs(outputs, save_dir, prefix = "")
```

## Arguments

- outputs:

  Named list containing data frames to save. Expected names are
  `read_classes` and `nonadenosine_residues`.

- save_dir:

  Character string. Full path of the directory where output files should
  be stored.

- prefix:

  Character string (optional, default `""`). If provided, inserted into
  file names between the timestamp and the file-type suffix.

## Value

Invisible `NULL`. Called for the side effect of writing files.

## Details

This function is called internally by
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
at the end of the pipeline. It is not intended to be called directly by
the user.

## See also

[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
which calls this function.
