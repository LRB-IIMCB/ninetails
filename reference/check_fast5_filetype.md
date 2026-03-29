# Check if the provided directory contains Fast5 files in the correct format

Analyzes the structure of the first Fast5 file in the given directory
and checks whether it fulfills the analysis requirements: the file must
be a multi-Fast5, basecalled by Guppy basecaller, and contain the
provided basecall group and RNA reads. Otherwise the function throws a
descriptive error.

## Usage

``` r
check_fast5_filetype(workspace, basecall_group)
```

## Arguments

- workspace:

  Character string. Full path of the directory containing the basecalled
  multi-Fast5 files.

- basecall_group:

  Character string. Name of the level in the Fast5 file hierarchy from
  which data should be extracted (e.g., `"Basecall_1D_000"`).

## Value

Prints to console: data type (RNA), Fast5 file type (multi-Fast5),
basecaller used, basecaller version, and basecalling model. Called for
its side effect; returns invisibly.

## Acknowledgements

This lookup function is inspired by adnaniazi's
`explore-basecaller-and-fast5type.R` from tailfindr:
<https://github.com/adnaniazi/tailfindr/blob/master/R/explore-basecaller-and-fast5type.R>.

## See also

[`is_multifast5`](https://LRB-IIMCB.github.io/ninetails/reference/is_multifast5.md)
for the multi-Fast5 format check,
[`is_RNA`](https://LRB-IIMCB.github.io/ninetails/reference/is_RNA.md)
for the RNA content check,
[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
where this function is called during pipeline initialization

## Examples

``` r
if (FALSE) { # \dontrun{

check_fast5_filetype(
  workspace = '/path/to/guppy/workspace',
  basecall_group = 'Basecall_1D_000'
)

} # }
```
