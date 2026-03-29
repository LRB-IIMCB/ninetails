# Check if fast5 file is multi-read format

Checks whether the provided fast5 file structure corresponds to a
multi-read fast5 file format. Multi-read fast5 files contain multiple
reads stored in groups named with "read\_" prefix.

## Usage

``` r
is_multifast5(fast5_file_structure)
```

## Arguments

- fast5_file_structure:

  A data frame returned by rhdf5::h5ls() containing the HDF5 structure
  of a fast5 file. Must contain a 'name' column.

## Value

Logical. TRUE if the file is multi-read fast5 format, FALSE otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
# Check if a fast5 file is multi-read format
fast5_path <- "path/to/file.fast5"
file_structure <- rhdf5::h5ls(fast5_path, recursive = FALSE)
is_multifast5(file_structure)
} # }
```
