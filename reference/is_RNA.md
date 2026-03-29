# Check if fast5 file contains RNA reads

Checks whether the provided fast5 file contains RNA sequencing reads by
examining the experiment_type attribute in the context_tags group.

## Usage

``` r
is_RNA(fast5_file, read_id)
```

## Arguments

- fast5_file:

  Character string. Path to the fast5 file.

- read_id:

  Character string. The read identifier within the fast5 file (e.g.,
  "read_abc123"). Can be obtained from rhdf5::h5ls() output.

## Value

Logical. TRUE if the file contains RNA reads, FALSE otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
# Check if a fast5 file contains RNA reads
fast5_path <- "path/to/file.fast5"
file_structure <- rhdf5::h5ls(fast5_path, recursive = FALSE)
read_id <- file_structure$name[1]
is_RNA(fast5_path, read_id)
} # }
```
