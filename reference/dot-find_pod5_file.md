# Find POD5 file in directory based on filename

Internal helper function to locate a POD5 file within a directory.
Searches recursively and matches by filename if provided.

## Usage

``` r
.find_pod5_file(filename, pod5_dir)
```

## Arguments

- filename:

  Character string. Filename from dorado summary (optional).

- pod5_dir:

  Character string. Directory containing POD5 files.

## Value

Full path to POD5 file or NULL if not found.
