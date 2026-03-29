# Split BAM file into parts based on read IDs from summary file

This function splits a large BAM file into smaller parts based on read
IDs from corresponding Dorado summary files. This is essential for
memory management when processing large cDNA datasets. The function
filters the BAM file to include only reads present in the summary file
and creates appropriately sized output files.

## Usage

``` r
split_bam_file_cdna(
  bam_file,
  dorado_summary,
  part_size = 1e+05,
  save_dir,
  part_number,
  cli_log = message
)
```

## Arguments

- bam_file:

  Character string. Path to input BAM file to be split.

- dorado_summary:

  Character string. Path to corresponding Dorado summary file containing
  read IDs to include in this part.

- part_size:

  Integer. Target number of reads per output file part.

- save_dir:

  Character string. Directory where split BAM files will be saved.

- part_number:

  Integer. Part number for naming output files.

- cli_log:

  Function for logging messages and progress.

## Value

Character vector of output BAM file paths created.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_files <- split_bam_file_cdna(
  bam_file = "large_dataset.bam",
  dorado_summary = "summary_part1.txt",
  part_size = 40000,
  save_dir = "bam_parts/",
  part_number = 1,
  cli_log = message
)
} # }
```
