# Extract poly(A) tail signal segments from POD5 files using parallel Python processing

This function extracts raw nanopore signal data corresponding to poly(A)
tail regions from POD5 format files. It leverages a Python subprocess
for efficient parallel extraction, avoiding R-Python serialization
bottlenecks inherent in reticulate. The function processes multiple
reads across multiple POD5 files simultaneously, applies winsorization
to reduce noise, and performs signal interpolation to standardize tail
representations.

## Usage

``` r
extract_tails_from_pod5(polya_data, pod5_dir, num_cores = 1)
```

## Arguments

- polya_data:

  Data frame containing poly(A) tail coordinate information. Must
  include the following columns:

  read_id

  :   Character. Unique identifier for each nanopore read

  filename

  :   Character. Name of the POD5 file containing the read

  poly_tail_start

  :   Integer. Start position of poly(A) tail in the raw signal

  poly_tail_end

  :   Integer. End position of poly(A) tail in the raw signal

  poly_tail_length

  :   Integer. (Optional) Length of the poly(A) tail in bases. If
      present, only tails \> 10 bases are processed

- pod5_dir:

  Character string. Full path to the directory containing POD5 files.
  All POD5 files referenced in the `filename` column of `polya_data`
  must be present in this directory.

- num_cores:

  Integer. Number of CPU cores to use for parallel processing. Default
  is 1. If set to NULL, the function automatically uses all available
  cores minus 1 to maintain system responsiveness. Values \> 1 enable
  parallel processing of POD5 files.

## Value

A named list of numeric vectors, where:

- Names correspond to read_id values from the input

- Each vector contains the winsorized and interpolated signal values

- Empty vectors indicate reads where extraction failed

Returns an empty list if no valid reads are found.

## Details

The function performs the following operations:

1.  **Input validation**: Checks for required columns and valid data

2.  **Read filtering**: Excludes reads with:

    - poly_tail_length â‰¤ 10 (if column exists)

    - poly_tail_start â‰¤ 0 (invalid coordinates)

3.  **Python subprocess execution**: Delegates extraction to an
    optimized Python script that:

    - Groups reads by POD5 file for efficient I/O

    - Processes files in parallel using multiprocessing

    - Extracts signal segments based on provided coordinates

4.  **Signal processing**: For each valid tail region:

    - Applies winsorization (0.5

    - Interpolates to 20

    - Converts to integer values for downstream compatibility

The Python subprocess approach bypasses reticulate's Global Interpreter
Lock (GIL) limitations, providing true parallel processing for large
datasets. Temporary files are used for data exchange and automatically
cleaned up after processing.

## Note

- Requires Dorado \>=1.0.0 to retrieve polyA coordinates

- Requires Python 3.6+ with the 'pod5' module installed:
  `pip install pod5`

- The POD5 extraction script must be present at
  `system.file("extdata", "extract_pod5_signals.py", package = "ninetails")`

- Large datasets may require substantial temporary disk space

- Progress messages are printed directly from the Python subprocess

## See also

[`create_tail_features_list_dorado`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_features_list_dorado.md)
for computing pseudomoves from signals,
[`process_dorado_signal_files`](https://LRB-IIMCB.github.io/ninetails/reference/process_dorado_signal_files.md)
for complete signal processing pipeline

## Examples

``` r
if (FALSE) { # \dontrun{
# Load poly(A) coordinate data
polya_data <- read.table("polya_coords.txt", header = TRUE)

# Extract signals using 4 cores
signals <- extract_tails_from_pod5(
  polya_data = polya_data,
  pod5_dir = "/path/to/pod5/files/",
  num_cores = 4
)} # }
```
