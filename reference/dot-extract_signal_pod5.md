# Extract full signal from POD5 file using Python script

Internal helper function that calls the Python extraction script to
retrieve raw signal and calibration data from a POD5 file.

## Usage

``` r
.extract_signal_pod5(read_id, pod5_file, winsorize = TRUE)
```

## Arguments

- read_id:

  Character string. Read ID to extract.

- pod5_file:

  Character string. Path to POD5 file.

- winsorize:

  Logical. If TRUE, apply winsorization to signal.

## Value

List with signal data, calibration parameters, and sample_rate, or NULL
on error.
