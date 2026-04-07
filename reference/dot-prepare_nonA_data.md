# Prepare non-A data for overlay on signal plots

Loads and filters the nonadenosine_residues table for a single read.
Handles both file paths and data frames, and normalizes the read ID
column name (`read_id` vs `readname`).

## Usage

``` r
.prepare_nonA_data(residue_data, readname)
```

## Arguments

- residue_data:

  Character string (file path) or data frame. The nonadenosine_residues
  table.

- readname:

  Character string. The read ID to filter for.

## Value

Data frame with non-A data for the specified read, or NULL if no
matching rows found.
