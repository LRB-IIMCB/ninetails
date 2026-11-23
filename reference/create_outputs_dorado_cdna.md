# Create ninetails output tables for the cDNA pipeline

Combines the tail coordinates with the non-adenosine calls of
`predict_nonA_cdna` into the two standard ninetails outputs: per-read
classes and non-adenosine residues with their estimated positions.
Inputs are read from the intermediate directories written by
[`process_polya_reads_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_reads_cdna.md)
or
[`process_polyt_reads_cdna`](https://LRB-IIMCB.github.io/ninetails/reference/process_polyt_reads_cdna.md).

## Usage

``` r
create_outputs_dorado_cdna(
  dorado_summary_dir,
  nonA_temp_dir,
  polya_chunks_dir,
  num_cores = 1,
  qc = TRUE
)
```

## Arguments

- dorado_summary_dir:

  Character string. Path to a directory with the tail coordinate tables
  (.txt, .tsv or .csv) of the reads of one orientation. Required
  columns: `read_id`, `contig`, `polya_length`, `qc_tag` and
  `tail_type`.

- nonA_temp_dir:

  Character string. Path to a directory with the per-part RDS files of
  non-adenosine calls. Each file holds a list with `nonA_calls` (output
  of `predict_nonA_cdna`) and `classified_reads` (IDs of the reads whose
  signals were classified).

- polya_chunks_dir:

  Character string. Kept for compatibility with the DRS pipeline
  (currently unused; the cDNA model does not produce chunk files).

- num_cores:

  Integer \[1\]. Number of cores used for loading the files.

- qc:

  Logical \[TRUE\]. Whether to discard terminal calls and reads failing
  the Lachesis QC.

## Value

A named list:

- read_classes:

  Data frame with columns `readname`, `contig`, `polya_length`,
  `qc_tag`, `class`, `comments`, `tail_type`

- nonadenosine_residues:

  Data frame with columns `readname`, `contig`, `prediction` (residue in
  the RNA tail), `est_nonA_pos` (position from the 3' end of the tail),
  `polya_length`, `qc_tag`, `tail_type`

## Read classes

- `decorated`, `YAY`: at least one non-adenosine residue

- `blank`, `MAU`: no non-adenosine residue detected

- `blank`, `MPU`: only terminal calls, discarded by `qc` (closer than 2
  nt to either end of the tail)

- `unclassified`, `IRL`: tail shorter than 10 nt

- `unclassified`, `BAC`: invalid tail coordinates or no signal available

- `unclassified`, `QCF`: Lachesis QC tag other than `PASS` (only with
  `qc = TRUE`)

## See also

`predict_nonA_cdna`,
[`merge_cdna_results`](https://LRB-IIMCB.github.io/ninetails/reference/merge_cdna_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
outputs <- ninetails::create_outputs_dorado_cdna(
  dorado_summary_dir = "polya_processing_dir/dorado_summary_dir",
  nonA_temp_dir = "polya_processing_dir/nonA_temp_dir",
  polya_chunks_dir = "polya_processing_dir/polya_chunks_dir",
  num_cores = 4,
  qc = TRUE
)
} # }
```
