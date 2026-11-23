# Complete Oxford Nanopore poly(A)/poly(T) tail analysis pipeline for Dorado cDNA data.

This function extends the check_tails_dorado_DRS pipeline to handle cDNA
sequencing data by adding BAM file processing for sequence extraction
and Dorado-style read orientation classification. The pipeline
identifies and characterizes non-adenosine nucleotides within poly(A)
and poly(T) tails using the same signal processing and machine learning
approach as the DRS pipeline, but with automatic classification of read
orientations.

## Usage

``` r
check_tails_dorado_cDNA(
  bam_file,
  dorado_summary,
  pod5_dir,
  num_cores = 1,
  qc = TRUE,
  save_dir,
  prefix = "",
  part_size = 40000,
  cleanup = FALSE,
  lachesis = NULL
)
```

## Arguments

- bam_file:

  Character string. Path to the BAM file containing aligned cDNA reads
  with basecalled sequences. This file will be split into parts for
  memory management.

- dorado_summary:

  Character string or data frame. Path to Dorado summary file or data
  frame containing per-read summary information. Must include standard
  columns such as read_id, filename, etc.

- pod5_dir:

  Character string. Path to directory containing POD5 files with raw
  nanopore signal data corresponding to the reads in the summary file.

- num_cores:

  Integer \[1\]. Number of CPU cores to use for parallel processing.
  Recommend using \`parallel::detectCores() - 1\` for optimal
  performance while maintaining system responsiveness.

- qc:

  Logical \[TRUE\]. Whether to apply quality control filtering during
  analysis. When TRUE, applies standard ninetails QC filters including
  tail length filtering and coordinate validation.

- save_dir:

  Character string. Path to directory where all output files and
  intermediate results will be saved. Directory will be created if it
  doesn't exist.

- prefix:

  Character string \[""\]. Optional prefix to add to all output file
  names. Useful for distinguishing between different experimental
  conditions or samples.

- part_size:

  Integer \[40000\]. Number of reads to process in each chunk when
  splitting large input files. Larger values use more memory but may be
  faster. Adjust based on available system memory.

- cleanup:

  Logical \[FALSE\]. Whether to remove intermediate files after
  successful pipeline completion. When FALSE, all intermediate files are
  preserved for inspection.

- lachesis:

  Character string or data frame \[NULL\]. Path to the Lachesis output
  file (or the same table in memory). When provided, the tail
  coordinates and read orientation (`strand` column) are based on
  Lachesis predictions. If that is the case, `bam_file` and
  `dorado_summary` are not required.

## Value

A named list containing the final analysis results:

- read_classes:

  Data frame with per-read classification results including readname,
  contig, poly(A) length, QC tag, class, comments, and tail_type

- nonadenosine_residues:

  Data frame with predicted non-adenosine positions within
  poly(A)/poly(T) tails including readname, contig, prediction,
  estimated position, poly(A) length, QC tag, and tail_type. Predictions
  refer to the residue present in the original RNA tail (C, G or U) and
  positions are counted from the 3' end of the tail, in both
  orientations

## Pipeline Overview

1.  **Tail coordinates and orientation**: taken from Lachesis
    (`lachesis` param), or from Dorado summary with orientation
    classified from the basecalled sequences in the BAM file (based on
    SSP/VNP primer matching)

2.  **Signal extraction**: raw tail signals with flanks are extracted
    from POD5 files

3.  **Classification**: sliding windows of each tail are classified with
    the poly(A) or poly(T) cDNA model (see `predict_nonA_cdna`)

4.  **Output with Tail Types**: Produces standard read_classes and
    nonadenosine_residues with tail_type information

## Input Requirements

This pipeline requires specific input formats:

- **Lachesis output**: tail coordinates and read orientation from
  Lachesis software (recommended) or

- **Dorado Summary**: Must contain standard columns for read information
  (`--estimate-poly-a`; Dorado \>=1.0.0)

- **BAM File**: Aligned cDNA reads with basecalled, untrimmed sequences
  (`--no-trim`)

- **POD5 Files**: Raw signal files corresponding to reads in summary

## Key Differences from DRS Pipeline

- **BAM Input**: Additional BAM file input for sequence extraction since
  cDNA data requires basecalled sequences

- **Dorado-Style Read Orientation Classification**: matching of SSP/VNP
  primers to classify reads as polyA, polyT, or unidentified before
  processing

- **Separate models**: poly(A) and poly(T) reads are classified with
  dedicated cDNA models working on the raw (not downsampled, nor
  winsorized) signal

- **RNA-level reporting**: predictions refer to the residues of the
  original RNA tail in both orientations

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic cDNA analysis
results <- ninetails::check_tails_dorado_cDNA(
  bam_file = "path/to/aligned_cdna.bam",
  dorado_summary = "path/to/dorado_summary.txt",
  pod5_dir = "path/to/pod5_files/",
  num_cores = 4,
  save_dir = "path/to/output/"
)

# Access results
head(results$read_classes)
head(results$nonadenosine_residues)

# Analysis with custom settings
results <- ninetails::check_tails_dorado_cDNA(
  bam_file = "large_dataset.bam",
  dorado_summary = summary_df,  # Can pass data frame
  pod5_dir = "/data/pod5/",
  num_cores = 8,
  qc = TRUE,
  save_dir = "/results/experiment1/",
  prefix = "exp1_sample_A",
  part_size = 20000,  # Smaller chunks for limited memory
  cleanup = TRUE,  # Remove intermediate files
  lachesis = NULL # use Dorado input
)
} # }
```
