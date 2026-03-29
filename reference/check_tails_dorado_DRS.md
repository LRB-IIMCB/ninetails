# Complete Oxford Nanopore poly(A) tail analysis pipeline for Dorado DRS data.

This comprehensive wrapper function orchestrates the complete analysis
of Oxford Nanopore direct RNA sequencing (DRS) data processed with
Dorado basecaller (\>=1.0.0) using POD5 file format. The pipeline
identifies and characterizes non-adenosine nucleotides within poly(A)
tails through advanced signal processing, machine learning-based
classification, and statistical analysis.

## Usage

``` r
check_tails_dorado_DRS(
  dorado_summary,
  pod5_dir,
  num_cores = 1,
  qc = TRUE,
  save_dir,
  prefix = "",
  part_size = 40000,
  cleanup = FALSE
)
```

## Arguments

- dorado_summary:

  Character string or data frame. Either the full path to a Dorado
  summary file (.txt, .tsv, .csv) or an in-memory data frame containing
  summary information. **Required columns**: `read_id` (unique
  identifier), `filename` (POD5 file name), `poly_tail_length` (tail
  length in nucleotides), `poly_tail_start` (start coordinate),
  `poly_tail_end` (end coordinate). Additional columns such as
  `alignment_genome`, `alignment_mapq` are automatically included if
  present.

- pod5_dir:

  Character string. Full path to the directory containing POD5 files.
  The directory should contain POD5 files referenced in the `filename`
  column of the Dorado summary. Files can be organized in subdirectories
  (recursive search is performed). Typical ONT directory structures are
  automatically handled.

- num_cores:

  Integer \[1\]. Number of physical CPU cores to use for parallel
  processing. Should not exceed `parallel::detectCores() - 1` to
  maintain system responsiveness. Higher core counts significantly
  reduce processing time for large datasets. Memory usage scales
  approximately linearly with core count.

- qc:

  Logical \[TRUE\]. Enable comprehensive quality control filtering. When
  `TRUE`, applies stringent filters including:

  - Minimum poly(A) tail length (10 nucleotides)

  - Coordinate validation and sanitization

  - Terminal position masking (reduces false positives)

  - Statistical outlier detection and removal

  Set to `FALSE` only for preliminary analysis or when using
  pre-filtered data.

- save_dir:

  Character string. Full path to the output directory where all results,
  intermediate files, and logs will be saved. The directory will be
  created if it doesn't exist. If the directory contains existing files,
  the user will be prompted to choose whether to abort or overwrite.
  Requires sufficient disk space (depending on dataset size).

- prefix:

  Character string \[""\]. Optional prefix for output filenames. When
  provided, this string is prepended to all output files between the
  timestamp and file type identifier. Useful for organizing multiple
  analyses or experimental conditions. Should contain only
  filesystem-safe characters (alphanumeric, underscore, hyphen).

- part_size:

  Integer \[40000\]. Maximum number of reads to process in each file
  partition. Must be \>= 1. Larger values increase memory usage but may
  improve processing efficiency. Smaller values reduce memory footprint
  and enable processing of very large datasets on memory-constrained
  systems. Optimal values typically range from 10,000 to 100,000
  depending on available RAM.

- cleanup:

  Logical \[FALSE\]. Controls removal of intermediate files after
  successful analysis completion. When `TRUE`, removes all temporary
  subdirectories (`dorado_summary_dir`, `polya_signal_dir`,
  `nonA_temp_dir`, `polya_chunks_dir`) keeping only final results and
  log files. When `FALSE`, preserves all intermediate files for
  debugging or detailed inspection. Recommended to keep `TRUE` to save
  disk space.

## Value

A named list containing comprehensive analysis results:

- `read_classes`:

  Data frame with per-read classification results including:

  - `readname`: Unique read identifier

  - `contig`: Reference genome/transcript name (if aligned)

  - `polya_length`: Estimated poly(A) tail length

  - `qc_tag`: Quality control status

  - `class`: Classification result ("decorated", "blank",
    "unclassified")

  - `comments`: Detailed classification rationale using standard codes

- `nonadenosine_residues`:

  Data frame with detailed information about detected non-A nucleotides
  including:

  - `readname`: Read identifier

  - `contig`: Reference information

  - `prediction`: Predicted nucleotide type (C, G, U)

  - `est_nonA_pos`: Estimated position within poly(A) tail

  - `polya_length`: Total tail length

  - `qc_tag`: Quality metrics

## Details

Reads are classified into three main categories with specific biological
interpretations. The classification system is hierarchical, with quality
control failures taking precedence over biological interpretation.

The Dorado and Guppy pipelines share the same core algorithm but differ
in several technical details that affect output interpretation.

- Dorado:

  Summary files include alignment_direction field explicitly. Unmapped
  reads have direction = "\*"

- Dorado:

  Integer positions: `round(position, 0)`. Reflects the discrete nature
  of nucleotide positions unlike nanopolish-based predictions

- Dorado:

  Does provides move data in BAM files, however iteration through them
  is time-consuming, so it was omitted in dorado pipelines. Pseudomoves
  are computed from raw signal using z-score peak detection algorithm.

- Dorado:

  BAC code specifically checks poly_tail_start = 0. In DRS, this is
  definitively an error as adapter sequences precede poly(A) signal

Despite these technical differences, both pipelines produce compatible
output tables with identical column names and interpretable values.
Users can compare results across pipelines by accounting for the integer
vs decimal position difference. Classification categories
(decorated/blank/unclassified) have identical biological meanings across
both pipelines.

**Recommendation:** Use Dorado pipeline for new analyses (modern format,
actively maintained). Use Guppy pipeline only for legacy data or when
POD5 conversion is not feasible.

## Note

**Important considerations**:

- Ensure sufficient disk space (typically 2-5x input size) in `save_dir`

- The function generates detailed log files for troubleshooting

- Results should always be assigned to a variable to prevent console
  overflow

- POD5 files must correspond exactly to reads in the Dorado summary

- For datasets \>1M reads, consider batch processing or increased
  `part_size`

- Position estimates are approximate (Â±2-3 nt accuracy); validate
  critical findings

## Pipeline Overview

The analysis pipeline consists of several integrated stages:

1.  **Input Validation**: Validates Dorado summary files and POD5
    directories

2.  **Data Preprocessing**: Splits large datasets and applies quality
    filters

3.  **Signal Extraction**: Extracts poly(A) tail signals from POD5 files

4.  **Feature Engineering**: Computes pseudomoves and signal
    characteristics

5.  **Segmentation**: Identifies candidate modification regions

6.  **GAF Creation**: Creates Gramian Angular Fields for CNN input

7.  **Classification**: Applies trained neural networks for prediction

8.  **Output Creation**: Produces comprehensive results and statistics

9.  **Cleanup**: Optionally removes intermediate files to save disk
    space

## Input Requirements

This pipeline requires specific input formats:

- **Dorado Summary**: Must contain poly(A) information columns:
  `read_id`, `filename`, `poly_tail_length`, `poly_tail_start`,
  `poly_tail_end`

- **POD5 Files**: Raw signal files corresponding to reads in summary

## Output table explanation

**1. read_classes Table** This table provides a complete accounting of
ALL reads in the analysis, with each read assigned to one of three
biological categories (class) and given a specific technical code
(comments) explaining the classification.

**2. nonadenosine_residues Table**

This table provides modification-level detail for decorated reads only.
Each row represents a single predicted non-adenosine residue within a
poly(A) tail. Reads can have multiple rows if multiple modifications
detected.

## Output Structure

The function creates several output subdirectories:

- `dorado_summary_dir/`:

  Split summary files for parallel processing

- `polya_signal_dir/`:

  Extracted poly(A) signals in RDS format

- `nonA_temp_dir/`:

  Intermediate non-A prediction files

- `polya_chunks_dir/`:

  Signal segmentation data

## Performance Characteristics

- **Scalability**: Efficient parallel processing across multiple cores

- **Memory Management**: Smart data partitioning prevents memory
  overflow

- **Progress Tracking**: Comprehensive logging and progress indicators

- **Error Handling**: Robust error recovery and informative diagnostics

## Quality Control

When `qc = TRUE`, the pipeline applies several quality filters:

- Poly(A) tail length filtering (minimum 10 nucleotides)

- Coordinate validation (start \< end positions)

- Terminal position masking to reduce false positives

- Duplicate read detection and handling

## Algorithm Details

The pipeline employs several sophisticated algorithms:

- **Pseudomove Computation**: Z-score based outlier detection with
  adaptive windowing

- **Signal Segmentation**: Run-length encoding for modification region
  identification

- **GAF Transformation**: Gramian Angular Field generation for CNN
  compatibility

- **Neural Classification**: Pre-trained CNNs for nucleotide type
  prediction

## Classification Codes

The `comments` column uses standardized 3-letter codes for precise
technical documentation. These codes are essential for filtering
datasets and understanding pipeline behavior:

**Comments:**

- YAY:

  Non-A residue successfully detected and classified

- MAU:

  No signal deviations detected; pure poly(A) signal

- MPU:

  Signal deviations present but predicted as adenosine only

- IRL:

  Poly(A) tail too short (\< 10 nt) for reliable analysis

- UNM:

  Read unmapped to reference genome

- BAC:

  Invalid coordinate system (poly_tail_start = 0)

## System Requirements

- **R Version**: \>= 3.5.0

- **Python**: \>= 3.6 with pod5 module (`pip install pod5`)

- **Memory**: \>= 8GB RAM (16GB+ recommended for large datasets)

- **Storage**: Temporary space ~2-5x input file size

- **Dependencies**: Tensorflow/Keras for neural network inference

## Error Handling

- Input validation with informative error messages

- Graceful handling of corrupted or missing files

- Automatic retry mechanisms for transient failures

- Detailed logging of all errors and warnings

- Safe cleanup of temporary files on failure

## Performance Tips

- Use SSDs for `save_dir` to improve I/O performance

- Adjust `part_size` based on available RAM (larger = faster, more
  memory)

- Use `num_cores = parallel::detectCores() - 1` for maximum speed

- Ensure POD5 files and output directory are on fast storage

- Consider preprocessing very large datasets (\>1M reads) in batches

## See also

[`check_tails_guppy`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
for legacy Guppy-based analysis,
[`preprocess_inputs`](https://LRB-IIMCB.github.io/ninetails/reference/preprocess_inputs.md)
for input preprocessing,
[`process_dorado_signal_files`](https://LRB-IIMCB.github.io/ninetails/reference/process_dorado_signal_files.md)
for signal processing,
[`create_outputs_dorado`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado.md)
for output generation,
[`filter_signal_by_threshold`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)
for signal analysis

## Examples

``` r
if (FALSE) { # \dontrun{
results <- check_tails_dorado_DRS(
  dorado_summary = "path/to/alignment_summary.txt",
  pod5_dir = "path/to/pod5/",
  num_cores = 2,
  qc = TRUE,
  save_dir = "path/to/output/",
  prefix = "experiment1",
  part_size = 40000,
  cleanup = TRUE
)
} # }
```
