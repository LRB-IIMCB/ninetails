# Detection of non-adenosines

This vignette describes the three analysis pipelines available in
**ninetails** for detecting non-adenosine residues in poly(A) tails. For
working with the output data, see
[`vignette("postprocessing")`](https://LRB-IIMCB.github.io/ninetails/articles/postprocessing.md).
For visualization, see
[`vignette("plotting")`](https://LRB-IIMCB.github.io/ninetails/articles/plotting.md)
and
[`vignette("shiny_app")`](https://LRB-IIMCB.github.io/ninetails/articles/shiny_app.md).

## Dorado DRS Pipeline (Recommended)

[`check_tails_dorado_DRS()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_DRS.md)
is the recommended function for direct RNA sequencing (DRS) data
basecalled with Dorado ≥ 1.0.0 using POD5 format.

``` r

results <- ninetails::check_tails_dorado_DRS(
  dorado_summary = "path/to/dorado_alignment_summary.txt",
  pod5_dir       = "path/to/pod5_dir/",
  num_cores      = 2,
  qc             = TRUE,
  save_dir       = "~/output/",
  prefix         = "experiment1",
  part_size      = 40000,
  cleanup        = FALSE
)
```

### Parameters

| Parameter | Type | Default | Description |
|----|----|----|----|
| `dorado_summary` | character/data.frame | *required* | Path to Dorado summary file. Must contain: `read_id`, `filename`, `poly_tail_length`, `poly_tail_start`, `poly_tail_end` |
| `pod5_dir` | character | *required* | Path to directory containing POD5 files |
| `num_cores` | integer | `1` | Number of CPU cores for parallel processing |
| `qc` | logical | `TRUE` | Apply quality control filtering (see below) |
| `save_dir` | character | *required* | Output directory for results files |
| `prefix` | character | `""` | Optional prefix for output filenames (e.g., `"experiment1_"`) |
| `part_size` | integer | `40000` | Maximum reads per processing chunk (reduce if memory limited) |
| `cleanup` | logical | `FALSE` | Remove intermediate files after completion |

### Quality control (`qc = TRUE`)

When `qc = TRUE`, the pipeline applies the following filters before
signal analysis:

- Reads with `poly_tail_start = 0` are classified as `BAC` (bad
  coordinates) and excluded from signal processing
- Reads with `poly_tail_length < 10` are classified as `IRL`
  (insufficient read length)
- Unmapped reads (no `contig` assignment) are classified as `UNM`

These reads are still included in the output `read_classes` table but
are marked as `unclassified`.

### Generating the Dorado summary file

The Dorado summary file is produced by running `dorado summary` on your
aligned BAM file:

``` bash
dorado summary aligned_reads.bam > dorado_summary.txt
```

The summary must include poly(A) tail estimates. Ensure that Dorado was
run with `--estimate-poly-a` during basecalling.

### Pipeline steps

The pipeline executes the following steps in order:

1.  **Input validation**
    ([`preprocess_inputs()`](https://LRB-IIMCB.github.io/ninetails/reference/preprocess_inputs.md)):
    Validates file paths, checks required columns, verifies POD5
    directory structure
2.  **Summary processing**
    ([`process_dorado_summary()`](https://LRB-IIMCB.github.io/ninetails/reference/process_dorado_summary.md)):
    Reads and standardizes the Dorado summary file
3.  **Quality filtering**
    ([`filter_dorado_summary()`](https://LRB-IIMCB.github.io/ninetails/reference/filter_dorado_summary.md)):
    Applies QC filters and assigns preliminary classifications
4.  **Signal extraction**
    ([`extract_tails_from_pod5()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tails_from_pod5.md)):
    Extracts raw poly(A) tail signals from POD5 files via the Python
    `pod5` module
5.  **Feature engineering**
    ([`create_tail_features_list_dorado()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_features_list_dorado.md)):
    Computes pseudomoves, signal statistics, and identifies candidate
    modification positions
6.  **Segmentation**
    ([`create_tail_chunk_list_dorado()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_dorado.md)):
    Splits tail signals into 100-point segments centered on candidate
    positions
7.  **GAF creation**
    ([`process_dorado_signal_files()`](https://LRB-IIMCB.github.io/ninetails/reference/process_dorado_signal_files.md)):
    Converts segments to Gramian Angular Fields (GASF + GADF two-channel
    images)
8.  **CNN classification**
    ([`predict_gaf_classes()`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md)):
    Applies the trained convolutional neural network to classify each
    segment as A, C, G, or U
9.  **Output creation**
    ([`create_outputs_dorado()`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado.md)):
    Aggregates per-segment predictions into per-read classifications and
    produces final output tables

### Output files

Two tab-separated files are saved in `save_dir`:

- `{prefix}read_classes.txt` - one row per read
- `{prefix}nonadenosine_residues.txt` - one row per detected non-A
  residue (decorated reads only)

The function also returns both tables as a named list:

``` r

class_data <- results$read_classes
residue_data <- results$nonadenosine_residues
```

------------------------------------------------------------------------

## Dorado cDNA Pipeline (Experimental)

> **Warning**
>
> The cDNA pipeline is experimental. The models were trained and tested
> on synthetic spike-ins only (see *Validation* below), so results on
> biological samples should be interpreted with caution.

[`check_tails_dorado_cDNA()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_cDNA.md)
detects non-adenosine residues in ONT cDNA reads. A cDNA read carries
the tail either as poly(A) at its 3’ end or as poly(T) at its 5’ end,
depending on the strand that was sequenced. Each orientation is
classified with a dedicated model, and the predictions are always
reported as residues of the original RNA tail (C, G, or U).

Tail coordinates can be provided in two ways. We recommend the Lachesis
output, because the models were trained on Lachesis coordinates.

### Lachesis input (recommended)

The Lachesis table provides the tail coordinates, the tail length, the
QC tag, and the read orientation (`strand`), so no BAM file is needed.

``` r

test_dir <- system.file('extdata', 'test_data', 'pod5_cDNA',
                        package = 'ninetails')

results <- ninetails::check_tails_dorado_cDNA(
  lachesis  = file.path(test_dir, 'lachesis_polya.tsv'),
  pod5_dir  = test_dir,
  num_cores = 2,
  qc        = TRUE,
  save_dir  = '~/output/',
  prefix    = 'experiment1',
  part_size = 40000,
  cleanup   = FALSE
)
```

### Dorado input

Without Lachesis, the coordinates are taken from the Dorado summary and
the orientation is determined from the basecalled sequences in the BAM
file (SSP/VNP primer matching).

``` r

results <- ninetails::check_tails_dorado_cDNA(
  bam_file       = file.path(test_dir, 'example_cDNA_dorado.bam'),
  dorado_summary = file.path(test_dir, 'dorado_summary.txt'),
  pod5_dir       = test_dir,
  num_cores      = 2,
  qc             = TRUE,
  save_dir       = '~/output/',
  prefix         = 'experiment1',
  part_size      = 40000,
  cleanup        = FALSE
)
```

Both examples run on the test data shipped with the package (26 reads of
the beta-gal spike-ins, 13 per orientation). The expected output tables
are stored next to the data as `test_cDNA_*.tsv` (Lachesis input) and
`test_cDNA_dorado_*.tsv` (Dorado input).

### Parameters

| Parameter | Type | Default | Description |
|----|----|----|----|
| `lachesis` | character/data.frame | `NULL` | Lachesis poly(A) output. Must contain: `readname`, `contig`, `adapter_end`, `polya_end`, `polya_length`, `qc_tag`, `strand`. If provided, `bam_file` and `dorado_summary` are not used |
| `bam_file` | character | *required without* `lachesis` | Path to BAM file with basecalled cDNA reads |
| `dorado_summary` | character/data.frame | *required without* `lachesis` | Path to Dorado summary file with poly(A) tail estimates |
| `pod5_dir` | character | *required* | Path to directory containing POD5 files |
| `num_cores` | integer | `1` | Number of CPU cores for parallel processing |
| `qc` | logical | `TRUE` | Apply quality control filtering (see below) |
| `save_dir` | character | *required* | Output directory for results files |
| `prefix` | character | `""` | Optional prefix for output filenames |
| `part_size` | integer | `40000` | Maximum reads per processing chunk (reduce if memory limited) |
| `cleanup` | logical | `FALSE` | Remove intermediate files after completion |

### Quality control (`qc = TRUE`)

- Reads with a Lachesis QC tag other than `PASS` are classified as `QCF`
  and excluded from signal processing
- Reads with a tail shorter than 10 nt are classified as `IRL`
- Reads with invalid tail coordinates, or without a signal, are
  classified as `BAC`
- Calls closer than 2 nt to either end of the tail are discarded; a read
  left without any call is reported as `blank` with code `MPU`

### Pipeline steps

1.  **Tail coordinates**
    ([`read_tail_coordinates_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/read_tail_coordinates_cdna.md)):
    Converts the Lachesis output or the Dorado summary to a common
    table. With Dorado input, the BAM file is split
    ([`split_bam_file_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/split_bam_file_cdna.md)),
    sequences are extracted
    ([`extract_data_from_bam()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_data_from_bam.md)),
    and the orientation is detected
    ([`detect_orientation_multiple()`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_multiple.md))
2.  **Signal extraction**
    ([`extract_tail_signals_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_signals_cdna.md)):
    Extracts the raw tail signal with 100 samples of flank on each side.
    Unlike the DRS pipeline, the signal is neither winsorised nor
    downsampled, because a nucleotide spans only about 11 samples in
    cDNA
3.  **Normalisation**
    ([`normalise_tail_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/normalise_tail_cdna.md)):
    Centres the signal on the tail median, scales it by a noise
    estimate, and compresses it to the interval (-1, 1)
4.  **Sliding windows**
    ([`create_tail_windows_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_windows_cdna.md)):
    Cuts 100-sample windows every 25 samples along the whole tail. No
    candidate positions are preselected
5.  **GAF creation**
    ([`create_gaf_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_cdna.md)):
    Converts each window to a two-channel Gramian Angular Field
6.  **CNN classification**
    ([`predict_nonA_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/predict_nonA_cdna.md)):
    Applies the poly(A) or the poly(T) model to every window
7.  **Calling**
    ([`call_nonA_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/call_nonA_cdna.md)):
    Reports a non-adenosine where at least two consecutive windows
    exceed the threshold; stricter thresholds apply near the tail
    boundaries of poly(A) reads
8.  **Output creation**
    ([`create_outputs_dorado_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado_cdna.md),
    [`merge_cdna_results()`](https://LRB-IIMCB.github.io/ninetails/reference/merge_cdna_results.md)):
    Aggregates the calls into per-read classifications for both
    orientations

As in the DRS pipeline, reads are processed in parts of `part_size`
reads and the intermediate results of each part are written to
`save_dir`, so memory usage does not grow with the size of the run.

### Output

The output tables have the same structure as in the DRS pipeline (see
*Output structure* below), with an additional `tail_type` column
(`polyA` or `polyT`) indicating the read orientation.

``` r

class_data <- results$read_classes
residue_data <- results$nonadenosine_residues

head(residue_data, 3)
#>                               readname   contig prediction est_nonA_pos polya_length qc_tag tail_type
#> 1 00520674-3559-4dd2-a7dc-d5449c15806c beta_gal          G        31.72        69.94   PASS     polyA
#> 2 144f3bab-764b-492f-a235-f9787e5ce905 beta_gal          U        35.26        64.64   PASS     polyA
#> 3 16624e6e-c846-4f5f-b28e-26e6646d15b9 beta_gal          C        33.61        66.25   PASS     polyA
```

`est_nonA_pos` is computed from the position of the call in the signal
and the provided tail length, so it inherits any error of the tail
length estimate. Positions are resolved in steps of about 2 nt.

### Validation

The models were trained on in vitro transcribed beta-gal spike-ins
carrying a single C, G, or U in the middle of a 60 nt poly(A) tail
(SQK-PCB114-24, 5 kHz, Dorado hac v6.0, Lachesis coordinates). On reads
from sequencing channels not used for training:

| Orientation | Recall | Residue accuracy | False-positive reads on pure poly(A) tails |
|----|----|----|----|
| poly(A) | 0.92 | 0.90 | 0.039 |
| poly(T) | 0.86 | 0.80 | 0.023 |

Keep in mind the limits of this validation:

- Only one tail length (60 nt) and one non-A position (the middle of the
  tail) were available for training
- False-positive rates were measured on a pure poly(A) control with
  tails of 10–200 nt; for poly(T) reads, tails longer than about 65 nt
  could not be assessed
- Calls close to the ends of tails shorter than 45 nt are less reliable
- Processing is considerably slower than in the DRS pipeline, because
  every window of every tail is classified

------------------------------------------------------------------------

## Guppy Legacy Pipeline

> **Warning**
>
> This pipeline works with Guppy basecaller 6.0.0 and lower only. It is
> maintained for backward compatibility and will not be further
> optimised. For new analyses, use the Dorado DRS pipeline.

[`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
is the legacy function for DRS data basecalled with Guppy using fast5
format and Nanopolish poly(A) coordinates.

``` r

results <- ninetails::check_tails_guppy(
  polya_data = system.file('extdata', 'test_data',
                           'nanopolish_output.tsv',
                           package = 'ninetails'),
  sequencing_summary = system.file('extdata', 'test_data',
                                   'sequencing_summary.txt',
                                   package = 'ninetails'),
  workspace = system.file('extdata', 'test_data',
                          'basecalled_fast5',
                          package = 'ninetails'),
  num_cores = 2,
  basecall_group = 'Basecall_1D_000',
  pass_only = TRUE,
  save_dir = '~/output/'
)
```

### Parameters

| Parameter | Type | Default | Description |
|----|----|----|----|
| `polya_data` | character/data.frame | *required* | Path to Nanopolish polya output or tailfindr output |
| `sequencing_summary` | character/data.frame | *required* | Path to sequencing summary file |
| `workspace` | character | *required* | Path to directory with multi-fast5 files |
| `num_cores` | integer | `1` | Number of CPU cores |
| `basecall_group` | character | `"Basecall_1D_000"` | Fast5 hierarchy level |
| `pass_only` | logical | `TRUE` | Include only “PASS” reads from Nanopolish QC |
| `qc` | logical | `TRUE` | Label terminal positions with “-WARN” suffix |
| `save_dir` | character | *required* | Output directory |
| `part_size` | integer | `1000000` | Max rows per processing chunk |

> **Note:** For tailfindR compatibility, see
> [`vignette("tailfindr")`](https://LRB-IIMCB.github.io/ninetails/articles/tailfindr.md).

------------------------------------------------------------------------

## Output structure

All pipelines return a named list with two data frames:

### read_classes

Complete accounting of all reads in the analysis.

| Column         | Description                                             |
|----------------|---------------------------------------------------------|
| `readname`     | Unique read identifier                                  |
| `contig`       | Reference sequence/transcript                           |
| `polya_length` | Estimated poly(A) tail length (nt)                      |
| `qc_tag`       | Quality tag from basecaller QC                          |
| `class`        | Classification: `decorated`, `blank`, or `unclassified` |
| `comments`     | 3-letter code explaining the classification (see below) |
| `tail_type`    | `polyA` or `polyT` (cDNA pipeline only)                 |

### nonadenosine_residues

Modification-level detail for `decorated` reads only. Each row
represents one detected non-adenosine residue.

| Column | Description |
|----|----|
| `readname` | Unique read identifier |
| `contig` | Reference sequence/transcript |
| `prediction` | Predicted nucleotide: `C`, `G`, or `U` |
| `est_nonA_pos` | Estimated position within the poly(A) tail (nucleotides from the 3’ end) |
| `polya_length` | Total tail length (nt) |
| `qc_tag` | Quality tag |
| `tail_type` | `polyA` or `polyT` (cDNA pipeline only) |

> **Note:** A single read can have multiple rows in
> `nonadenosine_residues` if it contains more than one non-A residue.

------------------------------------------------------------------------

## Classification codes

### Dorado pipelines (DRS and cDNA)

| Class | Code | Explanation |
|----|----|----|
| `decorated` | `YAY` | Non-adenosine residue detected |
| `blank` | `MAU` | No signal deviation; pure poly(A) tail |
| `blank` | `MPU` | Signal deviation present but predicted as adenosine |
| `unclassified` | `IRL` | Poly(A) tail too short (\< 10 nt) |
| `unclassified` | `UNM` | Read unmapped |
| `unclassified` | `BAC` | Invalid coordinates (`poly_tail_start = 0`) |
| `unclassified` | `QCF` | Lachesis QC failed (cDNA pipeline with Lachesis input only) |

In the cDNA pipeline, `MPU` marks reads whose only calls were terminal
and discarded by `qc`.

### Guppy legacy pipeline

| Class          | Code  | Explanation                               |
|----------------|-------|-------------------------------------------|
| `decorated`    | `YAY` | Move transition present, non-A detected   |
| `blank`        | `MAU` | Move transition absent, non-A undetected  |
| `blank`        | `MPU` | Move transition present, non-A undetected |
| `unclassified` | `QCF` | Nanopolish QC failed                      |
| `unclassified` | `NIN` | Not included (`pass_only = TRUE`)         |
| `unclassified` | `IRL` | Insufficient read length                  |

------------------------------------------------------------------------

## Memory considerations

> **Warning**
>
> Signal transformations can place a heavy load on memory, especially
> for full sequencing runs with millions of reads.

To manage memory:

- Adjust `part_size` to control reads per chunk. Default is 40,000
  (Dorado) or 1,000,000 (Guppy). Reduce to 10,000–20,000 if experiencing
  memory issues.
- Use `cleanup = TRUE` to remove intermediate files after processing
  completes.
- Monitor system memory during the GAF creation and CNN prediction
  steps, which are the most memory-intensive.
- For very large runs (\>1M reads), consider splitting input files
  manually and processing in batches.

------------------------------------------------------------------------

## Performance tips

- **Parallelization**: Set `num_cores` to the number of available CPU
  cores minus one. The pipeline parallelizes signal extraction, feature
  engineering, and GAF creation across chunks.
- **Disk space**: Intermediate files (signal data, feature lists, GAF
  matrices) can be large. Ensure sufficient free disk space in
  `save_dir`. Use `cleanup = TRUE` to remove intermediates
  automatically.
- **GPU acceleration**: If TensorFlow is configured with GPU support,
  the CNN classification step benefits from GPU acceleration. This is
  the single largest speedup available.
- **SSD storage**: Using SSD storage for `pod5_dir` and `save_dir`
  significantly reduces I/O bottlenecks during signal extraction.

------------------------------------------------------------------------

## Next steps

After running the pipeline:

1.  **Postprocess**: Load, annotate, merge, and summarize results -
    [`vignette("postprocessing")`](https://LRB-IIMCB.github.io/ninetails/articles/postprocessing.md)
2.  **Visualize**: Generate static plots -
    [`vignette("plotting")`](https://LRB-IIMCB.github.io/ninetails/articles/plotting.md)
3.  **Inspect signals**: Examine raw squiggles -
    [`vignette("signal_inspection")`](https://LRB-IIMCB.github.io/ninetails/articles/signal_inspection.md)
4.  **Explore interactively**: Launch the Shiny dashboard -
    [`vignette("shiny_app")`](https://LRB-IIMCB.github.io/ninetails/articles/shiny_app.md)
