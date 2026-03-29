# Get started with ninetails

**Ninetails** detects and characterises non-adenosine nucleotides
embedded within poly(A) tails of Oxford Nanopore sequencing reads. This
vignette provides a quick introduction to installation and running the
main analysis pipelines.

For complete documentation, see the [Ninetails
Wiki](https://github.com/LRB-IIMCB/ninetails/wiki).

## Installation

**Ninetails** is not currently available on CRAN or Bioconductor.
Install it from GitHub using `devtools`:

``` r
install.packages("devtools")
devtools::install_github('LRB-IIMCB/ninetails')
library(ninetails)
```

> **Note for Windows users:** Before installing `devtools` on Windows,
> install Rtools so packages compile correctly:
> <https://cran.r-project.org/bin/windows/Rtools/>

### Additional dependencies

Depending on which pipeline you use, additional components are required:

**For Dorado pipelines (DRS and cDNA):**

- Python with `pod5` module: `pip install pod5`
- Keras/TensorFlow for R:
  [`keras::install_keras()`](https://rdrr.io/pkg/keras/man/install_keras.html)

**For Guppy legacy pipeline:**

- `rhdf5` from Bioconductor: `BiocManager::install("rhdf5")`
- Keras/TensorFlow for R
- VBZ compression plugin (for newer MinKNOW data)

See the
[Wiki](https://github.com/LRB-IIMCB/ninetails/wiki/3.-Additional-requirements)
for detailed instructions.

## Available pipelines

Ninetails provides three analysis pipelines:

| Pipeline        | Basecaller     | Input format         | Function                                                                                                  | Status            |
|-----------------|----------------|----------------------|-----------------------------------------------------------------------------------------------------------|-------------------|
| **Dorado DRS**  | Dorado ≥ 1.0.0 | POD5 + summary       | [`check_tails_dorado_DRS()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_DRS.md)   | Recommended       |
| **Dorado cDNA** | Dorado ≥ 1.0.0 | POD5 + BAM + summary | [`check_tails_dorado_cDNA()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_cDNA.md) | Under development |
| **Guppy**       | Guppy ≤ 6.0.0  | fast5 + Nanopolish   | [`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)             | Legacy            |

## Quick start: Dorado DRS pipeline

The recommended pipeline for direct RNA sequencing (DRS) data:

``` r
results <- ninetails::check_tails_dorado_DRS(
  dorado_summary = "path/to/dorado_summary.txt",
  pod5_dir       = "path/to/pod5_dir/",
  num_cores      = 2,
  qc             = TRUE,
  save_dir       = "~/output/",
  prefix         = "my_experiment"
)
```

### Required inputs

- **dorado_summary**: Dorado summary file with columns `read_id`,
  `filename`, `poly_tail_length`, `poly_tail_start`, `poly_tail_end`
- **pod5_dir**: Directory containing POD5 files

### Output

The function returns a list with two data frames:

- **read_classes**: Per-read classification
  (decorated/blank/unclassified)
- **nonadenosine_residues**: Detailed non-A positions for decorated
  reads

Results are also saved as TSV files in `save_dir`.

## Next steps

- [`vignette("detection")`](https://LRB-IIMCB.github.io/ninetails/articles/detection.md)
  — Detailed pipeline documentation
- [`vignette("postprocessing")`](https://LRB-IIMCB.github.io/ninetails/articles/postprocessing.md)
  — Working with results
- [`vignette("plotting")`](https://LRB-IIMCB.github.io/ninetails/articles/plotting.md)
  — Visualization functions
- [Ninetails Wiki](https://github.com/LRB-IIMCB/ninetails/wiki) —
  Complete documentation

## Citation

Please cite Ninetails as:

Gumińska, N., Matylla-Kulińska, K., Krawczyk, P.S. et al. Direct
profiling of non-adenosines in poly(A) tails of endogenous and
therapeutic mRNAs with Ninetails. *Nat Commun* **16**, 2664 (2025).
<https://doi.org/10.1038/s41467-025-57787-6>
