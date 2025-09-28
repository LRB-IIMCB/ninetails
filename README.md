
<a href="https://github.com/nemitheasura/https://github.com/LRB-IIMCB/ninetails/releases/"><img src="https://img.shields.io/github/tag/LRB-IIMCB/ninetails?include_prereleases=&sort=semver&color=blue" alt="GitHub tag"></a>
<a href="https://github.com/nemitheasura/https://github.com/LRB-IIMCB/ninetails/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue" alt="License"></a>
<img src="https://img.shields.io/badge/maintained-yes-blue" alt="maintained - yes">
<a href="https://www.linux.org/" title="Go to Linux homepage"><img src="https://img.shields.io/badge/OS-Linux-blue?logo=linux&logoColor=white" alt="OS - Linux"></a>
<a href="https://www.microsoft.com/" title="Go to Microsoft homepage"><img src="https://img.shields.io/badge/OS-Windows-blue?logo=windows&logoColor=white" alt="OS - Windows"></a>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13309819.svg)](https://doi.org/10.5281/zenodo.13309819)

# Ninetails

**An R package for finding non-adenosine poly(A) residues in Oxford Nanopore direct RNA sequencing reads**

<img src="https://user-images.githubusercontent.com/68285258/168554098-a5a5dee9-2c8f-4351-86b4-e6420a5b8ced.png" align="right" width="200" height="220"/>

## Introduction

-   It works on Oxford Nanopore direct RNA sequencing reads basecalled by Guppy software
-   It requires tail delimitation data produced by Nanopolish software
-   It allows both for the detection of non-adenosine residues within the poly(A) tails and visual inspection of read signals

**Important notes!** 

-   **Since version 1.0.4 Ninetails is compatible with Dorado basecaller & pod5 format**
-   **Since version 1.0.2 Ninetails is compatible also with tailfindR**


Currently, **Ninetails** can distinguish characteristic signatures of four types of nucleotides: adenosines (A), cytosines (C), guanosines (G), and uridines (U).
<div>


![schemat_sieci_presentation_smol](https://github.com/LRB-IIMCB/ninetails/assets/68285258/5e110a7a-b00f-4e64-903f-1e5c11db6172)


> **Note**
> 
> **For detailed documentation including explanation of additional dataprocessing and datavis features see <a href="https://github.com/LRB-IIMCB/ninetails/wiki">Ninetails' Wiki</a>**
>
</div>


The software is still under development, so all suggestions to improving it are welcome. Please note that the code contained herein may change frequently, so use it with caution.

**Ninetails** was tested on Linux Mint 20.3, Ubuntu 20.04.3 and Windows 11 operating systems with R 4.1.2, R 4.2.0 and R 4.2.1.


## Installation

Currently, **Ninetails** is not available on CRAN/Bioconductor, so you need to install it using `devtools`.

If you do not have `devtools` installed already, you can do this with the following command in R/RStudio:

``` r
install.packages("devtools")
```

<div>

> **Note**
>
> **For Windows users:**
> 
> Before installation of `devtools` on Windows, you should install `Rtools`, so the packages would be correctly compiled: <https://cran.r-project.org/bin/windows/Rtools/>

</div>

Once you have `devtools` installed, you can install **Ninetails** using the command below in R/RStudio:

``` r
devtools::install_github('LRB-IIMCB/ninetails')
library(ninetails)
```
The installation of the repo takes approx. 20 seconds on typical PC. Additional time is required to install and configure additional components.


**Important info: Ninetails requires additional components/third party tools to operate.**
**For further info, read [Wiki](https://github.com/LRB-IIMCB/ninetails/wiki)**


## General usage

### Classification of reads using wrapper function

The pipeline wrappers available in **Ninetails** begin with `check_tails` prefix. 

In **Ninetails 1.0.4+**, `check_tails_dorado_DRS()` is the pipeline wrapper designed to classify reads/non-adenosines in ONT DRS data. 
``` r

results <- ninetails::check_tails_dorado_DRS(
   dorado_summary = "path/to/dorado_alignment_summary.txt",
   pod5_dir = "path/to/pod5_dir/",
   num_cores = 2,
   qc = TRUE,
   save_dir = "~/Downloads",
   prefix = "experiment1")
```

In **Ninetails <1.0.4** , `check_tails()` is the main function which allows to classify sequencing reads based on presence/absence of non-adenosine residues within their poly(A) tails (and additional conditions, such as minimal read length and qc_tag assigned by Nanopolish polya function). In newer versions (1.0.4+) this was replaced by `check_tails_guppy()`. This is a legacy pipeline which wouldnot be further developed (with exception of critical bug fixes).

Below is an example of how to use `check_tails_guppy()` function:

``` r
results <- ninetails::check_tails_guppy(
  nanopolish = system.file('extdata', 
                           'test_data', 
                           'nanopolish_output.tsv', 
                           package = 'ninetails'),
  sequencing_summary = system.file('extdata', 
                                   'test_data', 
                                   'sequencing_summary.txt', 
                                   package = 'ninetails'),
  workspace = system.file('extdata', 
                          'test_data', 
                          'basecalled_fast5', 
                          package = 'ninetails'),
  num_cores = 2,
  basecall_group = 'Basecall_1D_000',
  pass_only=TRUE,
  save_dir = '~/Downloads')
```

This function returns a list consisting of two tables: **read_classes** and **nonadenosine_residues**. In addition, the function saves results to text files in the user-specified directory.

Moreover, the function also creates a log file in the directory specified by the user.

The runtime depends on the hardware resources and sequencing depth. The processing of built-in test dataset should take around 1 minute.


### Output explanation

#### The **read_classes** dataframe (file) contains following columns:

| column name  | content                                                                  |
|------------------------------------|------------------------------------|
| readname     | an identifier of a given read (36 characters)                            |
| contig       | reference to which the given read was mapped (inherited from nanopolish) |
| polya_length | tail length estimation provided by nanopolish polya function             |
| qc_tag       | quality tag assigned by nanopolish polya function                        |
| class        | the crude result of classification                                       |
| comments     | a code indicating whether the classification criteria were met/unmet     |

The `class` column contains information whether the given read was recognized as modified (containing non-adenosine residue) or not. Whereas the `comment` column contains details underlying the classification outcome. The content of these columns is explained below:

| class        | comments | explanation                                      |
|--------------|----------|--------------------------------------------------|
| modified     | YAY      | move transition present, nonA residue detected   |
| unmodified   | MAU      | move transition absent, nonA residue undetected  |
| unmodified   | MPU      | move transition present, nonA residue undetected |
| unclassified | QCF      | nanopolish qc failed                             |
| unclassified | NIN      | not included in the analysis (pass only = T)     |
| unclassified | IRL      | insufficient read length                         |

#### The **nonadenosine_residues** dataframe (file) contains following columns:

| column name  | content                                                                                        |
|------------------------------------|------------------------------------|
| readname     | an identifier of a given read (36 characters)                                                  |
| prediction   | the result of classification (basic model: C, G, U assignment)                                 |
| est_nonA_pos | the approximate nucleotide position where nonadenosine is to be expected; reported from 5' end |
| polya_length | the tail length estimated according to Nanopolish polya function                               |
| qc_tag       | quality tag assigned by nanopolish polya function                                              |

## Important notes

<div>

> **Warning**
>
> **Current pre-release versions of the package work with Guppy basecaller 6.0.0 and lower. Please be aware to use compatible version of basecaller.**
>
</div>

Before running the program, it is recommended to ascertain that the given arguments (nanopolish, sequencing summary and directory with fast5 files) correspond with each other. In other words, that the records contained in the `Nanopolish` polyA output file correspond to the records contained in the sequencing summary file and in the fast5 files stored in the declared directory (workspace). If a complete discrepancy is detected, the program will not perform the analysis. Instead, it will throw an error. In case of the presence of incompatible records - they will be omitted from the result files and the pipeline will end with warning.

<div>

> **Warning**
> 
> Please be aware that signal transformations performed during analysis can place a heavy load on memory. This is especially true if your data covers the entire sequencing run. 
> 
</div>

<div>

> **Note**
>
> Currently, **Ninetails** does not support single fast5 files as this format is deprecated by ONT. Before running the program on single fast5 files, you should convert them to multifast5 with another tool, for instance with `ont-fast5-api`.
> 
</div>

<div>

> **Note**
> 
>**Ninetails** relies on external signal segmentation (Nanopolish/Tailfindr/Dorado) and therefore may underestimate terminal modifications (last and penultimate nucleotides of the tail).
>
</div>

## Citation

Please cite **Ninetails** as:

Gumińska, N., Matylla-Kulińska, K., Krawczyk, P.S. et al. Direct profiling of non-adenosines in poly(A) tails of endogenous and therapeutic mRNAs with Ninetails. Nat Commun 16, 2664 (2025). https://doi.org/10.1038/s41467-025-57787-6


## Troubleshooting

If you encounter a bug, please post it on github. To help diagnose the problem, send a minimal reproducible example (required inputs covering around 5-10 reads + corresponding nanopolish output & sequencing summary), so I will be able to reproduce the error and fix it for you.

## Maintainer

Any issues regarding the **Ninetails** should be addressed to Natalia Gumińska (nguminska (at) iimcb.gov.pl).

**Ninetails** has beed developed in the <a href="https://www.iimcb.gov.pl/en/research/41-laboratory-of-rna-biology-era-chairs-group">Laboratory of RNA Biology</a> (Dziembowski Lab) at the <a href="https://www.iimcb.gov.pl/en/">International Institute of Molecular and Cell Biology</a> in Warsaw. 
