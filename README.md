
<a href="https://github.com/nemitheasura/https://github.com/LRB-IIMCB/ninetails/releases/"><img src="https://img.shields.io/github/tag/LRB-IIMCB/ninetails?include_prereleases=&sort=semver&color=blue" alt="GitHub tag"></a>
<a href="https://github.com/nemitheasura/https://github.com/LRB-IIMCB/ninetails/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue" alt="License"></a>
<img src="https://img.shields.io/badge/maintained-yes-blue" alt="maintained - yes">
<a href="https://www.linux.org/" title="Go to Linux homepage"><img src="https://img.shields.io/badge/OS-Linux-blue?logo=linux&logoColor=white" alt="OS - Linux"></a>
<a href="https://www.microsoft.com/" title="Go to Microsoft homepage"><img src="https://img.shields.io/badge/OS-Windows-blue?logo=windows&logoColor=white" alt="OS - Windows"></a>


# Ninetails

**An R package for finding non-adenosine poly(A) residues in Oxford Nanopore direct RNA sequencing reads**

<img src="https://user-images.githubusercontent.com/68285258/168554098-a5a5dee9-2c8f-4351-86b4-e6420a5b8ced.png" align="right" width="200" height="220"/>

## Introduction

-   It works on Oxford Nanopore direct RNA sequencing reads basecalled by Guppy software
-   It requires tail delimitation data produced by Nanopolish software
-   It allows both for the detection of non-adenosine residues within the poly(A) tails and visual inspection of read signals

Currently, **Ninetails** can distinguish characteristic signatures of four types of nucleotides: adenosines (A), cytosines (C), guanosines (G), and uridines (U).

**Ninetails** relies on Nanopolish segmentation and therefore may underestimate terminal modifications (last and penultimate nucleotides of the tail).

<div>

> **Note**
>
> The software is still under development, so all suggestions to improving it are welcome. Please note that the code contained herein may change frequently, so use it with caution.
>
</div>


**Ninetails** was tested on Linux Mint 20.3, Ubuntu 20.04.3 and Windows 11 operating systems with R 4.1.2, R 4.2.0 and R 4.2.1.

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

For the moment, **Ninetails** does not offer the possibility of processing large data sets in chunks behind the scenes (under development). Therefore, to minimise the risk of unexpected crashes, it is highly recommended to split the output of the `Nanopolish` polyA function into smaller files to make it easier to process the data in subsets and then merge the final results.

## Prerequisites

### Input data & preprocessing tools

**Ninetails** requires the following input data to operate: \* multifast5 files basecalled by `Guppy` - for the signal data extraction \* sequencing_summary.txt file - for file ID extraction \* an output of `Nanopolish` polya function (tsv file) - to obtain the tail segmentation data

Therefore, please make sure that the third-party software necessary for the steps preceding the use of **Ninetails** is installed (`Nanopolish`, `Guppy`) and/or that you have all the required input files.

### Tensorflow & keras for R

The neural network in **ninetails** uses the tensorflow backend, so it is necessary to install it before running the program.

Instructions for installing `tensorflow` & `keras` can be found here: <https://tensorflow.rstudio.com/install/>

### HDF resources

Since fast5 is a binary format based on the HDF5, handling it requires installation of certain utilities in the operating system. Most of the necessary resources can be downloaded from the HDF5 Group website: <https://www.hdfgroup.org/downloads/hdf5/>.

**Ninetails** requires also `rhdf5` package and its dependencies for accessing & browsing files in fast5 format in R. It can be installed from Bioconductor (version available on CRAN is incompatible with newer R versions). The complete guide is available here: <https://bioconductor.org/packages/release/bioc/html/rhdf5.html>

``` r
install.packages("BiocManager")
BiocManager::install("rhdf5")
```

As a result, 3 packages are installed: `rhdf5`, `rhdf5filters` and `Rhdf5lib`. The `rhdf5filters` has to be reinstalled from the Github repo: <https://github.com/grimbough/rhdf5filters>, which can be done with following command in R/RStudio:

``` r
devtools::install_github('grimbough/rhdf5filters')
```

### VBZ compression plugin

Nanopore sequencing data produced with newer versions of the `MinKNOW` software are compressed with the VBZ algorithm which allows for more efficient compression than the formerly used GZIP. In order for the ninetails package to work with VBZ-compressed data, it is necessary to install an `ont-vbz-hdf-plugin` from the following repo: <https://github.com/nanoporetech/vbz_compression/releases>

Then it is necessary to navigate to the folder where the ONT VBZ compression plug-in was installed. On Linux operating systems, the default path is /usr/local/hdf5/lib/plugin/, whereas on Windows operating systems it might be C:\Program Files\OxfordNanopore\ont-vbz-hdf-plugin\hdf5\lib\plugin. The file(s) found there (libvbz_hdf_plugin.so on Linux and vbz_hdf_plugin.dll, vbz_hdf_plugin.lib on Windows) must be copied to the `rhdf5filters` folder (rhdf5filters/lib/) where the R libraries are stored.

### single fast5 support

Currently, **Ninetails** does not support single fast5 files as this format is deprecated by ONT. Before running the program on single fast5 files, you should convert them to multifast5 with another tool, for instance with `ont-fast5-api`.

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

## Usage

### Classification of reads using wrapper function

`check_tails()` is the main function which allows to classify sequencing reads based on presence/absence of non-adenosine residues within their poly(A) tails (and additional conditions, such as minimal read length and qc_tag assigned by Nanopolish polya function).

Below is an example of how to use `check_tails()` function:

``` r
results <- ninetails::check_tails(
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

### Classification of reads using standalone functions

The **Ninetails** pipeline may be also launched without the wrapper - as sometimes it might be useful, especially if the input files are large and/or you would like to plot some produced matrices.

The first function in processing pipeline is `create_tail_feature_list()`. It extracts the read data from the provided outputs and merges them based on read identifiers (readnames). This function works as follows:

``` r
tfl <- ninetails::create_tail_feature_list(
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
  pass_only=TRUE)
```

The second function, `create_tail_chunk_list()`, segments the reads and produces a list of segments in which a change of state (move = 1) along with significant local signal anomaly (so-called "pseudomove") has been recorded, possibly indicating the presence of a non-adenosine residue.

``` r
tcl <- ninetails::create_tail_chunk_list(tail_feature_list = tfl, 
                                         num_cores = 2)
```

The list of fragments should be then passed to the function `create_gaf_list()`, which transforms the signals into gramian angular fields (GAFs). The function outputs a list of arrays (100,100,2). First channel of each array consists of gramian angular summation field (GASF), while the second channel consists of gramian angular difference field (GADF).

``` r
gl <- ninetails::create_gaf_list(tail_chunk_list = tcl, 
                                 num_cores = 2)
```

The penultimate function, `predict_gaf_classes()`, launches the neural network to classify the input data. This function uses the tensorflow backend.

``` r
pl <- ninetails::predict_gaf_classes(gl)
```

The last function, `create_outputs()`, allows to obtain the final output: a list composed of **read_classes** (reads are labelled accordingly as "modified", "unmodified" and "unclassified" based on applied criteria) and **nonadenosine_residues** (detailed positional info regarding detected nonadenosine residues) data frames. Note that in this form the function does not automatically save data to files.

``` r
out <- ninetails::create_outputs(
  tail_feature_list = tfl,
  tail_chunk_list = tcl,
  nanopolish = system.file('extdata', 
                           'test_data', 
                           'nanopolish_output.tsv', 
                           package = 'ninetails'),
  predicted_list = pl,
  num_cores = 2,
  pass_only=TRUE)
```

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


## Citation

Please cite **Ninetails** as: Gumińska N., Matylla-Kulińska K., Krawczyk P., Orzeł W., Maj M., Mroczek S., Dziembowski A. Direct detection of non-adenosine nucleotides within poly(A) tails -- a new tool for the analysis of post-transcriptional mRNA tailing, 27th Annual Meeting of the RNA Society, Boulder, Colorado, May 31 to June 5, 2022.

Preprint is in the preparation.


## Troubleshooting

If you encounter a bug, please post it on github. To help diagnose the problem, send a minimal reproducible example (required inputs covering around 5-10 reads + corresponding nanopolish output & sequencing summary), so I will be able to reproduce the error and fix it for you.

## Maintainer

Any issues regarding the **Ninetails** should be addressed to Natalia Gumińska (nguminska (at) iimcb.gov.pl).
