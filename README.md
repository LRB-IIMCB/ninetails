# ninetails

**An R package for finding non-adenosine poly(A) residues in Oxford Nanopore direct RNA sequencing reads**

<img align="right" width="200" height="220" src="https://user-images.githubusercontent.com/68285258/168554098-a5a5dee9-2c8f-4351-86b4-e6420a5b8ced.png">

## Introduction
* It works on Oxford Nanopore direct RNA sequencing reads basecalled by Guppy software
* It requires tail delimitation data produced by Nanopolish software
* It allows both for the detection of non-adenosine residues within the poly(A) tails and visual inspection of  read signals

The software is still under development, so all suggestions to improving it are welcome. Please note that the code contained herein may change frequently, so use it with caution.

## Prerequisites

Ninetails requires the following input data to operate:
* multifast5 files basecalled by ```Guppy``` - for the signal data extraction 
* sequencing_summary.txt file - for file ID extraction
* an output of ```Nanopolish``` polya function (tsv file) - to obtain the tail segmentation data

Therefore, please make sure that the third-party software necessary for the steps preceding the use of **ninetails** is installed (```Nanopolish```, ```Guppy```) and/or that you have all the required input files.

Currently, **ninetails** does not support single fast5 files as this format is deprecated by ONT. Before running the program on single fast5 files, you should convert them to multifast5 with another tool, for instance with ```ont-fast5-api```.


## Installation

Currently, **ninetails** is not available on CRAN/Bioconductor, so you need to install it using ```devtools```.

If you do not have ```devtools``` installed already, you can do this with the following command in R/RStudio:

```r
install.packages("devtools")
```
Once you have ```devtools``` installed, you can install **ninetails** using the command below in R/RStudio:

```r
devtools::install_github('LRB-IIMCB/ninetails')
library(ninetails)
```

## Usage

### Classification of reads using wrapper function

```check_tails()``` is the main function which allows to classify sequencing reads based on presence/absence of non-adenosine residues within their poly(A) tails (and additional conditions, such as minimal read length and qc_tag assigned by Nanopolish polya function). 

Below is an example of how to use ```check_tails()``` function:

```r
library(ninetails)

results <- check_tails(nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', 
                                                package = 'ninetails'),
                       sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt', 
                                                        package = 'ninetails'),
                       workspace = system.file('extdata', 'test_data', 'basecalled_fast5', 
                                               package = 'ninetails'),
                       num_cores = 2,
                       basecall_group = 'Basecall_1D_000',
                       pass_only=TRUE,
                       save_dir = '~/Downloads')
```

This function returns a list consisting of two tables: **binary_classified_reads** and **detailed_positional_nonadenosine_residues**. In addition, the function saves results to text files in the user-specified directory. 

Moreover, the function also creates a log file in the directory specified by the user.

### Classification of reads using standalone functions

### Output explanation

The **binary_classified_reads** dataframe (file) contains following columns:

| column name  | content |
| ------------- | ------------- |
| readname  | an identifier of a given read  (36 characters)|
| class  | the crude result of classification  |

The ```class``` column contains information whether the given read was recognized as modified (containing non-adenosine residue) or not. Additionally, it contains information about reads discarded from classification - the reason for this is also included. This column may contain a following content:

| class column content  | explanation |
| ------------- | ------------- |
| modified | read recognized as containing non-adenosine residue within poly(A) tail|
|  unmodified | read containing only adenosines in poly(A) tail  |
|  unclassified - insufficient read length | read length estimated by Nanopolish as <10 nt (insufficient for proper tail segmentation)|
|  unclassified - nanopolish qc failed (adapter) | insufficient read quality assigned by Nanopolish |
|  unclassified - nanopolish qc failed (suffclip) | insufficient read quality assigned by Nanopolish. This category shows only if ```pass_only=TRUE``` was passed to the processing function, so reads tagged by Nanopolish as "suffclip" were excluded from the analysis|

The **detailed_positional_nonadenosine_residues** dataframe (file) contains following columns:




### Visual inspection of reads of interest

**Ninetails** has built-in functions ```plot_squiggle()``` and ```plot_tail_range()``` for plotting whole reads and the poly(A) tail region, respectively.

With their help you can visualise:
* raw signal (```rescale=FALSE```)
* signal scaled to picocamperes [pA] per second [s] (```rescale=TRUE```)

In addition, the graphs can depict only signal data (```moves=FALSE```) or they can also contain informations regarding the change of state between the subsequent k-mers (moves) (```moves=TRUE```).

#### Plotting an entire signal

The following example demonstrates how to use the ```plot_squiggle()``` function:

```r
library(ninetails)

plot <- plot_squiggle(readname = "cd6c7f1d-6ef4-45a0-a67a-0a2853967e5b",
                      nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', 
                                                package = 'ninetails'),
                      sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt', 
                                                        package = 'ninetails'),
                      workspace = system.file('extdata', 'test_data', 'basecalled_fast5', 
                                               package = 'ninetails'),
                      basecall_group = 'Basecall_1D_000',
                      moves = FALSE,
                      rescale = TRUE)

print(plot)
```

The output of the above command is the following graph:
![moves_false_whole_sig](https://user-images.githubusercontent.com/68285258/170456526-e4b05d2a-1fda-45f4-b1f2-f4dfa82b751c.png)

If the (```moves=TRUE```) argument is passed, then the graph also contains information regarding moves, which looks as follows:
![moves_true_whole_sig](https://user-images.githubusercontent.com/68285258/170457349-d3adcf55-37e3-4a70-b8e1-689ed9f05182.png)



#### Plotting tail range

The ```plot_tail_range()``` function accepts the same arguments as the abovementioned function. An example usage looks as follows:

```r

library(ninetails)

plot <- plot_tail_range(readname = "cd6c7f1d-6ef4-45a0-a67a-0a2853967e5b",
                        nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', 
                                                 package = 'ninetails'),
                        sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt', 
                                                         package = 'ninetails'),
                        workspace = system.file('extdata', 'test_data', 'basecalled_fast5', 
                                                package = 'ninetails'),
                        basecall_group = 'Basecall_1D_000',
                        moves = FALSE,
                        rescale = TRUE)

print(plot)
```
Which outputs:
![moves_false_tail_sig](https://user-images.githubusercontent.com/68285258/170458080-1c550bc3-488c-4c14-8390-688d613f2ca4.png)

Or below, if the (```moves=TRUE```) argument is passed:
![moves_true_tail_sig](https://user-images.githubusercontent.com/68285258/170458256-c850e981-b8b0-4c29-bce0-2094ab9136e2.png)




## Citation

Please cite **ninetails** as: Gumińska N et al., Direct detection of non-adenosine nucleotides within poly(A) tails – a new tool for the analysis of post-transcriptional mRNA tailing

Preprint is in the preparation.

## Future plans

## Troubleshooting

If you encounter a bug, please post it on github. To help diagnose the problem, send a minimal reproducible example (required inputs covering around 5-10 reads), so I will be able to reproduce the error and fix it for you.

## Maintainer

Any issues regarding the **ninetails** should be addressed to Natalia Gumińska (nguminska (at) iimcb.gov.pl).

