# ninetails

**An R package for finding non-adenosine poly(A) residues in Oxford Nanopore direct RNA sequencing reads**

<img align="right" width="200" height="220" src="https://user-images.githubusercontent.com/68285258/168554098-a5a5dee9-2c8f-4351-86b4-e6420a5b8ced.png">

## Introduction
* It works on Oxford Nanopore direct RNA sequencing reads basecalled by Guppy software
* It requires tail delimitation data produced by Nanopolish software
* It allows both for the detection of non-adenosine residues within the poly(A) tails and visual inspection of  read signals

Currently, **ninetails** can distinguish characteristic signatures of four types of nucleotides: adenosines (A), cytosines (C), guanosines (G), and uridines (U).

**Ninetails** relies on Nanopolish segmentation and therefore may underestimate terminal modifications (last and penultimate nucleotides of the tail).

The software is still under development, so all suggestions to improving it are welcome. Please note that the code contained herein may change frequently, so use it with caution.

## Important note

Please be aware that signal transformations performed during analysis can place a heavy load on memory. This is especially true if your data covers the entire sequencing run. For the moment, **ninetails** does not offer the possibility of processing large data sets in chunks behind the scenes (under development). Therefore, to minimise the risk of unexpected crashes, it is highly recommended to split the output of the ```Nanopolish``` polyA function into smaller files to make it easier to process the data in subsets and then merge the final results.

## Prerequisites

**Ninetails** requires the following input data to operate:
* multifast5 files basecalled by ```Guppy``` - for the signal data extraction 
* sequencing_summary.txt file - for file ID extraction
* an output of ```Nanopolish``` polya function (tsv file) - to obtain the tail segmentation data

Therefore, please make sure that the third-party software necessary for the steps preceding the use of **ninetails** is installed (```Nanopolish```, ```Guppy```) and/or that you have all the required input files.

Currently, **ninetails** does not support single fast5 files as this format is deprecated by ONT. Before running the program on single fast5 files, you should convert them to multifast5 with another tool, for instance with ```ont-fast5-api```.

The neural network in **ninetails** uses the tensorflow backend, so it is necessary to install this package before running the program. 

Instructions for installing tensorflow can be found here:
https://tensorflow.rstudio.com/installation/


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

The **ninetails** pipeline may be also launched without the wrapper - as sometimes it might be useful, especially if the input files are large and/or you would like to plot some produced matrices. 

The first function in processing pipeline is ```create_tail_feature_list()```. It extracts the read data from the provided outputs and merges them based on read identifiers (readnames).  This function works as follows:

```r
tfl <- create_tail_feature_list(nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', 
                                                         package = 'ninetails'),
                                sequencing_summary = system.file('extdata', 'test_data', 'sequencing_summary.txt',
                                                                 package = 'ninetails'),
                                workspace = system.file('extdata', 'test_data', 'basecalled_fast5', 
                                                        package = 'ninetails'),
                                num_cores = 10,
                                basecall_group = 'Basecall_1D_000',
                                pass_only=TRUE)
```

The second function, ```create_tail_chunk_list_moved()```,  segments the reads and produces a list of segments in which a change of state (move = 1) has been recorded, possibly indicating the presence of a non-adenosine residue.

```r
tcl <- create_tail_chunk_list_moved(tail_feature_list = tfl, 
                                    num_cores = 10)
```
The list of fragments should be then passed to the function ```create_gasf_list()```, which transforms the signals into gramian angular summation fields. The function outputs a list of GASF matrices. 

```r
gl <- create_gasf_list(tail_chunk_list = tcl, 
                       num_cores = 10)
```

The list of matrices should then be passed to the ```predict_classes()``` function, which launches the neural network to classify the input data. This function uses the tensorflow backend.

```r
pl <- predict_classes(gl)
```

The penultimate function, ```create_coordinate_dataframe()```,  creates a table with data necessary to estimate the position of individual modifications within the tails harboring them. 

```r
cdf <- create_coordinate_dataframe(tail_feature_list=tfl,
                                   num_cores=10)
```
The last function, ```analyze_results()```,  allows to obtain the final output: a list composed of **binary_classified_reads** and **detailed_positional_nonadenosine_residues** data frames. Note that in this form the function does not automatically save data to files.

```r
results <- analyze_results(nanopolish = system.file('extdata', 'test_data', 'nanopolish_output.tsv', 
                                                    package = 'ninetails'), 
                           coordinate_df=cdf, 
                           predicted_list=pl, 
                           pass_only=TRUE)
```






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

| column name  | content |
| ------------- | ------------- |
| readname  | an identifier of a given read  (36 characters)|
| chunk	  | number of tail segment in which the given non-adenosine residue was spotted; reported from 3' end  |
| prediction  | the result of classification (basic model: C, G, U assignment)  |
| total_chunk  | total number of segments per given tail |
| tail_length  | the tail length estimated according to Nanopolish polya function  |
| interval  | the approximate region of the poly(A) tail where the modification is to be expected |
| centered_pos  | the approximate nucleotide position in the centre of a given tail segment where modifications are to be expected; reported from 5' end |



### Visual inspection of reads of interest

**Ninetails** has built-in functions ```plot_squiggle()``` and ```plot_tail_range()``` for plotting whole reads and the poly(A) tail region, respectively.

With their help you can visualise:
* raw signal (```rescale=FALSE```)
* signal scaled to picoamperes [pA] per second [s] (```rescale=TRUE```)

In addition, the graphs can depict only signal data (```moves=FALSE```) or they can also contain informations regarding the change of state between the subsequent k-mers (moves) (```moves=TRUE```).

#### Plotting an entire signal

The following example demonstrates how to use the ```plot_squiggle()``` function:

```r
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

#### Plotting the tail segment of interest

Ninetails allows you to visualise a particular fragment among the list of fragments generated by the ```create_tail_chunk_list()``` function. This is what the function ```plot_tail_chunk()``` is for. This function only allows to preview the raw signal, currently there is no built-in scaling to picoamperes [pA]. 

An example of how the ```plot_tail_chunk()``` function works:
```r

example <- plot_tail_chunk(chunk_name = "1625f71d-287e-42b0-ae37-dc14c2f4ed8e_5",
                          tail_chunk_list = tcl)

print(example)
```

And an example output:

![00001a](https://user-images.githubusercontent.com/68285258/170530349-fd36cf05-776d-4236-989c-ca7c22492c1e.png)


### Plotting the gramian angular summation fields

The package allows to create a visual representation of gramian angular summation fields (GASFs) using ```ggplot2```. 

#### Plotting single GASF of interest

The ```plot_gasf()``` function draws a single matrix of interest. It requires the name of a particular segment and a list of matrices produced by the ```create_gasf_list()``` function as an input. 

Below is an example of the usage of the ```plot_gasf()``` function. Please note that in order for this example to work properly, one must first execute the 3 first commands from the **Classification of reads using standalone functions** section.  

```r

example_gasf <- plot_gasf(gasf_name = "1625f71d-287e-42b0-ae37-dc14c2f4ed8e_5",
                          gasf_list = gl, 
                          save_file = TRUE)

print(example_gasf)
```
And here is an example output:


![1625f71d-287e-42b0-ae37-dc14c2f4ed8e_5](https://user-images.githubusercontent.com/68285258/170500508-b27289fa-e62c-4a14-9e50-ff1a8798dfa7.png)



#### Plotting multiple GASFs

**Ninetails** also allows the user to plot the entire list of matrices produced by the ```create_gasf_list()``` function at once. The files will be saved in the current working directory. An example of usage is given below:

```r
plot_multiple_gasf(gasf_list = gl, num_cores = 10)

```
And this results in multiple plots, like this: 

![drawing](https://user-images.githubusercontent.com/68285258/170512100-ddfb03dc-63c8-449e-8339-57a4ec16ce43.png)





However, it is advisable to use this function with caution. The data contained in a ```gasf_list``` object tends to be large. Drawing multiple graphs at once may overload the system and cause it to crash.


## Citation

Please cite **ninetails** as: Gumińska N., Matylla-Kulińska K., Krawczyk P., Orzeł W., Maj M., Mroczek S., Dziembowski A. (2022), Direct detection of non-adenosine nucleotides within poly(A) tails – a new tool for the analysis of post-transcriptional mRNA tailing, 

Preprint is in the preparation.

## Future plans
* model finetuning
* optimized positioning (position calibration)
* port to Python

## Troubleshooting

If you encounter a bug, please post it on github. To help diagnose the problem, send a minimal reproducible example (required inputs covering around 5-10 reads + corresponding nanopolish output & sequencing summary), so I will be able to reproduce the error and fix it for you.

## Maintainer

Any issues regarding the **ninetails** should be addressed to Natalia Gumińska (nguminska (at) iimcb.gov.pl).

