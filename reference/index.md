# Package index

## Pipeline wrappers

Top-level functions for running the complete analysis pipeline. Choose
the appropriate wrapper based on your basecaller and data type.

- [`check_tails_dorado_DRS()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_DRS.md)
  : Complete Oxford Nanopore poly(A) tail analysis pipeline for Dorado
  DRS data.
- [`check_tails_dorado_cDNA()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_dorado_cDNA.md)
  : Complete Oxford Nanopore poly(A)/poly(T) tail analysis pipeline for
  Dorado cDNA data.
- [`check_tails_guppy()`](https://LRB-IIMCB.github.io/ninetails/reference/check_tails_guppy.md)
  : Wrapper function for complete DRS processing by ninetails package
  (legacy mode).

## Dorado DRS pipeline

Functions for processing direct RNA sequencing (DRS) data basecalled
with Dorado ≥ 1.0.0 in POD5 format.

- [`preprocess_inputs()`](https://LRB-IIMCB.github.io/ninetails/reference/preprocess_inputs.md)
  : Preprocess Dorado inputs for ninetails analysis (no BAM processing)
- [`process_dorado_summary()`](https://LRB-IIMCB.github.io/ninetails/reference/process_dorado_summary.md)
  : Process and split Dorado summary file into smaller parts
- [`filter_dorado_summary()`](https://LRB-IIMCB.github.io/ninetails/reference/filter_dorado_summary.md)
  : Filter Dorado summary for reads fulfilling ninetails quality
  criteria
- [`extract_tails_from_pod5()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tails_from_pod5.md)
  : Extract poly(A) tail signal segments from POD5 files using parallel
  Python processing
- [`create_tail_features_list_dorado()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_features_list_dorado.md)
  : Creates a nested list of Dorado tail features (raw signal +
  pseudomoves).
- [`create_tail_chunk_list_dorado()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_dorado.md)
  : Creates list of poly(A) tail chunks (Dorado mode) centered on
  significant signal deviations.
- [`split_tail_centered_dorado()`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_dorado.md)
  : Extracts fragments of poly(A) tail signal (Dorado mode) containing
  potential modifications along with their delimitation (positional
  indices; coordinates) within the tail.
- [`process_dorado_signal_files()`](https://LRB-IIMCB.github.io/ninetails/reference/process_dorado_signal_files.md)
  : Process Dorado poly(A) signal files for non-A prediction and tail
  chunk extraction
- [`create_outputs_dorado()`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado.md)
  : Create Ninetails output tables for Dorado DRS pipeline

## Dorado cDNA pipeline

Functions for processing cDNA sequencing data, including BAM file
processing, basecalled sequence extraction, and read orientation
classification (polyA vs polyT).

- [`preprocess_inputs_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/preprocess_inputs_cdna.md)
  : Preprocess Dorado inputs for ninetails cDNA analysis
- [`split_bam_file_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/split_bam_file_cdna.md)
  : Split BAM file into parts based on read IDs from summary file
- [`extract_data_from_bam()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_data_from_bam.md)
  : Extract data from BAM file for cDNA analysis
- [`detect_orientation_single()`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_single.md)
  : Detect poly tail type for a single sequence using Dorado-style
  algorithm
- [`detect_orientation_multiple()`](https://LRB-IIMCB.github.io/ninetails/reference/detect_orientation_multiple.md)
  : Classify multiple cDNA read orientations using Dorado-style poly
  tail detection
- [`process_polya_reads_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_reads_cdna.md)
  : Process polyA reads using standard ninetails pipeline
- [`process_polyt_reads_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/process_polyt_reads_cdna.md)
  : Process polyT reads using ninetails pipeline
- [`create_outputs_dorado_cdna()`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs_dorado_cdna.md)
  : Create Ninetails output tables for Dorado cDNA pipeline
- [`merge_cdna_results()`](https://LRB-IIMCB.github.io/ninetails/reference/merge_cdna_results.md)
  : Merge polyA and polyT processing results for cDNA analysis
- [`save_cdna_outputs()`](https://LRB-IIMCB.github.io/ninetails/reference/save_cdna_outputs.md)
  : Save cDNA pipeline outputs in standard ninetails format

## Guppy legacy pipeline

Functions for processing DRS data basecalled with Guppy ≤ 6.0.0 using
fast5 format and Nanopolish poly(A) coordinates. This pipeline is no
longer actively developed.

- [`extract_polya_data()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_polya_data.md)
  : Extract poly(A) data from nanopolish output and sequencing summary
- [`extract_tail_data()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data.md)
  : Extract tail features of a single RNA read from a multi-Fast5 file
- [`create_tail_feature_list()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list.md)
  : Create list of poly(A) tail features from multi-Fast5 files
- [`create_tail_chunk_list()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list.md)
  : Create list of poly(A) tail chunks centered on significant signal
  deviations
- [`split_tail_centered()`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered.md)
  : Extract modification-centered signal fragments from a poly(A) tail
- [`create_gaf()`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf.md)
  : Convert ONT signal to Gramian Angular Field
- [`create_gaf_list()`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list.md)
  : Create list of Gramian Angular Field matrices from tail chunks
- [`process_polya_complete()`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_complete.md)
  : Process a single (unsplit) poly(A) data file through the Guppy
  pipeline.
- [`process_polya_parts()`](https://LRB-IIMCB.github.io/ninetails/reference/process_polya_parts.md)
  : Process poly(A) data split into multiple parts through the Guppy
  pipeline.
- [`split_polya_data()`](https://LRB-IIMCB.github.io/ninetails/reference/split_polya_data.md)
  : Split large poly(A) data file into smaller parts.
- [`create_outputs()`](https://LRB-IIMCB.github.io/ninetails/reference/create_outputs.md)
  : Create ninetails output tables (Guppy legacy pipeline)
- [`save_outputs()`](https://LRB-IIMCB.github.io/ninetails/reference/save_outputs.md)
  : Save pipeline outputs to files.

## Training dataset production

Functions for preparing training and validation datasets for the
convolutional neural network (CNN) model.

- [`prepare_trainingset()`](https://LRB-IIMCB.github.io/ninetails/reference/prepare_trainingset.md)
  : Filters out signals of a given nucleotide type for neural network
  training-set preparation.
- [`extract_tail_data_trainingset()`](https://LRB-IIMCB.github.io/ninetails/reference/extract_tail_data_trainingset.md)
  : Extracts tail features of single RNA read from respective basecalled
  multi-fast5 file.
- [`create_tail_feature_list_trainingset()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_trainingset.md)
  : Extracts features of poly(A) tails of ONT RNA reads required for
  finding non-A nucleotides within the given tails.
- [`create_tail_feature_list_A()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_feature_list_A.md)
  : Extracts features of poly(A) tails containing only A nucleotides for
  training-set preparation.
- [`create_tail_chunk_list_trainingset()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_trainingset.md)
  : Extracts decoration-centered fragments of poly(A) tails for all
  reads and appends positional data to a nested list.
- [`create_tail_chunk_list_A()`](https://LRB-IIMCB.github.io/ninetails/reference/create_tail_chunk_list_A.md)
  : Creates list of tail chunks containing only A nucleotides.
- [`split_tail_centered_trainingset()`](https://LRB-IIMCB.github.io/ninetails/reference/split_tail_centered_trainingset.md)
  : Extracts decoration-centered fragments of poly(A) tail signal along
  with positional coordinates.
- [`split_with_overlaps()`](https://LRB-IIMCB.github.io/ninetails/reference/split_with_overlaps.md)
  : Splits signal to overlapping fragments of equal length.
- [`filter_nonA_chunks_trainingset()`](https://LRB-IIMCB.github.io/ninetails/reference/filter_nonA_chunks_trainingset.md)
  : Filters read chunks containing non-adenosine nucleotides of interest
  for neural network training-set preparation.
- [`filter_signal_by_threshold_trainingset()`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold_trainingset.md)
  : Detection of outliers (peaks & valleys) in ONT signal using
  z-scores.
- [`create_gaf_list_A()`](https://LRB-IIMCB.github.io/ninetails/reference/create_gaf_list_A.md)
  : Produces list of GAFs containing exclusively A-nucleotides for
  neural network training.

## Data postprocessing

Functions for correcting, reclassifying, and reshaping ninetails output
tables after the pipeline has run.

- [`correct_class_data()`](https://LRB-IIMCB.github.io/ninetails/reference/correct_class_data.md)
  : Corrects the classification of reads contained in the class_data
  table.
- [`correct_residue_data()`](https://LRB-IIMCB.github.io/ninetails/reference/correct_residue_data.md)
  : Marks uncertain positions of non-A residues in ninetails output
  data.
- [`correct_labels()`](https://LRB-IIMCB.github.io/ninetails/reference/correct_labels.md)
  : Correct read class labels for backward compatibility
- [`reclassify_ninetails_data()`](https://LRB-IIMCB.github.io/ninetails/reference/reclassify_ninetails_data.md)
  : Reclassifies ambiguous non-A residues to mitigate potential errors
  inherited from nanopolish segmentation.
- [`read_class_single()`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_single.md)
  : Reads ninetails read_classes data frame from file.
- [`read_class_multiple()`](https://LRB-IIMCB.github.io/ninetails/reference/read_class_multiple.md)
  : Reads multiple ninetails read_classes outputs at once.
- [`read_residue_single()`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_single.md)
  : Reads ninetails nonadenosine_residues data from file.
- [`read_residue_multiple()`](https://LRB-IIMCB.github.io/ninetails/reference/read_residue_multiple.md)
  : Reads multiple ninetails nonadenosine_residues outputs at once.
- [`merge_nonA_tables()`](https://LRB-IIMCB.github.io/ninetails/reference/merge_nonA_tables.md)
  : Merges ninetails tabular outputs (read classes and nonadenosine
  residue data) to produce one concise table.
- [`spread_nonA_residues()`](https://LRB-IIMCB.github.io/ninetails/reference/spread_nonA_residues.md)
  : Reshapes nonadenosine_residues data frame to wide format.

## Annotation

Functions for biological annotation of ninetails results using external
databases.

- [`annotate_with_biomart()`](https://LRB-IIMCB.github.io/ninetails/reference/annotate_with_biomart.md)
  : Annotate ninetails output data with biomaRt

## Statistics

Functions for statistical analysis and quantification of non-adenosine
residues across reads and conditions.

- [`calculate_fisher()`](https://LRB-IIMCB.github.io/ninetails/reference/calculate_fisher.md)
  : Perform Fisher's exact test per transcript with BH p-value
  adjustment
- [`nonA_fisher()`](https://LRB-IIMCB.github.io/ninetails/reference/nonA_fisher.md)
  : Perform Fisher's exact test on a single transcript in ninetails
  output
- [`count_class()`](https://LRB-IIMCB.github.io/ninetails/reference/count_class.md)
  : Counts read classes found in a read_classes data frame produced by
  the ninetails pipeline.
- [`count_nonA_abundance()`](https://LRB-IIMCB.github.io/ninetails/reference/count_nonA_abundance.md)
  : Counts reads by number of non-A occurrence instances.
- [`count_residues()`](https://LRB-IIMCB.github.io/ninetails/reference/count_residues.md)
  : Counts non-A residues found in a nonadenosine_residues data frame
  produced by the ninetails pipeline.
- [`summarize_nonA()`](https://LRB-IIMCB.github.io/ninetails/reference/summarize_nonA.md)
  : Produces summary table of non-A occurrences within an analyzed
  dataset.
- [`nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/nanopolish_qc.md)
  : Aggregates nanopolish polya quality control information.

## Visualisation

Plotting functions for inspection of raw signals, GAF images,
classification results, and statistical summaries.

- [`plot_class_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_class_counts.md)
  : Plotting read classes data per category assigned to the analyzed
  reads.
- [`plot_gaf()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_gaf.md)
  : Creates a visual representation of gramian angular field
  corresponding to the given poly(A) tail fragment (chunk).
- [`plot_multiple_gaf()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_multiple_gaf.md)
  : Creates a visual representation of multiple gramian angular fields
  based on provided gaf_list (plots all gafs from the given list).
- [`plot_nanopolish_qc()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nanopolish_qc.md)
  : Plots qc data (qc_tag) inherited from nanopolish polya function.
- [`plot_nonA_abundance()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_nonA_abundance.md)
  : Plot abundances of reads with given amount of non-A residues per
  read
- [`plot_panel_characteristics()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_panel_characteristics.md)
  : Plot panel characteristics of ninetails output
- [`plot_residue_counts()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_residue_counts.md)
  : Plot counts of nonadenosine residues found in ninetails output data
- [`plot_rug_density()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_rug_density.md)
  : Scatterplot of nonA residue positions within poly(A) tail
- [`plot_squiggle_fast5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_squiggle_fast5.md)
  : Draws an entire squiggle for given read.
- [`plot_squiggle_pod5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_squiggle_pod5.md)
  : Draws an entire squiggle for given read from POD5 file.
- [`plot_tail_chunk()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_chunk.md)
  : Draws a portion of poly(A) tail squiggle (chunk) for given read.
- [`plot_tail_distribution()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_distribution.md)
  : Plots poly(A) tail length (or estimated non-A position) distribution
  in analyzed sample(s).
- [`plot_tail_range_fast5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_range_fast5.md)
  : Draws tail range squiggle for given read.
- [`plot_tail_range_pod5()`](https://LRB-IIMCB.github.io/ninetails/reference/plot_tail_range_pod5.md)
  : Draws tail range squiggle for given read from POD5 file.

## tailfindr compatibility

Functions for converting tailfindr output into a format compatible with
the ninetails Guppy legacy pipeline.

- [`convert_tailfindr_output()`](https://LRB-IIMCB.github.io/ninetails/reference/convert_tailfindr_output.md)
  : Converts tailfindr results to format compatible with ninetails
- [`check_polya_length_filetype()`](https://LRB-IIMCB.github.io/ninetails/reference/check_polya_length_filetype.md)
  : Check and convert poly(A) length file format

## Signal processing

Core signal processing utilities and CNN-related helpers used internally
by the pipeline functions.

- [`filter_signal_by_threshold()`](https://LRB-IIMCB.github.io/ninetails/reference/filter_signal_by_threshold.md)
  : Detect outliers (peaks and valleys) in ONT signal using z-scores
- [`winsorize_signal()`](https://LRB-IIMCB.github.io/ninetails/reference/winsorize_signal.md)
  : Winsorize nanopore signal
- [`substitute_gaps()`](https://LRB-IIMCB.github.io/ninetails/reference/substitute_gaps.md)
  : Substitute short zero-gaps surrounded by nonzero pseudomoves
- [`combine_gafs()`](https://LRB-IIMCB.github.io/ninetails/reference/combine_gafs.md)
  : Combine GASF and GADF into a two-channel array
- [`predict_gaf_classes()`](https://LRB-IIMCB.github.io/ninetails/reference/predict_gaf_classes.md)
  : Classify Gramian Angular Field matrices with a pretrained CNN
- [`load_keras_model()`](https://LRB-IIMCB.github.io/ninetails/reference/load_keras_model.md)
  : Load Keras model for multiclass signal prediction

## Sequence helpers

Helper functions for primer matching and DNA sequence manipulation used
by the cDNA orientation classification step.

- [`reverse_complement()`](https://LRB-IIMCB.github.io/ninetails/reference/reverse_complement.md)
  : Generate reverse complement of a DNA sequence
- [`edit_distance_hw()`](https://LRB-IIMCB.github.io/ninetails/reference/edit_distance_hw.md)
  : Calculate edit distance with sliding window (HW mode)
- [`count_trailing_chars()`](https://LRB-IIMCB.github.io/ninetails/reference/count_trailing_chars.md)
  : Count trailing occurrences of a character in a string

## Input validation

Internal assertion and type-checking utilities used throughout the
package for input validation.

- [`assert_condition()`](https://LRB-IIMCB.github.io/ninetails/reference/assert_condition.md)
  : Assert condition is TRUE, stop with message if FALSE
- [`assert_dir_exists()`](https://LRB-IIMCB.github.io/ninetails/reference/assert_dir_exists.md)
  : Assert directory exists with informative error
- [`assert_file_exists()`](https://LRB-IIMCB.github.io/ninetails/reference/assert_file_exists.md)
  : Assert file exists with informative error
- [`check_fast5_filetype()`](https://LRB-IIMCB.github.io/ninetails/reference/check_fast5_filetype.md)
  : Check if the provided directory contains Fast5 files in the correct
  format
- [`check_output_directory()`](https://LRB-IIMCB.github.io/ninetails/reference/check_output_directory.md)
  : Check and handle existing output directory for ninetails analysis
- [`is_RNA()`](https://LRB-IIMCB.github.io/ninetails/reference/is_RNA.md)
  : Check if fast5 file contains RNA reads
- [`is_multifast5()`](https://LRB-IIMCB.github.io/ninetails/reference/is_multifast5.md)
  : Check if fast5 file is multi-read format
- [`is_string()`](https://LRB-IIMCB.github.io/ninetails/reference/is_string.md)
  : Test if x is a single non-empty character string
- [`no_na()`](https://LRB-IIMCB.github.io/ninetails/reference/no_na.md)
  : Check for no NA values
- [`get_mode()`](https://LRB-IIMCB.github.io/ninetails/reference/get_mode.md)
  : Calculate the statistical mode of a numeric vector

## Package

- [`ninetails`](https://LRB-IIMCB.github.io/ninetails/reference/ninetails-package.md)
  [`ninetails-package`](https://LRB-IIMCB.github.io/ninetails/reference/ninetails-package.md)
  : ninetails: Nonadenosine Nucleotides in Poly(A) Tails
- [`` `%>%` ``](https://LRB-IIMCB.github.io/ninetails/reference/pipe.md)
  : Pipe operator
