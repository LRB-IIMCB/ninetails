#' Reads ninetails read_classes dataframe from file
#'
#' This is the basic function used to import read_classes output from check_tails
#' to R
#'
#' @param class_path a character string. Path to ninetails output file
#'
#' @param sample_name a sample identifier (optional), provided as a character
#' string. If specified will be included as an additional column sample_name.
#'
#' Function based on PK (smaegol) NanoTail read_polya_single
#' https://github.com/LRB-IIMCB/nanotail/
#' Many thanks to the NanoTail developer for help and advice!
#'
#' @return a [tibble] with read_classes predictions
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#'class_path <- "/directory/with/ninetails/read_class_output.tsv"
#'class_data <- ninetails::read_class_single(class_path)
#'}
read_class_single <- function(class_path, sample_name = NA) {

  # assertions
  if (missing(class_path)) {
    stop("The path to class predictions (argument class_path) is missing",
         call. = FALSE)
  }
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(class_path),
                          msg = "Empty string provided as an input. Please provide a class_path as a string")
  assertthat::assert_that(assertive::is_existing_file(class_path),
                          msg=paste("File ",class_path," does not exist",sep=""))
  assertthat::assert_that(assertive::is_non_empty_file(class_path),
                          msg=paste("File ",class_path," is empty",sep=""))

  # load class data
  message(paste0("Loading data from ",class_path))
  class_data <- vroom::vroom(class_path, show_col_types = FALSE) %>% dplyr::as_tibble()

  if(!is.na(sample_name)) {
    # set sample_name (if was set)
    if (! "sample_name" %in% colnames(class_data)) {
      warning("The sample_name was provided in the input file. Overwriting with the provided one")
    }
    class_data$sample_name = sample_name
    class_data$sample_name <- as.factor(class_data$sample_name)
  }

  return(class_data)
}



#' Reads multiple ninetails read_classes outputs at once
#'
#' This function can be used to load any number of files with read_classes
#' predictions with single invocation, allowing for metadata specification.
#'
#' @param samples_table data.frame or tibble containing samples metadata
#' and paths to files.
#' Should have at least two columns: \itemize{
#' \item class_path - containing path to the read_classes predictions file
#' \item sample_name - unique sample identifier
#' }
#'
#' @param ... - additional parameters to pass to read_class_single();
#' currently accepts sample_name argument
#'
#' Function based on PK (smaegol) NanoTail read_polya_multiple
#' https://github.com/LRB-IIMCB/nanotail/
#' Many thanks to the NanoTail developer for help and advice!
#'
#' @return a [tibble] containing read_classes data for
#' all specified samples, with metadata provided in samples_table
#' stored as separate columns
#'
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' samples_table <- data.frame(class_path = c(path_1,path_2,path_3,path_4),
#'                             sample_name =c("wt_1","mut_1","wt_2","mut_2"),
#'                             group = c("wt","mut","wt","mut"))
#'
#' classes_data <- ninetails::read_class_multiple(samples_table)
#'}
read_class_multiple <- function(samples_table,...) {

  # var binding
  sample_name <- class_path <- class_contents <- NULL

  if (missing(samples_table)) {
    stop("Samples table argument is missing",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(samples_table),
                          msg = "Empty data frame provided as an input (samples_table). Please provide samples_table describing data to load")
  assertthat::assert_that("class_path" %in% colnames(samples_table),
                          msg = "The samples_table should contain at least class_path and sample_name columns")
  assertthat::assert_that("sample_name" %in% colnames(samples_table),
                          msg = "The samples_table should contain at least class_path and sample_name columns")

  samples_data <- samples_table %>%
    dplyr::as.tbl() %>%
    dplyr::mutate_if(is.character,as.factor) %>%
    dplyr::mutate(class_path = as.character(class_path)) %>%
    dplyr::group_by(sample_name) %>% dplyr::mutate(class_contents=purrr::map(class_path, function(x) ninetails::read_class_single(x))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-class_path)

  class_data <- tidyr::unnest(samples_data, cols = c(class_contents))

  return(class_data)
}

#' Counts read classes found in read_classes dataframe produced
#' by ninetails pipeline
#'
#' Process the information returned by ninetails found in the read_classes
#' dataframe in the prediction column
#'
#' @param class_data A dataframe or tibble containig read_classes predictions
#' made by ninetails pipeline
#'
#' @param grouping_factor character string. A grouping variable
#' (e.g. "sample_name")
#'
#' @param detailed logical [TRUE/FALSE]. If TRUE, the counts will be provided
#' based on the "comments" column, which contains detailed information on the
#' assigned class. If FALSE, the counts will be provided based on "class" column
#' which gives more crude glimpse on the classification - i.e. provides an info
#' whether the reads were considered as "modified", "unmodified" and
#' "unclassified" only. By default, the TRUE option is set.
#'
#'
#' @return A tibble with counts for each non-A residue
#'
#' @export
#'
#' @examples
#'\dontrun{
#' class_counted <- ninetails::count_class(class_data=out[[1]],
#'                                         grouping_factor=NA,
#'                                         detailed=TRUE)
#'}
count_class <- function(class_data, grouping_factor=NA, detailed=TRUE) {
  # variable binding
  comments <- NULL

  #assertions
  assertthat::assert_that(assertive::has_rows(class_data),
                          msg = "Empty dataframe provided as an input")

  if (detailed==TRUE){
    if(!is.na(grouping_factor)) {
      assertthat::assert_that(grouping_factor %in% colnames(class_data),
                              msg=paste0(grouping_factor," is not a column of input dataset"))
      class_counts <- class_data %>%
        dplyr::mutate(comments=forcats::fct_relevel(comments,"move transition present, nonA residue detected", after = Inf)) %>%
        dplyr::group_by(!!rlang::sym(grouping_factor),comments) %>%
        dplyr::count()
    } else {
      class_counts <- class_data %>%
        dplyr::mutate(comments=forcats::fct_relevel(comments,"move transition present, nonA residue detected", after = Inf)) %>%
        dplyr::group_by(comments) %>%
        dplyr::count()
    }

  } else {
    if(!is.na(grouping_factor)) {
      assertthat::assert_that(grouping_factor %in% colnames(class_data),
                              msg=paste0(grouping_factor," is not a column of input dataset"))
      class_counts <- class_data %>%
        dplyr::mutate(class=forcats::fct_relevel(class,"modified", after = Inf)) %>%
        dplyr::group_by(!!rlang::sym(grouping_factor),class) %>%
        dplyr::count()
    } else {
      class_counts <- class_data %>%
        dplyr::mutate(class=forcats::fct_relevel(class,"modified", after = Inf)) %>%
        dplyr::group_by(class) %>%
        dplyr::count()
    }
  }

  return (class_counts)
}


#' Reads ninetails nonadenosine_residues data from file
#'
#' This is the basic function used to import output from check_tails to R:
#' it imports the dataframe containing non-A detailed position info
#' (nonadenosine_residues).
#'
#' @param residue_path path to ninetails nonadenosine_residues output file
#'
#' @param sample_name sample name (optional), provided as a string.
#' If specified will be included as an additional column sample_name.
#'
#' Function based on PK (smaegol) NanoTail read_polya_single
#' https://github.com/LRB-IIMCB/nanotail/
#' Many thanks to the NanoTail developer for help and advice!
#'
#' @export
#'
#' @return a [tibble] with non-A predictions
#'
#' @examples
#'\dontrun{
#' residue_path <- "/directory/with/ninetails/nonadenosine_residues_output.tsv"
#' residue_data <- ninetails::read_residue_single(residue_path)
#'}
read_residue_single <- function(residue_path, sample_name = NA) {

  #assertions
  if (missing(residue_path)) {
    stop("The path to residue predictions (argument residue_path) is missing",
         call. = FALSE)
  }
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(residue_path),
                          msg = "Empty string provided as an input. Please provide a residue_path as a string")
  assertthat::assert_that(assertive::is_existing_file(residue_path),
                          msg=paste("File ",residue_path," does not exist",sep=""))
  assertthat::assert_that(assertive::is_non_empty_file(residue_path),
                          msg=paste("File ",residue_path," is empty",sep=""))

  # load the data
  message(paste0("Loading non-A residue data from ", residue_path))
  residue_data <- vroom::vroom(residue_path, show_col_types = FALSE) %>% dplyr::as_tibble()


  if(!is.na(sample_name)) {
    # set sample_name (if was set)
    if (! "sample_name" %in% colnames(residue_data)) {
      warning("sample_name was provided in the input file. Overwriting with the provided one")
    }
    residue_data$sample_name = sample_name
    residue_data$sample_name <- as.factor(residue_data$sample_name)
  }

  return(residue_data)
}


#' Reads multiple ninetails nonadenosine_residues outputs at once
#'
#' This function can be used to load any number of files
#' with nonadenosine_residues predictions with single invocation,
#' allowing for metadata specification.
#'
#' It requires providing a samples_table with metadata necessary
#' to identify each sample.
#'
#' @param samples_table data.frame or tibble containing samples metadata and
#' paths to files.
#' It should contain at least two columns: \itemize{
#' \item residue_path - containing path to the nonadenosine_residues predictions file
#' \item sample_name - unique sample identifier
#' }
#' @param ... - additional parameters to pass to read_residue_single();
#' currently accepts sample_name argument
#'
#' Function based on PK (smaegol) NanoTail read_polya_multiple
#' https://github.com/LRB-IIMCB/nanotail/
#' #' Many thanks to the NanoTail developer for help and advice!
#'
#' @return a [tibble] containing non-A residue data for
#' all specified samples, with metadata provided in samples_table
#' stored as separate columns
#'
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' samples_table <- data.frame(residue_path = c(path_1,path_2,path_3,path_4),
#'                             sample_name =c("wt_1","mut_1","wt_2","mut_2"),
#'                             group = c("wt","mut","wt","mut"))
#'
#' residues_data <- ninetails::read_residue_multiple(samples_table)
#'
#'}
read_residue_multiple <- function(samples_table,...) {

  # var binding
  sample_name <- residue_path <- residue_contents <- NULL

  #assertions
  if (missing(samples_table)) {
    stop("The samples_table argument is missing",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(samples_table),
                          msg = "Empty data frame provided as an input (samples_table). Please provide samples_table describing data to load.")
  assertthat::assert_that("residue_path" %in% colnames(samples_table),
                          msg = "The samples_table should contain at least residue_path and sample_name columns.")
  assertthat::assert_that("sample_name" %in% colnames(samples_table),
                          msg = "The samples_table should contain at least residue_path and sample_name columns.")

  samples_data <- samples_table %>%
    dplyr::as.tbl() %>%
    dplyr::mutate_if(is.character,as.factor) %>%
    dplyr::mutate(residue_path = as.character(residue_path)) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(residue_contents=purrr::map(residue_path, function(x) ninetails::read_residue_single(x))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-residue_path)
  residue_data <- tidyr::unnest(samples_data,cols = c(residue_contents))

  return(residue_data)
}

#' Counts non-A residues found in nonadenosine_residues dataframe produced
#' by ninetails pipeline.
#'
#' Process the information returned by ninetails found in the nonadenosine_residues
#' dataframe in the prediction column. It is important to highlight that this function
#' returns counts as summary of hits, it does not provide a summary per read.
#'
#' @param residue_data A dataframe or tibble containig non-A residue predictions
#' made by ninetails pipeline
#'
#' @param grouping_factor grouping variable (e.g. "sample_name")
#'
#' @return A tibble with counts for each non-A residue
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' residue_counted <- ninetails::count_residues(residue_data=out[[2]],
#'                                              grouping_factor=NA)
#'}
#'
count_residues <- function(residue_data, grouping_factor=NA) {
  # variable binding
  prediction <- NULL

  # assertions
  assertthat::assert_that(assertive::has_rows(residue_data),
                          msg = "Empty dataframe provided as an input")

  if(!is.na(grouping_factor)) {
    assertthat::assert_that(grouping_factor %in% colnames(residue_data),
                            msg=paste0(grouping_factor," is not a column of input dataset"))

    residue_counts <- residue_data %>%
      dplyr::mutate(prediction=forcats::fct_relevel(prediction,"C", after = Inf)) %>%
      dplyr::group_by(!!rlang::sym(grouping_factor),prediction) %>%
      dplyr::count()
  }
  else {
    residue_counts <- residue_data %>%
      dplyr::mutate(prediction=forcats::fct_relevel(prediction,"C", after = Inf)) %>%
      dplyr::group_by(prediction) %>%
      dplyr::count()
  }

  return (residue_counts)
}


#' Spreads counts of the nonA residues in modified reads in
#' *_nonadenosine_residues dataframe (residue data) produced by ninetails pipeline
#' to 3 separate columns, each containing respective C, G, U prediction
#' per given read.
#'
#' In addition, an extra "nonA_residues" column is located at the end of
#' the output table. It contains all non-A residues positions summarized (per read),
#' given from the 5' to 3' end, separated by ":".
#'
#' @param residue_data A dataframe or tibble containig non-A residue predictions
#' made by ninetails pipeline (filename ends with "_nonadenosine_residues.tsv")
#'
#' @return a tibble with spread prediction column (3 cols instead)
#' @export
#'
#' @examples
#'\dontrun{
#'
#' spread_table <- ninetails::spread_nonA_residues(residue_data=residue_data)
#'
#'}
spread_nonA_residues <- function(residue_data){

  # var binding
  group <- readname <- est_nonA_pos <- prediction <- NULL

  #Assertions
  if (missing(residue_data)) {
    stop("A datafraeme with non-A residue data is missing. Please provide a valid residue_data argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(residue_data),
                          msg = "Empty dataframe provided as an input (residue_data)")


  #create column for non-A residue summary (cigar-like string)
  cigar <- residue_data %>% dplyr::group_by(group, readname) %>%
    dplyr::arrange(est_nonA_pos, .by_group = TRUE) %>%
    dplyr::summarise(nonA_residues = paste0(prediction, est_nonA_pos, collapse = ':'), .groups = 'drop')

  # create contingency table with C, G, U counts per read
  contingency <- residue_data %>%
    dplyr::select(-est_nonA_pos) %>%
    tidyr::pivot_wider(names_from = prediction,
                       names_sort = TRUE,
                       names_prefix = 'prediction_',
                       values_from = prediction,
                       values_fn = length,
                       values_fill = 0)

  #merge both tables
  spread_table <- contingency %>% dplyr::left_join(cigar, by=c("readname", "group"))


  return(spread_table)
}

#' Merges ninetails tabular outputs (read classes and nonadenosine residue data)
#' to produce one concise table for all data.
#'
#' Each read is represented by single row.
#' The output of this function is a tibble containing additional columns.
#' The "prediction" column from the residue_data is spread to 3 separate columns
#' (named with "prediction_" prefix), containing either C, G or U hits per read,
#' respectively. "Hit" is understood as a single presence of a given non-A residue.
#'
#' In addition, an extra "nonA_residues" column is located at the end of
#' the output table. It contains all non-A residues positions summarized (per read),
#' given from the 5' to 3' end, separated by ":".
#'
#' In this table, only reads that have been classified by ninetails
#' are included (reads marked "unclassified" are omitted from the analysis).
#'
#' @param class_data A dataframe or tibble containing read_classes predictions
#' made by ninetails pipeline
#'
#' @param residue_data A dataframe or tibble containig non-A residue predictions
#' made by ninetails pipeline
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#' IMPORTANT NOTE! This option must be set the same as it was during the
#' pipeline output production.
#'
#' @return a tibble with summarized info from both ninetails outputs
#' @export
#'
#' @examples
#'\dontrun{
#'
#' merged_tables <- ninetails::merge_nonA_tables(class_data=class_data,
#'                                               residue_data=residue_data,
#'                                               pass_only=TRUE)
#'
#'}
merge_nonA_tables <- function(class_data, residue_data, pass_only=TRUE){


  #Assertions
  if (missing(class_data)) {
    stop("A datafraeme with class data is missing. Please provide a valid class_data argument",
         call. = FALSE)
  }
  if (missing(residue_data)) {
    stop("A datafraeme with non-A residue data is missing. Please provide a valid residue_data argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(class_data),
                          msg = "Empty dataframe provided as an input (class_data)")
  assertthat::assert_that(assertive::has_rows(residue_data),
                          msg = "Empty dataframe provided as an input (residue_data)")
  assertthat::assert_that(assertive::is_a_bool(pass_only),
                          msg="Please provide TRUE/FALSE values for pass_only parameter")


  #drop all unclassified reads
  class_data <- class_data[!(class_data$class=="unclassified"),]

  # filter class_data according to predefined condition
  if(pass_only==TRUE){
    class2 <- class_data[class_data$qc_tag == "PASS", ]
  } else {
    class2 <- class_data[class_data$qc_tag %in% c("PASS", "SUFFCLIP"), ]
  }

  #spread residue_data
  spread_table <- ninetails::spread_nonA_residues(residue_data)

  #merge the data (spreaded residue + classes)
  merged_tables <- class2 %>% dplyr::full_join(spread_table)

  # replace NA with 0 in numeric columns
  # tidyselect::where() issue solved as in
  # https://stackoverflow.com/questions/62459736/how-do-i-use-tidyselect-where-in-a-custom-package
  merged_nonA_tables <- merged_tables %>%
    dplyr::mutate(dplyr::across(tidyselect::vars_select_helpers$where(is.numeric), ~ ifelse(is.na(.), 0, .)))

  return(merged_nonA_tables)

}

#' Produces summary table of nonA occurrences within analyzed dataset.
#'
#' Creates a table with statistics for non-A residues occurrences
#' for each transcript (contig) per each analyzed sample in given dataset.
#'
#' In the table, "counts" are understood as the number of reads in total
#' or containing a given type of non-A residue (see column headers for details).
#' Whereas "hits" are understood as the number of occurrences of a given
#' modification in total (see column headers for details). Please be aware that
#' there may be several "hits" in one read.
#'
#' The function also reports the mean and median poly(A) tail length
#' by transcript.
#'
#' @param merged_nonA_tables an output of merge_nonA_tables() function
#'
#' @param summary_factors character string or vector of strings;
#' column(s) used for grouping (default: "group")
#'
#' @param transcript_id_column character string; column with transcript id data
#' (default: "contig", as inherited from nanopolish; can be changed by the user)
#'
#' @return a summary table of nonA occurrences (tibble)
#' @export
#'
#' @examples
#'\dontrun{
#'
#' summarized <- ninetails::summarize_nonA(merged_nonA_tables=merged_nonA_tables,
#'                                         summary_factors="group",
#'                                         transcript_id_column="contig")
#'
#'}
summarize_nonA <- function(merged_nonA_tables,
                           summary_factors = c("group"),
                           transcript_id_column = c("contig")) {
  # var binding
  starts_with <- sum_nonA <- prediction_C <- prediction_G <- prediction_U <- median <- polya_length <- counts_total <- counts_unmod <- zeromod_summarized <- NULL

  #Assertions
  if (missing(merged_nonA_tables)) {
    stop("Ninetails' merged_nonA_tables output is missing. Please provide a valid merged_nonA_tables argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(merged_nonA_tables),
                          msg = "Empty dataframe provided as an input")
  assertthat::assert_that(assertive::is_character(summary_factors),
                          msg = "Non-character argument is not alowed for `summary_factors`. Please provide either string or vector of strings")
  assertthat::assert_that(all(summary_factors %in% colnames(merged_nonA_tables)),
                          msg="Non-existent column name provided as the argument (summary_factors)")


  # this function is slow; TODO: rethink the syntax
  # this syntax is slightly faster than without creating new variable according to microbenchmark
  nonA_data_summarized <- merged_nonA_tables %>%
    # drop all previous grouping vars
    dplyr::ungroup() %>%
    # count all nonA reads
    dplyr::mutate(sum_nonA = rowSums(dplyr::across(dplyr::starts_with('prediction_')))) %>%
    # group by provided vars
    dplyr::group_by(!!!rlang::syms(c(transcript_id_column,summary_factors))) %>%
    #provide summaries
    dplyr::summarise(
      polya_median = stats::median(polya_length), #add median polya length
      polya_mean = mean(polya_length), # add mean polya length
      counts_total = dplyr::n(), # add transcript count
      #summarize counts (number of reads with given non-A residues)
      #and hits (occurences of given residues)
      dplyr::across(c(sum_nonA, prediction_C, prediction_G, prediction_U),
                    list(counts = ~ sum(.x != 0), hits = ~ sum(.x))), .groups= 'drop') %>%
    # rename columns according to desired convention
    dplyr::rename_with(~stringr::str_replace(.x, '^\\w+_(\\w+)_(\\w+)', '\\2_\\1'), 4:dplyr::last_col())

  # add counts of unmodified reads - it has to be assigned to new var, because otherwise throws an error
  zeromod_summarized <- merged_nonA_tables %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!rlang::syms(c(transcript_id_column,summary_factors))) %>%
    dplyr::summarise(counts_unmod = sum(dplyr::if_all(tidyselect::starts_with('prediction_'), ~ .x == 0)), .groups = 'drop')

  nonA_data_summarized <- nonA_data_summarized%>%
    # add counts of unmodified reads - this syntax is because otherwise it throws an error
    # also: avoiding assigning new variable to temporary table
    dplyr::left_join(zeromod_summarized, by=c(transcript_id_column,summary_factors)) %>%
    # move unmodified reads count near total counts - for table clarity
    dplyr::relocate(counts_unmod, .after = counts_total)

  return(nonA_data_summarized)
}

#' Performs Fisher's exact test for testing the null of independence of rows and
#' columns in a contingency table representing given transcript in ninetails
#' output data. This is a wrapper for fisher.test function from stats package
#' with additional features to facilitate data wrangling.
#'
#' It is suitable only for the pairwise comparisons (i.e. for 2x2 contingency
#' table), where 2 conditions (e.g. WT vs KO) are compared at once.
#'
#' The function allows the user to set a cutoff number of reads required for
#' the analysis.
#'
#' This function is intended to work under the calculate_fisher function.
#'
#' The function was inspired by the Nanotail package written & maintained by
#' Pawel Krawczyk (smaegol): https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R
#'
#' Many thanks to the developer of original source code.
#'
#'
#' @param ninetails_data dataframe - the output of ninetails::merge_nonA_tables
#' function (merged tabular output containing read classification &
#' non-A position data).
#'
#' @param grouping_factor [character string] the name of factor variable defining
#' groups/conditions (needs to have 2 levels!)
#'
#' @param base [character string] letter representing particular non-A nucleotide,
#' for which the statistics are meant to be computed. Currently function accepts
#' C, G, U arguments. The "C" value is set by default.
#'
#' @param min_reads [numeric] minimum number of reads representing given
#' transcript to include it in the analysis
#'
#' @param transcript_id_column [character string] name of the column in which
#' the identifiers of the transcripts are stored. It is set to NA
#' by default.
#'
#' @return a tibble with results for given transcript, including pvalue,
#' adjusted pvalue, stats_code (the variable describing whether conditions
#' are met) and significance (FDR based on padj).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' test <- ninetails::nonA_fisher(ninetails_data=merged_nonA_tables,
#'                                grouping_factor = "sample_name",
#'                                base="C",
#'                                min_reads=100)
#'
#' }
#'
nonA_fisher <- function(ninetails_data,grouping_factor, base, min_reads=0, transcript_id_column=NA) {

  # var binding
  counts_unmod <- NULL

  # Assertions
  if (missing(ninetails_data)) {
    stop("Ninetails data are missing. Please provide a valid ninetails_data argument",
         call. = FALSE)
  }
  if (missing(base)) {
    stop("Base is missing. Please provide 'base' argument as character string (C, G or U).",
         call. =FALSE)
  }
  if (missing(transcript_id_column)) {
    stop("Transcript_id_column is missing. Please provide 'transcript_id_column' argument as character string.",
         call. =FALSE)
  }

  assertthat::assert_that(assertive::has_rows(ninetails_data),
                          msg = "Empty data.frame provided as an input")
  assertthat::assert_that(assertive::is_numeric(min_reads),
                          msg = "Non-numeric parameter provided (min_reads)")
  assertthat::assert_that(grouping_factor %in% colnames(ninetails_data),
                          msg=paste0(grouping_factor," is not a column of input dataset"))


  # if grouping factor has more than two levels
  if (length(levels(ninetails_data[[grouping_factor]]))>2) {
    if(is.na(condition1) && is.na(condition2)) {
      #throw error when no conditions for comparison are specified
      stop(paste0("grouping_factor ",grouping_factor," has more than 2 levels. Please specify condtion1 and condition2 to select comparison pairs"))
    } else {
      # filter input data leaving only specified conditions, dropping other factor levels
      assertthat::assert_that(condition1 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition1," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition2," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 != condition1,
                              msg="condition2 should be different than condition1")
      ninetails_data <- ninetails_data %>% dplyr::filter(!!rlang::sym(grouping_factor) %in% c(condition1,condition2)) %>% droplevels()
    }
  } else if (length(levels(ninetails_data[[grouping_factor]]))==1) {
    stop("Only 1 level present for grouping factor. Choose another groping factor for comparison")
  } else {
    condition1 = levels(ninetails_data[[grouping_factor]])[1]
    condition2 = levels(ninetails_data[[grouping_factor]])[2]
  }

  # initial status code
  stats_code = codes_stats = "OK"
  # calculate group counts
  group_counts = ninetails_data %>% dplyr::group_by(!!!rlang::syms(c(grouping_factor))) %>% dplyr::count()

  stats <- NA

  if (base=="C") {
    count_column <- "counts_C"
  } else if (base=="G") {
    count_column <- "counts_G"
  } else if (base=="U") {
    count_column <- "counts_U"
  } else {
    stop("Wrong non-A nucleotide defined. To compute statistics, please provide 'base' argument as character string (C, G or U).")
  }


  if (nrow(group_counts)==2) {
    if (group_counts[1,]$n < min_reads) {
      if (group_counts[2,]$n < min_reads) {
        stats_code = "B_LC"
      } else {
        stats_code = "G_LC"
      }
    } else if (group_counts[2,]$n < min_reads) {
      stats_code = "G_LC"
    } else {
      options(scipen = 999)

      # summarize nonAs
      contingency_table <- ninetails::summarize_nonA(merged_nonA_tables = ninetails_data,
                                                     summary_factors=grouping_factor,
                                                     transcript_id_column=transcript_id_column) %>%
        dplyr::select(!!rlang::sym(grouping_factor),
                      counts_unmod,
                      !!rlang::sym(count_column))
      contingency_table <- as.data.frame(contingency_table) # coerce tibble to df as setting names to tibble is deprecated
      row.names(contingency_table) <- contingency_table[[grouping_factor]] # set rownames
      contingency_table[[grouping_factor]] <- NULL # drop grouping col
      stats <- suppressWarnings(stats::fisher.test(contingency_table))$p.value

    }
  } else if (nrow(group_counts)==1) {
    stats_code = "G_NA"
  } else if (nrow(group_counts)==0) {
    stats_code = "B_NA"
  } else {
    stats_code = "ERR"
  }

  # create output
  stats <- tibble::tibble(p.value=stats,stats_code=as.character(stats_code))

  return(stats)

}

stat_codes_list = list(OK = "OK",
                       G1_NA = "GROUP1_NA",
                       G2_NA = "GROUP2_NA",
                       G1_LC = "G1_LOW_COUNT",
                       G2_LC = "G2_LOW_COUNT",
                       B_NA = "DATA FOR BOTH GROUPS NOT AVAILABLE",
                       B_LC = "LOW COUNTS FOR BOTH GROUPS",
                       G_LC = "LOW COUNT FOR ONE GROUP",
                       G_NA = "DATA FOR ONE GROUP NOT AVAILABLE",
                       ERR = "OTHER ERROR")


#' Performs Fisher's exact test for each transcript in ninetails
#' output data. Then, it performs the Benjamini-Hochberg procedure
#' (BH step-up procedure) to control the FDR.
#'
#' This is a wrapper for fisher.test function from stats package and p.adjust
#' functions with additional features to facilitate data wrangling.
#'
#' The function was inspired by the Nanotail package written & maintained by
#' Pawel Krawczyk (smaegol): https://github.com/LRB-IIMCB/nanotail/blob/dev/R/polya_stats.R
#'
#' Many thanks to the developer of original source code.
#'
#' @param ninetails_data dataframe - the output of ninetails::merge_nonA_tables
#' function (merged tabular output containing read classification &
#' non-A position data).
#'
#' @param transcript_id_column [character string] column with transcript id data
#' (default: "contig", as inherited from nanopolish; can be changed by the user)
#'
#' @param min_reads [numeric] minimum number of reads representing given
#' transcript to include it in the analysis. This parameter is set by default
#' to 0. Please keep in mind that taking into account many transcripts with low
#' coverage increases the risk of reject true null hypothesis
#' (Benjamini-Hochberg procedure).
#'
#' @param min_nonA_reads [numeric] minimum number of reads containing
#' nonadenosine residues (summary for C, G, U alltogether) per given
#' transcript to include it in the analysis. This parameter prevents from
#' considering too many observations as nonsignificant in the further pvalue
#' adjustation. In general, the non-A containing reads are a small fraction
#' of the total pool of reads. As a rule of thumb, additional filtering
#' criteria can provide more valuable information regarding the samples
#' (prevent from rejecting true null hypothesis).
#' This is set by default to 0.
#'
#' @param grouping_factor [character string] grouping variable
#' (e.g. "sample_name" - default)
#'
#' @param condition1 [character string] first level of `grouping_factor`
#' to use for comparison
#'
#' @param condition2 [character string] second level of `grouping_factor`
#' to use for comparison
#'
#' @param alpha [numeric] an alpha value to consider a hit significant.
#' Default: 0.05.
#'
#' @param base [character string] letter representing particular non-A nucleotide,
#' for which the statistics are meant to be computed. Currently function accepts
#' C, G, U arguments. The "C" value is set by default.
#'
#' @param ... additional parameters to pass to nonA_fisher (under development)
#'
#' @return a summary table with pvalues, padj and significance levels
#' for each transcript  (tibble)
#'
#' @importFrom rlang :=
#'
#' @export
#'
#' @examples
#' \dontrun{
#' test <- ninetails::calculate_fisher(ninetails_data=merged_nonA_tables,
#'                                     transcript_id_column = "contig",
#'                                     min_reads=100,
#'                                     min_nonA_reads=10,
#'                                     grouping_factor = "sample_name",
#'                                     condition1="WT",
#'                                     condition2="KO",
#'                                     alpha=0.05,
#'                                     base="C")
#' }
calculate_fisher <- function(ninetails_data,
                             transcript_id_column = "contig",
                             min_reads = 0,
                             min_nonA_reads=0,
                             grouping_factor = "sample_name",
                             condition1=NA,
                             condition2=NA,
                             alpha=0.05,
                             base="C",
                             ...)
{

  # vr binding
  merged_nonA_tables<- sum_nonA <- counts_nonA <- contig <- data<- stats<- padj<- NULL

  # Assertions
  if (missing(ninetails_data)) {
    stop("Ninetails data are missing. Please provide a valid ninetails_data argument",
         call. = FALSE)
  }
  if (missing(transcript_id_column)) {
    stop("Transcript_id_column is missing. Please provide a valid transcript_id_column argument",
         call. = FALSE)
  }
  if (missing(min_reads)) {
    stop("Min_reads are missing. Please provide a valid min_reads argument",
         call. = FALSE)
  }
  if (missing(min_nonA_reads)) {
    stop("Min_nonA_reads data are missing. Please provide a valid min_nonA_reads argument",
         call. = FALSE)
  }
  if (missing(base)) {
    stop("Base definition is missing. Please provide a valid base argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(min_reads),
                          msg=paste0("Min_reads must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_numeric(min_nonA_reads),
                          msg=paste0("Min_nonA_reads must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_numeric(alpha),
                          msg=paste0("Alpha must be numeric. Please provide a valid argument."))

  assertthat::assert_that(assertive::has_rows(ninetails_data),
                          msg = "Empty data.frame provided as an input")


  # if grouping factor has more than two levels
  if (length(levels(ninetails_data[[grouping_factor]]))>2) {
    if(is.na(condition1) && is.na(condition2)) {
      #throw error when no conditions for comparison are specified
      stop(paste0("grouping_factor ",grouping_factor," has more than 2 levels. Please specify condtion1 and condition2 to select comparison pairs"))
    }
    else {
      # filter input data leaving only specified conditions, dropping other factor levels
      assertthat::assert_that(condition1 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition1," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 %in% levels(ninetails_data[[grouping_factor]]),
                              msg=paste0(condition2," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")

      ninetails_data <- ninetails_data %>% dplyr::filter(!!rlang::sym(grouping_factor) %in% c(condition1,condition2)) %>%
        dplyr::mutate() %>%
        droplevels()

    }
  }
  else if (length(levels(ninetails_data[[grouping_factor]]))==1) {
    stop("Only 1 level present for grouping factor. Choose another groping factor for comparison")
  }
  else {
    condition1 = levels(ninetails_data[[grouping_factor]])[1]
    condition2 = levels(ninetails_data[[grouping_factor]])[2]
  }

  # filter out transcripts with not enough amount of non-A reads among the whole pool of reads:
  mod_summarized <- ninetails_data %>% dplyr::ungroup() %>%
    dplyr::mutate(sum_nonA = rowSums(dplyr::across(dplyr::starts_with('prediction_')))) %>%
    dplyr::group_by(!!!rlang::syms(c(transcript_id_column,grouping_factor))) %>%
    dplyr::summarise(dplyr::across(c(sum_nonA), list(counts = ~ sum(.x != 0))), .groups= 'drop') %>%
    dplyr::rename_with(~stringr::str_replace(.x, '^\\w+_(\\w+)_(\\w+)', '\\2_\\1'), 3:dplyr::last_col())
  # apply filtering criterion (minimal nonA read content)
  mod_summarized_filtered <- mod_summarized %>% dplyr::filter(counts_nonA>=min_nonA_reads)
  #extract filtered readnames
  mod_summarized_filtered <- unique(mod_summarized_filtered$contig)

  ninetails_data <- ninetails_data[ninetails_data$contig %in% mod_summarized_filtered,]


  ninetails_data_stat <- ninetails_data %>%
    dplyr::mutate(transcript_id=get(c(transcript_id_column))) %>%
    dplyr::group_by(get(c(transcript_id_column))) %>%
    tidyr::nest()


  ninetails_data_stat <- ninetails_data_stat %>%
    dplyr::mutate(stats=purrr::map(data,ninetails::nonA_fisher,grouping_factor=grouping_factor,min_reads=min_reads, base=base,transcript_id_column = "transcript_id")) %>%
    dplyr::select(-data) %>% tidyr::unnest(cols = c(stats)) %>% dplyr::rename(!!transcript_id_column := "get(c(transcript_id_column))")
  message("calculating statistics")

  message("Finished")

  message("Adjusting p.value")
  ninetails_data_stat$padj <- stats::p.adjust(ninetails_data_stat$p.value, method = "BH")

  # create significance factor
  ninetails_data_stat <- ninetails_data_stat %>%
    dplyr::mutate(significance = dplyr::case_when(is.na(padj)  ~ "NotSig",
                                                  (padj < alpha) ~ paste0("FDR<", alpha),
                                                  TRUE ~ "NotSig"))
  ninetails_data_stat$stats_code <- sapply(ninetails_data_stat$stats_code,
                                           FUN = function(x) {stat_codes_list[[x]]},simplify = "vector",USE.NAMES = FALSE) %>% unlist()



  ninetails_data_stat <- ninetails_data_stat %>% dplyr::arrange(padj)

  return(ninetails_data_stat)
}


