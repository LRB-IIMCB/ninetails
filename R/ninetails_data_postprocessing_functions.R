#' Reads ninetails read_classes dataframe from file.
#'
#' This is the basic function used to import read_classes output from check_tails
#' to R.
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

  ## correct annotation <- if gencode format: create new columns with ensembl IDs
  # else created columns could be easily dropped
  # this code chunk was originally written by Paweł Krawczyk (smaegol) & incorporated in NanoTail package
  transcript_names <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", class_data$contig)
  class_data$transcript <- transcript_names
  ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", class_data$contig)
  ensembl_transcript_ids_short <- gsub("(.*)\\..*", "\\1", ensembl_transcript_ids) # without version number
  class_data$ensembl_transcript_id_full <- ensembl_transcript_ids
  class_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short


  return(class_data)
}



#' Reads multiple ninetails read_classes outputs at once.
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
#' @param ... - additional parameters to pass to \code{\link{read_class_single}};
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
    tibble::as_tibble() %>%
    dplyr::mutate_if(is.character,as.factor) %>%
    dplyr::mutate(class_path = as.character(class_path)) %>%
    dplyr::group_by(sample_name) %>% dplyr::mutate(class_contents=purrr::map(class_path, function(x) ninetails::read_class_single(x))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-class_path)

  class_data <- tidyr::unnest(samples_data, cols = c(class_contents))

  return(class_data)
}

#' Counts read classes found in read_classes dataframe produced
#' by ninetails pipeline.
#'
#' Process the information returned by ninetails found in the read_classes
#' dataframe in the prediction column.
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
        dplyr::mutate(comments=forcats::fct_relevel(comments,"YAY", after = Inf)) %>%
        dplyr::group_by(!!rlang::sym(grouping_factor),comments) %>%
        dplyr::count()
    } else {
      class_counts <- class_data %>%
        dplyr::mutate(comments=forcats::fct_relevel(comments,"YAY", after = Inf)) %>%
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


#' Reads ninetails nonadenosine_residues data from file.
#'
#' This is the basic function used to import output from \code{\link{check_tails}} to R:
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

  ## correct annotation <- if gencode format: create new columns with ensembl IDs
  # else created columns could be easily dropped
  # this code chunk was originally written by Paweł Krawczyk (smaegol) & incorporated in NanoTail package

  transcript_names <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", residue_data$contig)
  residue_data$transcript <- transcript_names
  ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", residue_data$contig)
  ensembl_transcript_ids_short <- gsub("(.*)\\..*", "\\1", ensembl_transcript_ids) # without version number
  residue_data$ensembl_transcript_id_full <- ensembl_transcript_ids
  residue_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short



  return(residue_data)
}


#' Reads multiple ninetails nonadenosine_residues outputs at once.
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
#' @param ... - additional parameters to pass to \code{\link{read_residue_single}};
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
    tibble::as_tibble() %>%
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
      dplyr::mutate(prediction=forcats::fct_relevel(prediction,"U", after = Inf)) %>%
      dplyr::group_by(!!rlang::sym(grouping_factor),prediction) %>%
      dplyr::count()
  }
  else {
    residue_counts <- residue_data %>%
      dplyr::mutate(prediction=forcats::fct_relevel(prediction,"U", after = Inf)) %>%
      dplyr::group_by(prediction) %>%
      dplyr::count()
  }

  return (residue_counts)
}


#' Reshapes nonadenosine_residues dataframe.
#'
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
#' for each transcript (ensembl_transcrit_id_short) per each analyzed sample in given dataset.
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
#' @param merged_nonA_tables an output of \code{\link{merge_nonA_tables}} function
#'
#' @param summary_factors character string or vector of strings;
#' column(s) used for grouping (default: "group")
#'
#' @param transcript_id_column character string; column with transcript id data
#' (default: "ensembl_transcript_id_short" as added during data preprocessing;
#' can be changed by the user)
#'
#' @return a summary table of nonA occurrences (tibble)
#' @export
#'
#' @examples
#'\dontrun{
#'
#' summarized <- ninetails::summarize_nonA(merged_nonA_tables=merged_nonA_tables,
#'                                         summary_factors="group",
#'                                         transcript_id_column="ensembl_transcript_id_short")
#'
#'}
summarize_nonA <- function(merged_nonA_tables,
                           summary_factors = c("group"),
                           transcript_id_column = c("ensembl_transcript_id_short")) {
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


#' Aggregates the quality control info produced by nanopolish polya function.
#'
#' This function returns read counts assigned to each of the qc_tags
#' per user-predefined grouping variable (grouping_factor) which might be either
#' a sample name or an experiment condition (the column of choice must be present
#' within the input table).
#'
#' This is the ninetails' implementation of
#' \code{\link[nanotail:get_nanopolish_processing_info]{name}}
#' function originally written by P. Krawczyk (smeagol) and incorporated within
#' the NanoTail package.
#'
#' For original source code, see:
#' https://github.com/LRB-IIMCB/nanotail/blob/master/R/polya_stats.R
#'
#' The variable names were adjusted according to the naming convention within
#' ninetails to avoid confusion.
#'
#' Many thanks to the author of original source code for help and advice.
#'
#' @param class_data A dataframe or tibble containig class_data output
#' from Ninetails
#'
#' CAUTION! Do not use \code{\link{merge_nonA_tables}} output, since it is cleaned from
#' poor quality reads by default!
#'
#' @param grouping_factor [character string] variable used for grouping the data
#' (e.g. by sample_name)
#'
#' @return A tibble with counts for each qc_tag present in the run
#' @export
#'
nanopolish_qc <- function(class_data,
                          grouping_factor=NA) {
  #var binding
  qc_tag <- NULL

  #assertions
  assertthat::assert_that(assertive::has_rows(class_data),
                          msg = "Empty dataframe provided as an input")


  if(!is.na(grouping_factor)) {
    assertthat::assert_that(grouping_factor %in% colnames(class_data),
                            msg=paste0(grouping_factor," is not a column of input dataset"))
    processing_info <- class_data %>%
      dplyr::mutate(qc_tag=forcats::fct_relevel(qc_tag,"PASS", after = Inf)) %>%
      dplyr::group_by(!!rlang::sym(grouping_factor),qc_tag) %>%
      dplyr::count()
  }
  else {
    processing_info <- class_data %>%
      dplyr::mutate(qc_tag=forcats::fct_relevel(qc_tag,"PASS", after = Inf)) %>%
      dplyr::group_by(qc_tag) %>%
      dplyr::count()
  }

  return (processing_info)
}

