#' Extracts features of multiple RNA reads from nanopolish output
#' & sequencing summary.
#'
#' This function extracts features of poly(A) tails of selected RNA reads
#' from the output table provided by nanopolish polya function and sequencing
#' summary provided by the sequencer. Filenames are taken from the sequencing
#' summary file.
#'
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function.
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @return A tibble/df containing read information organized by the read ID
#' is returned. Always assign this returned tibble to a variable. Otherwise
#' the long tibble will be printed to the console, which may crash your
#' R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' extract_polya_data(nanopolish = '/path/to/nanopolish/polya/output.tsv',
#'                    sequencing_summary = '/path/to/sequencing_summary.tsv',
#'                    pass_only = TRUE)
#'
#'}
#'
#'

extract_polya_data <- function(nanopolish, sequencing_summary, pass_only = TRUE){

  if (missing(nanopolish)) {
    stop("Nanopolish polya output is missing. Please provide a valid nanopolish argument.", .call = FALSE)
  }

  if (missing(sequencing_summary)) {
    stop("Sequencing summary file is missing. Please provide a valid sequencing_summary argument.", .call = FALSE)
  }


  nanopolish_polya_table <- vroom::vroom(nanopolish, col_select=c(readname, polya_start, transcript_start, adapter_start, leader_start, qc_tag), show_col_types = FALSE)
  sequencing_summary_table <- vroom::vroom(sequencing_summary, col_select = c(filename, read_id), show_col_types = FALSE)
  colnames(sequencing_summary_table)[2] <- "readname"

  #assertions
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(nanopolish),msg = "Empty string provided as an input. Please provide a nanopolish as a string")
  assertthat::assert_that(assertive::is_existing_file(nanopolish), msg=paste("File ",nanopolish," does not exist",sep=""))
  assertthat::assert_that(assertive::has_rows(nanopolish_polya_table), msg = "Empty data frame provided as an input (nanopolish). Please provide valid input")

  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(sequencing_summary),msg = "Empty string provided as an input. Please provide a sequencing_summary as a string")
  assertthat::assert_that(assertive::is_existing_file(sequencing_summary), msg=paste("File ",sequencing_summary," does not exist",sep=""))
  assertthat::assert_that(assertive::has_rows(sequencing_summary_table), msg = "Empty data frame provided as an input (sequencing_summary). Please provide valid input")

  assertthat::assert_that(assertive::is_a_bool(pass_only),msg="Please provide TRUE/FALSE values for pass_only parameter")


  # Add filtering criterion: select only pass or pass $ suffclip
  if(pass_only == TRUE){
    polya_summary <- dplyr::left_join(nanopolish_polya_table[which(nanopolish_polya_table$qc_tag=="PASS"),], sequencing_summary_table, by="readname")
  } else {
    polya_summary <- dplyr::left_join(nanopolish_polya_table[which(nanopolish_polya_table$qc_tag %in% c("PASS", "SUFFCLIP")),], sequencing_summary_table, by="readname")
  }

  names(polya_summary$filename) <- polya_summary$readname #named vec of filenames (names = readnames)
  attr(polya_summary, 'spec') <- NULL #drop attributes left by vroom

  #take only first n reads (n= value defined in "reads" argument)
  #polya_summary <- polya_summary[1:reads, ]

  return(polya_summary)
}
