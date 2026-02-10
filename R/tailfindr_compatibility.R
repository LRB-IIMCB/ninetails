#' Converts tailfindr results to format compatible with ninetails
#'
#' Reformats the output of the tailfindr pipeline so it can be passed
#' to the ninetails legacy (Guppy) pipeline in place of nanopolish polya
#' output. The function renames key columns, creates derived columns, and
#' introduces dummy quality tags to approximate the nanopolish output
#' schema.
#'
#' @details
#' The following column mappings are applied:
#' \describe{
#'   \item{\code{read_id} \eqn{\rightarrow} \code{readname}}{Read identifier
#'     renamed to match the ninetails naming convention.}
#'   \item{\code{tail_start} \eqn{\rightarrow} \code{polya_start}}{Start
#'     coordinate of the poly(A) tail.}
#'   \item{\code{tail_length} \eqn{\rightarrow} \code{polya_length}}{Estimated
#'     poly(A) tail length.}
#' }
#'
#' The following columns are created:
#' \describe{
#'   \item{\code{transcript_start}}{Set to \code{tail_end + 1}.}
#'   \item{\code{contig}}{Dummy value \code{"tailfindr_out"} (nanopolish
#'     returns mapping information, which is absent in tailfindr output).}
#'   \item{\code{qc_tag}}{Simplified quality tag. Because the R9.4.1 pore
#'     detection region spans 5 nucleotides, tails shorter than 10 nt
#'     (2 full adjacent 5-mers with no overlap) are unreliable. Therefore
#'     \code{qc_tag} is set to \code{"PASS"} for tails >= 10 nt and
#'     \code{"NOREGION"} otherwise.}
#' }
#'
#' @section Warning:
#' Ninetails is optimised to work with nanopolish polya output. The
#' HMM-based approach of nanopolish provides more robust tail boundary
#' predictions than the slope estimator used by tailfindr. In particular,
#' ninetails relies on the quality tags produced by nanopolish
#' (\code{"PASS"} / \code{"SUFFCLIP"}). When using tailfindr output, the
#' signal quality cannot be inferred, which may lead to poor-quality
#' signals being passed to the CNN and consequently misclassified.
#' \strong{Use tailfindr input at your own risk.}
#'
#' @param tailfindr_output Character string or data frame. Either the full
#'   path of the \code{.csv} file produced by tailfindr, or an in-memory
#'   data frame containing tailfindr result data. Required columns:
#'   \code{read_id}, \code{tail_start}, \code{tail_end},
#'   \code{tail_length}.
#'
#' @return A tibble containing tailfindr results reformatted to resemble
#'   the output of \code{nanopolish polya}, with columns:
#' \describe{
#'   \item{readname}{Character. Read identifier (renamed from
#'     \code{read_id}).}
#'   \item{polya_start}{Integer. Start coordinate of the poly(A) tail
#'     (renamed from \code{tail_start}).}
#'   \item{tail_end}{Integer. End coordinate of the poly(A) tail
#'     (retained from tailfindr).}
#'   \item{polya_length}{Numeric. Estimated tail length (renamed from
#'     \code{tail_length}).}
#'   \item{transcript_start}{Integer. Transcript start coordinate
#'     (\code{tail_end + 1}).}
#'   \item{contig}{Character. Dummy mapping target
#'     (\code{"tailfindr_out"}).}
#'   \item{qc_tag}{Character. Quality tag (\code{"PASS"} or
#'     \code{"NOREGION"}).}
#' }
#' Always assign this returned tibble to a variable; printing the full
#' tibble to the console may crash the R session.
#'
#' @seealso \code{\link{check_polya_length_filetype}} which uses this
#'   function internally,
#'   \code{\link{check_tails_guppy}} for the legacy pipeline that accepts the
#'   converted output as its \code{nanopolish} argument,
#'   \code{\link{extract_polya_data}} for how nanopolish output is
#'   normally processed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' df <- ninetails::convert_tailfindr_output(
#'   tailfindr_output = '/path/to/tailfindr_out.csv')
#'
#' # The output can be passed to the ninetails pipeline as the
#' # nanopolish argument:
#' results <- ninetails::check_tails(
#'   nanopolish = df,
#'   sequencing_summary = '/path/to/sequencing_summary.txt',
#'   workspace = '/path/to/workspace',
#'   num_cores = 2,
#'   basecall_group = 'Basecall_1D_000',
#'   pass_only = TRUE,
#'   save_dir = '~/Downloads')
#'
#' }
convert_tailfindr_output <- function(tailfindr_output){


  # Accept either path to file or in-memory file
  if (is_string(tailfindr_output)) {
    # if string provided as an argument, read from file
    assert_file_exists(tailfindr_output, "tailfindr_output")
    tailfindr_out <- vroom::vroom(tailfindr_output,
                                  col_select=c(read_id, tail_start, tail_end, tail_length),
                                  show_col_types = FALSE)
  } else {
    # make sure that tailfindr_output is an object with rows
    if (!is.data.frame(tailfindr_output) || nrow(tailfindr_output) == 0) {
      stop("Empty data frame provided as an input (tailfindr_output). Please provide valid input")
    }

    tailfindr_out <- tailfindr_output[,c("read_id","tail_start","tail_end","tail_length")]
  }

  # rename columns to achieve compatibility with ninetails' naming convention
  # adjust transcript start coordinate
  # provide dummy variables for refseq (nanopolish returns mapping info)
  # and quality tag (tailfindr does not provide quality assessment)

  converted_tailfindr_output <- tailfindr_out %>%
    dplyr::rename(readname = read_id,
                  polya_start = tail_start,
                  polya_length = tail_length) %>%
    dplyr::mutate(transcript_start = tail_end+1,
                  contig = "tailfindr_out",
                  qc_tag = ifelse(polya_length >= 10, "PASS", "NOREGION"))


  return(converted_tailfindr_output)

}


#' Check and convert poly(A) length file format
#'
#' Detects whether an input poly(A) length file originates from nanopolish
#' or tailfindr and, if necessary, converts it to the nanopolish-like
#' format expected by the ninetails legacy (Guppy/Fast5) pipeline.
#' Tailfindr cDNA output is explicitly rejected because the legacy
#' pipeline supports only DRS (direct RNA sequencing) data.
#'
#' @details
#' File-type detection is based on column names:
#' \describe{
#'   \item{nanopolish}{Identified by the presence of a \code{qc_tag}
#'     column. Returned as-is.}
#'   \item{tailfindr DRS}{Identified by the presence of \code{read_id},
#'     \code{tail_start}, \code{tail_end}, \code{samples_per_nt}, and
#'     \code{tail_length}. Converted via
#'     \code{\link{convert_tailfindr_output}}.}
#'   \item{tailfindr cDNA}{Identified by the presence of a
#'     \code{tail_is_valid} column. Raises an error because cDNA data are
#'     not compatible with this pipeline.}
#' }
#' If none of the above patterns match, an error is raised.
#'
#' This function is part of the legacy pipeline (Guppy basecaller,
#' multi-Fast5 files). It may be retired if the Fast5 format is
#' deprecated.
#'
#' @param input Character string or data frame. Either the path to a
#'   poly(A) length file (nanopolish \code{.tsv} or tailfindr \code{.csv})
#'   or an in-memory data frame. Only DRS data are accepted for tailfindr
#'   input.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{data}{Data frame / tibble. Standardised poly(A) length data in
#'     nanopolish-like format, ready for downstream ninetails processing.}
#'   \item{file_type}{Character. One of \code{"nanopolish"} or
#'     \code{"tailfindr_drs"}, indicating the detected input format.}
#' }
#'
#' @seealso \code{\link{convert_tailfindr_output}} for the tailfindr
#'   conversion logic,
#'   \code{\link{extract_polya_data}} for subsequent processing of the
#'   standardised output,
#'   \code{\link{check_tails_guppy}} for the legacy pipeline entry point.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # run on nanopolish output
#' test <- ninetails::check_polya_length_filetype(
#'   input = system.file('extdata', 'test_data',
#'                       'nanopolish_output.tsv',
#'                       package = 'ninetails'))
#'
#' # run on tailfindr output
#' test <- ninetails::check_polya_length_filetype(
#'   input = system.file('extdata', 'test_data',
#'                       'tailfindr_output.csv',
#'                       package = 'ninetails'))
#'
#' }
check_polya_length_filetype <- function(input) {

  # Initialize result list
  result <- list(
    data = NULL,
    file_type = NULL
  )

  # Accept either path to file or in-memory file
  if (is_string(input)) {
    assert_file_exists(input, "Input")
    input_data <- vroom::vroom(input, show_col_types = FALSE)
  } else {
    if (!is.data.frame(input) || nrow(input) == 0) {
      stop("Empty data frame provided as input. Please provide valid input")
    }
    input_data <- input
  }

  # Get column names
  cols <- colnames(input_data)

  # Check file type based on columns
  is_nanopolish <- "qc_tag" %in% cols
  is_tailfindr_drs <- all(c("read_id", "tail_start", "tail_end",
                            "samples_per_nt", "tail_length") %in% cols)
  is_tailfindr_cdna <- "tail_is_valid" %in% cols

  # Process based on file type
  if (is_nanopolish) {
    result$data <- input_data
    result$file_type <- "nanopolish"
  } else if (is_tailfindr_drs) {
    result$data <- ninetails::convert_tailfindr_output(input_data)
    result$file_type <- "tailfindr_drs"
  } else if (is_tailfindr_cdna) {
    stop("Tailfindr cDNA output is not compatible with this pipeline.
          Please provide DRS (direct RNA sequencing) data.", call. = FALSE)
  } else {
    stop("Unrecognized input format. This pipeline accepts only nanopolish poly(A)
          or tailfindr DRS output. Please check your input file.", call. = FALSE)
  }

  return(result)
}
