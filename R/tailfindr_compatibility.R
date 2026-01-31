#' Converts tailfindR results to format compatible with Ninetails
#'
#' This function allows to reformat the output of the tailfindR pipeline,
#' so it could be passed to the Ninetails pipeline.
#'
#' The function renames some crucial columns according to the naming convention
#' adopted by Ninetails:\itemize{
#'  \item read_id -> readname
#'  \item tail_start -> polya_start
#'  \item tail_length -> polya_length
#' }
#'
#' The function creates additional columns based on columns present in original
#' tailfindR output, to mimick the output of nanopolish polya function:\itemize{
#'  \item transcript_start - based on the coordinate in tail_end +1 position forward
#'  \item contig - nanopolish returns mapping info, which is absent in tailfindR
#'  output. Therefore, the dummy value "tailfindr_out" is inserted there.
#'  \item qc_tag - since tailfindR does not provide the quality tag, the simple
#'  dummy is introduced. As the detection region of the R9.4.1 pore spans 5
#'  nucleotides, which in turn leads to predictions shorter than 10 nt (2 full
#'  adjacent 5-mers with no overlaps) to be unreliable, this variable assumes 2
#'  possible values: PASS for tails longer or equal to 10 nt, and NOREGION
#'  for tails shorter than 10 nt. This is not the ideal solution, but the most
#'  plausible one, given the fact that tailfindR does not produce the signal
#'  quality metrics itself.
#' }
#'
#' IMPORTANT WARNING
#'
#' Ninetails was originally intended to work with the output produced by
#' nanopolish polya function.
#' According to our experience in Laboratory of RNA Biology, HMM provide more
#' robust prediction than the slope estimator. Therefore, we rely on the output
#' provided by nanopolish. Ninetails is optimized to work with data produced by
#' nanopolish. Hence, predictions made with tailfindR results should be treated
#' with high caution. Especially given that the Ninetails relies heavily on
#' quality tags produced by nanopolish polya function. It takes into consideration
#' only those signals, which fulfill nanopolish quality criteria (marked either
#' as "PASS" or "SUFFCLIP"). In case of tailfindR output, the quality of signal
#' can't be inferred based on the result provided. Therefore it may happen that
#' the signal of poor quality, insufficient to produce reliable predictions,
#' can be passed to the neural network, and - as a result - misclassified.
#'
#' TL;DR
#' Ninetails is optimized to work with nanopolish. If you use it with tailfindR
#' output - it is on your own risk.
#'
#' @param tailfindr_output character string. Either full path of the .csv file
#' produced by tailfindR or an environment object containing tailfindR result
#' data.
#'
#' @return A tibble/df containing tailfindR results reformatted to resemble the
#' output of nanopolish polya function compatible with the Ninetails pipeline.
#' Always assign this returned tibble to a variable. Otherwise
#' the long tibble will be printed to the console, which may crash your
#' R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' df <- ninetails::convert_tailfindr_output(tailfindr_output = '/tf/out.csv')
#'
#' #The output of this function shall be passed to the Ninetails pipeline as
#' nanopolish variable, for instance:
#'
#' results <- ninetails::check_tails(nanopolish = converted_tailfindr,
#'   sequencing_summary = '/path/to/sequencing_summary.txt',
#'   workspace = '/path/to/workspace',
#'   num_cores = 2,
#'   basecall_group = 'Basecall_1D_000',
#'   pass_only=TRUE,
#'   save_dir = '~/Downloads')
#'}
convert_tailfindr_output <- function(tailfindr_output){

  #var binding
  read_id <- tail_start <- tail_end <- tail_length <- polya_length <- NULL

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
#' This function checks whether the input poly(A) length files store data from
#' correct sequencing platform (Oxford Nanopore DRS sequencing) - if other data
#' types are loaded, an error occurs.
#'
#' Please keep in mind:
#' This function is an element of legacy pipeline which would probably
#' be retired if fast5 format would be deprecated.
#'
#' @param input Path to file or data frame containing poly-A length data:
#' either produced by nanopolish or tailfindr. In case of tailfindr, only DRS
#' data are accepted.
#'
#' @return Standardized data frame in nanopolish-like format compatible with
#' further ninetails processing steps
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' # run on nanopolish output
#' test <- ninetails::check_polya_length_filetype(
#' input = system.file('extdata',
#'                     'test_data',
#'                     'nanopolish_output.tsv',
#'                     package = 'ninetails'))
#'
#' # run on tailfindr output
#' test <- ninetails::check_polya_length_filetype(
#' input = system.file('extdata',
#'                     'test_data',
#'                     'tailfindr_output.csv',
#'                     package = 'ninetails'))
#'
#'}
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
