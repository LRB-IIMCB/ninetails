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
#' @examples
#'\dontrun{
#'
#' converted_tailfindr <- ninetails::convert_tailfindr_output(tailfindr_output = '/path/to/tailfindr/output.csv')
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
  if (checkmate::test_string(tailfindr_output)) {
    # if string provided as an argument, read from file
    checkmate::assert_file_exists(tailfindr_output)
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
