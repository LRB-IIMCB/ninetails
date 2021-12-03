#' Creating parameters for nnet.
#'
#' @param nanopolish character string. Full path of the .tsv file produced
#' by nanopolish polya function.
#'
#' @param sequencing_summary character string. Full path of the .txt file
#' with sequencing summary.
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#'
#' @param pass_only logical [TRUE/FALSE]. If TRUE, only reads tagged by
#' nanopolish as "PASS" would be taken into consideration. Otherwise, reads
#' tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
#' As a default, "TRUE" value is set.
#'
#' @return A list containing read information organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @examples
#' \dontrun{
#'
#' check_signals(nanopolish = '/path/to/file',
#'               sequencing_summary = '/path/to/file',
#'               workspace = '/path/to/guppy/workspace',
#'               num_cores = 10,
#'               basecalled_group = 'Basecall_1D_000',
#'               pass_only=TRUE)
#'
#' }


check_signals <- function(nanopolish, sequencing_summary, workspace, num_cores, basecall_group, pass_only=TRUE){

  cat('Welcome to Ninetails v.0.0.4 (beta)\n',
      'Pipeline initialized:', as.character(Sys.time()),'\n','\n')

  feature_list <- create_feature_list(nanopolish, sequencing_summary, workspace, num_cores, basecall_group, pass_only=TRUE)
  tail_chunk_list <- create_chunk_list(feature_list, num_cores)
  #gasf_list <- create_gasf_list(tail_chunk_list, num_cores)

  cat('Processing finished. ***** ***')

  return(tail_chunk_list)
}
