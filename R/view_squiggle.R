#' Creates segmented plot of given raw ONT RNA signal.
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
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
#' @param basecall_group character string ("Basecall_1D_000" is set
#' as a default). Name of the level in the Fast5 file hierarchy from
#' which the data should be extracted.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#'
#' @return
#' @export
#'
#' @examples
#''\dontrun{
#'
#' view_squiggle(nanopolish = '/path/to/file',
#'               sequencing_summary = '/path/to/file',
#'               workspace = '/path/to/guppy/workspace',
#'               basecalled_group = 'Basecall_1D_000',
#'               num_cores = 3)
#'
#'}


view_squiggle <- function(readname, nanopolish, sequencing_summary, workspace, basecall_group = "Basecall_1D_000", num_cores){

  #Assertions
  if (missing(readname)) {
    stop("Readname [string] is missing. Please provide a valid readname argument.", call. =FALSE)
  }

  if (missing(workspace)) {
    stop("Directory [string] with basecalled fast5s is missing. Please provide a valid workspace argument.", call. =FALSE)
  }

  if (missing(nanopolish)) {
    stop("Nanopolish polya output [string] is missing. Please provide a valid nanopolish argument.", .call = FALSE)
  }

  if (missing(sequencing_summary)) {
    stop("Sequencing summary file [string] is missing. Please provide a valid sequencing_summary argument.", .call = FALSE)
  }

  if (missing(num_cores)) {
    stop("Number of declared cores [numeric] is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(nanopolish),msg = "Empty string provided as an input. Please provide a nanopolish as a string")
  assertthat::assert_that(assertive::is_existing_file(nanopolish), msg=paste("File ",nanopolish," does not exist",sep=""))
  assertthat::assert_that(assertive::is_a_non_missing_nor_empty_string(sequencing_summary),msg = "Empty string provided as an input. Please provide a sequencing_summary as a string")
  assertthat::assert_that(assertive::is_existing_file(sequencing_summary), msg=paste("File ",sequencing_summary," does not exist",sep=""))
  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))

  # handle nanopolish
  nanopolish <- vroom::vroom(nanopolish, col_select=c(readname, polya_start, transcript_start, adapter_start, leader_start), show_col_types = FALSE)

  assertthat::assert_that(assertive::has_rows(nanopolish), msg = "Empty data frame provided as an input (nanopolish). Please provide valid input")

  colnames(nanopolish)[1] <- "read_id" #because there was conflict with the same vars
  nanopolish <- dplyr::filter(nanopolish,read_id==readname)
  adapter_start_position <- nanopolish$adapter_start
  leader_start_position <- nanopolish$leader_start
  polya_start_position <- nanopolish$polya_start
  transcript_start_position <-nanopolish$transcript_start
  #define polya end position
  polya_end_position <- transcript_start_position -1

  #handle sequencing summary
  sequencing_summary <- vroom::vroom(sequencing_summary, col_select = c(filename, read_id), show_col_types = FALSE)
  sequencing_summary <- dplyr::filter(sequencing_summary,read_id==readname)

  assertthat::assert_that(assertive::has_rows(sequencing_summary), msg = "Empty data frame provided as an input (sequencing_summary). Please provide valid input")

  # Extract data from fast5 file
  fast5_file <- sequencing_summary$filename
  fast5_readname <- paste0("read_",readname) # in fast5 file structure each readname has suffix "read_"
  fast5_file_path <-file.path(workspace, fast5_file) #this is path to browsed fast5 file

  signal <- rhdf5::h5read(file.path(fast5_file_path),paste0(fast5_readname,"/Raw/Signal"))
  signal <- as.vector(signal) # initially signal is stored as an array, I prefer to vectorize it for further manipulations with dframes

  #winsorize signal (remove cliffs) thanks to this, all signals are plotted without current jets
  signal_q <- stats::quantile(x=signal, probs=c(0.002, 0.998), na.rm=TRUE, type=7)
  minimal_val <- signal_q[1L]
  maximal_val <- signal_q[2L]
  signal[signal<minimal_val] <- minimal_val
  signal[signal>maximal_val] <- maximal_val

  signal <- as.integer(signal)


  #read parameters (attrs) stored in basecall_1d_template
  basecall_1d_template <- rhdf5::h5readAttributes(fast5_file_path,paste0(fast5_readname,"/Analyses/", basecall_group, "/Summary/basecall_1d_template")) # parent dir for attributes (within fast5)
  stride <- basecall_1d_template$block_stride #  this parameter allows to sample data elements along a dimension
  called_events <- basecall_1d_template$called_events # number of events (nanopore translocations) recorded by device for given read

  # close all handled instances (groups, attrs) of fast5 file
  rhdf5::h5closeAll()


  signal_length <- length(signal)

  #creating signal df
  signal_df <- data.frame(position=seq(1,signal_length,1), signal=signal[1:signal_length])
  signal_df$segment[signal_df$position < polya_start_position] <- "adapter"
  signal_df$segment[signal_df$position >= transcript_start_position] <- "transcript"
  signal_df$segment[signal_df$position >= polya_start_position & signal_df$position <= polya_end_position] <- "poly(A)"

  signal_df$segment <- as.factor(signal_df$segment)


  # plotting squiggle

  squiggle <- ggplot2::ggplot(data = signal_df, aes(x = position)) + geom_line(aes(y = signal, color = segment)) + ggplot2::theme_bw() + ggplot2::scale_colour_manual(values = c("#089bcc", "#f56042", "#3a414d"))
  g.line <- ggplot2::geom_vline(xintercept = polya_start_position, color = "#700f25")
  g.line2 <- ggplot2::geom_vline(xintercept = polya_end_position, color = "#0f3473")
  g.labs <- ggplot2::labs(title= paste0("Read ", readname),
                          x="position",
                          y= "signal")

  plot_squiggle <- squiggle + g.line + g.line2 + g.labs



  return(plot_squiggle)

}

