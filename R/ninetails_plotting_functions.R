#' Draws tail range squiggle for given read.
#'
#' Creates segmented plot of raw/rescaled ONT RNA signal (with or without moves).
#' A standalone function; does not rely on any other preprocessing, depends
#' solely on Nanopolish, Guppy and fast5 input.
#'
#' The output plot includes an entire tail region (orange) and +/- 150 positions
#' flanks of adapter (blue) and transcript body (black) regions. Vertical lines
#' mark the 5' (navy blue) and 3' (red) termini of polyA tail according to the
#' Nanopolish polyA function. In order to maintain readability of the graph
#' (and to avoid plotting high cliffs - e.g. jets of the signal caused
#' by a sudden surge of current in the sensor) the signal is winsorised.
#'
#' Moves may be plotted only for reads basecalled by Guppy basecaller.
#' Otherwise the function will throw an error.
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
#' @param moves logical [TRUE/FALSE]. If TRUE, moves would be plotted in
#' the background as vertical bars/gaps (for values 0/1, respectively)
#' and the signal (squiggle) would be plotted in the foreground.
#' Otherwise, only the signal would be plotted. As a default,
#' "FALSE" value is set.
#'
#' @param rescale logical [TRUE/FALSE]. If TRUE, the signal will be rescaled for
#' picoamps (pA) per second (s). If FALSE, raw signal per position will be
#' plotted. As a default, the "FALSE" value is set.
#'
#' @return ggplot2 object with squiggle plot centered on tail range.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' plot_tail_range(nanopolish = '/path/to/file',
#'                 sequencing_summary = '/path/to/file',
#'                 workspace = '/path/to/guppy/workspace',
#'                 basecalled_group = 'Basecall_1D_000',
#'                 num_cores = 3, moves=TRUE, rescale=TRUE)
#'
#'}
#
# TO DO : add rescaling module (ifelse)
plot_tail_range <- function(readname, nanopolish, sequencing_summary, workspace, basecall_group = "Basecall_1D_000", moves=FALSE, rescale=FALSE){

  # variable binding (suppressing R CMD check from throwing an error)
  polya_start <- transcript_start <- adapter_start <- leader_start <- filename <- read_id <- position <- time <- pA <- segment <- NULL


  #Assertions
  if (missing(readname)) {
    stop("Readname [string] is missing. Please provide a valid readname argument.", call. =FALSE)
  }

  if (missing(workspace)) {
    stop("Directory [string] with basecalled fast5s is missing. Please provide a valid workspace argument.", call. =FALSE)
  }
  if (missing(sequencing_summary)) {
    stop("Sequencing summary file [string] is missing. Please provide a valid sequencing_summary argument.", .call = FALSE)
  }

  if (missing(nanopolish)) {
    stop("Nanopolish polya output [string] is missing. Please provide a valid nanopolish argument.", .call = FALSE)
  }


  if (checkmate::test_string(nanopolish)) {
    # if string provided as an argument, read from file
    # handle nanopolish
    assertthat::assert_that(assertive::is_existing_file(nanopolish), msg=paste("File ",nanopolish," does not exist",sep=""))
    nanopolish <- vroom::vroom(nanopolish, col_select=c(readname, polya_start, transcript_start, adapter_start, leader_start), show_col_types = FALSE)
    colnames(nanopolish)[1] <- "read_id" #because there was conflict with the same vars
  }

  #else: make sure that nanopolish is an  object with rows
  assertthat::assert_that(assertive::has_rows(nanopolish), msg = "Empty data frame provided as an input (nanopolish). Please provide valid input")


  nanopolish <- nanopolish[nanopolish$read_id==readname,]
  adapter_start_position <- nanopolish$adapter_start
  leader_start_position <- nanopolish$leader_start
  polya_start_position <- nanopolish$polya_start
  transcript_start_position <-nanopolish$transcript_start
  #define polya end position
  polya_end_position <- transcript_start_position -1


  if (checkmate::test_string(sequencing_summary)) {
    #handle sequencing summary
    assertthat::assert_that(assertive::is_existing_file(sequencing_summary), msg=paste("File ",sequencing_summary," does not exist",sep=""))
    sequencing_summary <- vroom::vroom(sequencing_summary, col_select = c(filename, read_id), show_col_types = FALSE)
  }

  sequencing_summary <- sequencing_summary[sequencing_summary$read_id==readname,]

  assertthat::assert_that(assertive::has_rows(sequencing_summary), msg = "Empty data frame provided as an input (sequencing_summary). Please provide valid input")

  # Extract data from fast5 file
  fast5_file <- sequencing_summary$filename
  fast5_readname <- paste0("read_",readname) # in fast5 file structure each readname has suffix "read_"
  fast5_file_path <-file.path(workspace, fast5_file) #this is path to browsed fast5 file


  signal <- rhdf5::h5read(file.path(fast5_file_path),paste0(fast5_readname,"/Raw/Signal"))
  signal <- as.vector(signal) # initially signal is stored as an array, I prefer to vectorize it for further manipulations with dframes

  # retrieving sequence, traces & move tables
  BaseCalled_template <- rhdf5::h5read(fast5_file_path,paste0(fast5_readname,"/Analyses/", basecall_group, "/BaseCalled_template"))
  move <- BaseCalled_template$Move #how the basecaller "moves" through the called sequence, and allows for a mapping from basecall to raw data
  move <- as.numeric(move) # change data type as moves are originally stored as raw

  #read parameters stored in channel_id group
  channel_id <- rhdf5::h5readAttributes(fast5_file_path,paste0(fast5_readname,"/channel_id")) # parent dir for attributes (within fast5 file)
  digitisation <- channel_id$digitisation # number of quantisation levels in the analog to digital converter
  digitisation <- as.integer(digitisation)
  offset <- channel_id$offset # analog to digital signal error
  offset <- as.integer(offset)
  range <- channel_id$range # difference between the smallest and greatest values
  range <- as.integer(range)
  sampling_rate <- channel_id$sampling_rate # number of data points collected per second
  sampling_rate <- as.integer(sampling_rate)

  #read parameters (attrs) stored in basecall_1d_template
  basecall_1d_template <- rhdf5::h5readAttributes(fast5_file_path,paste0(fast5_readname,"/Analyses/", basecall_group, "/Summary/basecall_1d_template")) # parent dir for attributes (within fast5)
  stride <- basecall_1d_template$block_stride #  this parameter allows to sample data elements along a dimension
  called_events <- basecall_1d_template$called_events # number of events (nanopore translocations) recorded by device for given read

  # close all handled instances (groups, attrs) of fast5 file
  rhdf5::h5closeAll()

  #winsorize signal (remove cliffs) thanks to this, all signals are plotted without current jets
  signal_q <- stats::quantile(x=signal, probs=c(0.002, 0.998), na.rm=TRUE, type=7)
  minimal_val <- signal_q[1L]
  maximal_val <- signal_q[2L]
  signal[signal<minimal_val] <- minimal_val
  signal[signal>maximal_val] <- maximal_val
  signal <- as.integer(signal)


  # number of events expanded for whole signal vec (this is estimation of signal length, however keep in mind that decimal values are ignored)
  number_of_events <- called_events * stride
  signal_length <- length(signal)

  #handling move data
  moves_sample_wise_vector <- c(rep(move, each=stride), rep(NA, signal_length - number_of_events))

  #creating signal df
  signal_df <- data.frame(position=seq(1,signal_length,1), signal=signal[1:signal_length], moves=moves_sample_wise_vector) # this is signal converted to dframe


  # signal segmentation factor
  #adapter sequence
  signal_df$segment[signal_df$position < polya_start_position] <- "adapter"
  #transcript sequence
  signal_df$segment[signal_df$position >= transcript_start_position] <- "transcript"
  #polya sequence
  signal_df$segment[signal_df$position >= polya_start_position & signal_df$position <= polya_end_position] <- "poly(A)"
  signal_df$segment <- as.factor(signal_df$segment)

  # trim tail region:
  trim_position_upstream <- polya_start_position -150
  trim_position_downstream <- transcript_start_position +150

  signal_df <- signal_df[trim_position_upstream:trim_position_downstream,]

  # plotting squiggle

  if (rescale==TRUE) {
    #additional df cols need for rescaling data
    signal_df <- dplyr::mutate(signal_df, time = position/sampling_rate)
    signal_df <- dplyr::mutate(signal_df, pA = ((signal + offset) * (range/digitisation)))
    #plotting signal rescaled to picoamps per second
    squiggle <- ggplot2::ggplot(data = signal_df, ggplot2::aes(x = time)) +
      ggplot2::geom_line(ggplot2::aes(y = pA, color = segment)) + ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = c("#089bcc", "#f56042", "#3a414d"))
    g.line <- ggplot2::geom_vline(xintercept = polya_start_position/sampling_rate, color = "#700f25")
    g.line2 <- ggplot2::geom_vline(xintercept = polya_end_position/sampling_rate, color = "#0f3473")
    g.labs <- ggplot2::labs(title= paste0("Read ", readname),
                            x="time [s]",
                            y= "signal [pA]")
    g.moves <- ggplot2::geom_line(ggplot2::aes(y = moves * 150), alpha = 0.3)

  } else if (rescale==FALSE) {
    #plotting raw signal
    squiggle <- ggplot2::ggplot(data = signal_df, ggplot2::aes(x = position)) +
      ggplot2::geom_line(ggplot2::aes(y = signal, color = segment)) + ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = c("#089bcc", "#f56042", "#3a414d"))
    g.line <- ggplot2::geom_vline(xintercept = polya_start_position, color = "#700f25")
    g.line2 <- ggplot2::geom_vline(xintercept = polya_end_position, color = "#0f3473")
    g.labs <- ggplot2::labs(title= paste0("Read ", readname),
                            x="position",
                            y= "signal [raw]")
    g.moves <- ggplot2::geom_line(ggplot2::aes(y = moves * 1000), alpha = 0.3)
  }

  if (moves==TRUE) {
    plot_squiggle <- squiggle + g.line + g.line2 + g.labs + g.moves
  }
  else if (moves==FALSE) {
    plot_squiggle <- squiggle + g.line + g.line2 + g.labs
  }

  return(plot_squiggle)

}




