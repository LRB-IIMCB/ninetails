################################################################################
# VISUAL INSPECTION OF SIGNALS - FAST5-BASED
################################################################################
#' Draws an entire squiggle for given read.
#'
#' Creates segmented plot of raw/rescaled ONT RNA signal (with or without moves).
#' A standalone function; does not rely on any other preprocessing, depends
#' solely on Nanopolish, Guppy and fast5 input.
#'
#' The output plot includes an entire squiggle corresponding to the given ONT
#' read. Vertical lines mark the 5' (navy blue) and 3' (red) termini of polyA
#' tail according to the Nanopolish polyA function. In order to maintain
#' readability of the graph (and to avoid plotting high cliffs - e.g. jets
#' of the signal caused by a sudden surge of current in the sensor)
#' the signal is winsorised.
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
#' @return ggplot2 object with squiggle plot depicting nanopore read signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#'plot <- ninetails::plot_squiggle_fast5(
#'  readname = "0226b5df-f9e5-4774-bbee-7719676f2ceb",
#'  nanopolish = system.file('extdata',
#'                           'test_data',
#'                           'nanopolish_output.tsv',
#'                           package = 'ninetails'),
#'  sequencing_summary = system.file('extdata',
#'                                   'test_data',
#'                                   'sequencing_summary.txt',
#'                                   package = 'ninetails'),
#'  workspace = system.file('extdata',
#'                          'test_data',
#'                          'basecalled_fast5',
#'                          package = 'ninetails'),
#'  basecall_group = 'Basecall_1D_000',
#'  moves = FALSE,
#'  rescale = TRUE)
#'
#' print(plot)
#'
#'}
#

plot_squiggle_fast5 <- function(readname,
                                nanopolish,
                                sequencing_summary,
                                workspace,
                                basecall_group = "Basecall_1D_000",
                                moves = FALSE,
                                rescale = FALSE) {

  #Assertions
  if (missing(readname)) {
    stop(
      "Readname [string] is missing. Please provide a valid readname argument.",
      call. = FALSE
    )
  }

  if (missing(workspace)) {
    stop(
      "Directory [string] with basecalled fast5s is missing. Please provide a valid workspace argument.",
      call. = FALSE
    )
  }
  if (missing(sequencing_summary)) {
    stop(
      "Sequencing summary file [string] is missing. Please provide a valid sequencing_summary argument.",
      .call = FALSE
    )
  }

  if (missing(nanopolish)) {
    stop(
      "Nanopolish polya output [string] is missing. Please provide a valid nanopolish argument.",
      .call = FALSE
    )
  }

  if (is_string(nanopolish)) {
    # if string provided as an argument, read from file
    # handle nanopolish
    assert_file_exists(nanopolish)
    nanopolish <- vroom::vroom(
      nanopolish,
      col_select = c(
        readname,
        polya_start,
        transcript_start,
        adapter_start,
        leader_start
      ),
      show_col_types = FALSE
    )
    colnames(nanopolish)[1] <- "read_id" #because there was conflict with the same vars
  }

  #else: make sure that nanopolish is an  object with rows
  if (!is.data.frame(nanopolish) || nrow(nanopolish) == 0) {
    stop(
      "Empty data frame provided as an input (nanopolish). Please provide valid input"
    )
  }

  nanopolish <- nanopolish[nanopolish$read_id == readname, ]
  adapter_start_position <- nanopolish$adapter_start
  leader_start_position <- nanopolish$leader_start
  polya_start_position <- nanopolish$polya_start
  transcript_start_position <- nanopolish$transcript_start
  #define polya end position
  polya_end_position <- transcript_start_position - 1

  if (is_string(sequencing_summary)) {
    #handle sequencing summary
    assert_file_exists(sequencing_summary)
    sequencing_summary <- vroom::vroom(
      sequencing_summary,
      col_select = c(filename, read_id),
      show_col_types = FALSE
    )
  }

  sequencing_summary <- sequencing_summary[
    sequencing_summary$read_id == readname,
  ]

  if (!is.data.frame(sequencing_summary) || nrow(sequencing_summary) == 0) {
    stop(
      "Empty data frame provided as an input (sequencing_summary). Please provide valid input"
    )
  }

  # Extract data from fast5 file
  fast5_file <- sequencing_summary$filename
  fast5_readname <- paste0("read_", readname) # in fast5 file structure each readname has suffix "read_"
  fast5_file_path <- file.path(workspace, fast5_file) #this is path to browsed fast5 file

  signal <- rhdf5::h5read(
    file.path(fast5_file_path),
    paste0(fast5_readname, "/Raw/Signal")
  )
  signal <- as.vector(signal) # initially signal is stored as an array, I prefer to vectorize it for further manipulations with dframes

  # retrieving sequence, traces & move tables
  BaseCalled_template <- rhdf5::h5read(
    fast5_file_path,
    paste0(fast5_readname, "/Analyses/", basecall_group, "/BaseCalled_template")
  )
  move <- BaseCalled_template$Move #how the basecaller "moves" through the called sequence, and allows for a mapping from basecall to raw data
  move <- as.numeric(move) # change data type as moves are originally stored as raw

  #read parameters stored in channel_id group
  channel_id <- rhdf5::h5readAttributes(
    fast5_file_path,
    paste0(fast5_readname, "/channel_id")
  ) # parent dir for attributes (within fast5 file)
  digitisation <- channel_id$digitisation # number of quantisation levels in the analog to digital converter
  digitisation <- as.integer(digitisation)
  offset <- channel_id$offset # analog to digital signal error
  offset <- as.integer(offset)
  range <- channel_id$range # difference between the smallest and greatest values
  range <- as.integer(range)
  sampling_rate <- channel_id$sampling_rate # number of data points collected per second
  sampling_rate <- as.integer(sampling_rate)

  #read parameters (attrs) stored in basecall_1d_template
  basecall_1d_template <- rhdf5::h5readAttributes(
    fast5_file_path,
    paste0(
      fast5_readname,
      "/Analyses/",
      basecall_group,
      "/Summary/basecall_1d_template"
    )
  ) # parent dir for attributes (within fast5)
  stride <- basecall_1d_template$block_stride #  this parameter allows to sample data elements along a dimension
  called_events <- basecall_1d_template$called_events # number of events (nanopore translocations) recorded by device for given read

  # close all handled instances (groups, attrs) of fast5 file
  rhdf5::h5closeAll()

  #winsorize signal (remove cliffs) thanks to this, all signals are plotted without current jets
  signal_q <- stats::quantile(
    x = signal,
    probs = c(0.002, 0.998),
    na.rm = TRUE,
    type = 7
  )
  minimal_val <- signal_q[1L]
  maximal_val <- signal_q[2L]
  signal[signal < minimal_val] <- minimal_val
  signal[signal > maximal_val] <- maximal_val
  signal <- as.integer(signal)

  # number of events expanded for whole signal vec (this is estimation of signal length, however keep in mind that decimal values are ignored)
  number_of_events <- called_events * stride
  signal_length <- length(signal)

  #handling move data
  moves_sample_wise_vector <- c(
    rep(move, each = stride),
    rep(NA, signal_length - number_of_events)
  )

  #creating signal df
  signal_df <- data.frame(
    position = seq(1, signal_length, 1),
    signal = signal[1:signal_length],
    moves = moves_sample_wise_vector
  ) # this is signal converted to dframe

  # signal segmentation factor
  #adapter sequence
  signal_df$segment[signal_df$position < polya_start_position] <- "adapter"
  #transcript sequence
  signal_df$segment[
    signal_df$position >= transcript_start_position
  ] <- "transcript"
  #polya sequence
  signal_df$segment[
    signal_df$position >= polya_start_position &
      signal_df$position <= polya_end_position
  ] <- "poly(A)"
  signal_df$segment <- as.factor(signal_df$segment)

  # plotting squiggle

  if (rescale == TRUE) {
    #additional df cols need for rescaling data
    signal_df <- dplyr::mutate(signal_df, time = position / sampling_rate)
    signal_df <- dplyr::mutate(
      signal_df,
      pA = ((signal + offset) * (range / digitisation))
    )
    #plotting signal rescaled to picoamps per second
    squiggle <- ggplot2::ggplot(data = signal_df, ggplot2::aes(x = time)) +
      ggplot2::geom_line(ggplot2::aes(y = pA, color = segment)) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = c("#089bcc", "#f56042", "#3a414d"))
    g.line <- ggplot2::geom_vline(
      xintercept = polya_start_position / sampling_rate,
      color = "#700f25"
    )
    g.line2 <- ggplot2::geom_vline(
      xintercept = polya_end_position / sampling_rate,
      color = "#0f3473"
    )
    g.labs <- ggplot2::labs(
      title = paste0("Read ", readname),
      x = "time [s]",
      y = "signal [pA]"
    )
    g.moves <- ggplot2::geom_line(ggplot2::aes(y = moves * 150), alpha = 0.3)
  } else if (rescale == FALSE) {
    #plotting raw signal
    squiggle <- ggplot2::ggplot(data = signal_df, ggplot2::aes(x = position)) +
      ggplot2::geom_line(ggplot2::aes(y = signal, color = segment)) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = c("#089bcc", "#f56042", "#3a414d"))
    g.line <- ggplot2::geom_vline(
      xintercept = polya_start_position,
      color = "#700f25"
    )
    g.line2 <- ggplot2::geom_vline(
      xintercept = polya_end_position,
      color = "#0f3473"
    )
    g.labs <- ggplot2::labs(
      title = paste0("Read ", readname),
      x = "position",
      y = "signal [raw]"
    )
    g.moves <- ggplot2::geom_line(ggplot2::aes(y = moves * 1000), alpha = 0.3)
  }

  if (moves == TRUE) {
    plot_squiggle <- squiggle + g.line + g.line2 + g.labs + g.moves
  } else if (moves == FALSE) {
    plot_squiggle <- squiggle + g.line + g.line2 + g.labs
  }

  return(plot_squiggle)
}


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
#' plot <- ninetails::plot_tail_range_fast5(
#'   readname = "0226b5df-f9e5-4774-bbee-7719676f2ceb",
#'   nanopolish = system.file('extdata',
#'                            'test_data',
#'                            'nanopolish_output.tsv',
#'                            package = 'ninetails'),
#'   sequencing_summary = system.file('extdata',
#'                                    'test_data',
#'                                    'sequencing_summary.txt',
#'                                    package = 'ninetails'),
#'   workspace = system.file('extdata',
#'                           'test_data',
#'                           'basecalled_fast5',
#'                           package = 'ninetails'),
#'   basecall_group = 'Basecall_1D_000',
#'   moves = TRUE,
#'   rescale = TRUE)
#'
#' print(plot)
#'
#'}
#

plot_tail_range_fast5 <- function(readname,
                                  nanopolish,
                                  sequencing_summary,
                                  workspace,
                                  basecall_group = "Basecall_1D_000",
                                  moves = FALSE,
                                  rescale = FALSE) {

  #Assertions
  if (missing(readname)) {
    stop(
      "Readname [string] is missing. Please provide a valid readname argument.",
      call. = FALSE
    )
  }

  if (missing(workspace)) {
    stop(
      "Directory [string] with basecalled fast5s is missing. Please provide a valid workspace argument.",
      call. = FALSE
    )
  }
  if (missing(sequencing_summary)) {
    stop(
      "Sequencing summary file [string] is missing. Please provide a valid sequencing_summary argument.",
      .call = FALSE
    )
  }

  if (missing(nanopolish)) {
    stop(
      "Nanopolish polya output [string] is missing. Please provide a valid nanopolish argument.",
      .call = FALSE
    )
  }

  if (is_string(nanopolish)) {
    # if string provided as an argument, read from file
    # handle nanopolish
    assert_file_exists(nanopolish)
    nanopolish <- vroom::vroom(
      nanopolish,
      col_select = c(
        readname,
        polya_start,
        transcript_start,
        adapter_start,
        leader_start
      ),
      show_col_types = FALSE
    )
    colnames(nanopolish)[1] <- "read_id" #because there was conflict with the same vars
  }

  #else: make sure that nanopolish is an  object with rows
  if (!is.data.frame(nanopolish) || nrow(nanopolish) == 0) {
    stop(
      "Empty data frame provided as an input (nanopolish). Please provide valid input"
    )
  }

  colnames(nanopolish)[1] <- "read_id" #because there was conflict with the same vars

  nanopolish <- nanopolish[nanopolish$read_id == readname, ]
  adapter_start_position <- nanopolish$adapter_start
  leader_start_position <- nanopolish$leader_start
  polya_start_position <- nanopolish$polya_start
  transcript_start_position <- nanopolish$transcript_start
  #define polya end position
  polya_end_position <- transcript_start_position - 1

  if (is_string(sequencing_summary)) {
    #handle sequencing summary
    assert_file_exists(sequencing_summary)
    sequencing_summary <- vroom::vroom(
      sequencing_summary,
      col_select = c(filename, read_id),
      show_col_types = FALSE
    )
  }

  sequencing_summary <- sequencing_summary[
    sequencing_summary$read_id == readname,
  ]

  if (!is.data.frame(sequencing_summary) || nrow(sequencing_summary) == 0) {
    stop(
      "Empty data frame provided as an input (sequencing_summary). Please provide valid input"
    )
  }

  sequencing_summary <- sequencing_summary[
    sequencing_summary$read_id == readname,
  ]

  # Extract data from fast5 file
  fast5_file <- sequencing_summary$filename
  fast5_readname <- paste0("read_", readname) # in fast5 file structure each readname has suffix "read_"
  fast5_file_path <- file.path(workspace, fast5_file) #this is path to browsed fast5 file

  signal <- rhdf5::h5read(
    file.path(fast5_file_path),
    paste0(fast5_readname, "/Raw/Signal")
  )
  signal <- as.vector(signal) # initially signal is stored as an array, I prefer to vectorize it for further manipulations with dframes

  # retrieving sequence, traces & move tables
  BaseCalled_template <- rhdf5::h5read(
    fast5_file_path,
    paste0(fast5_readname, "/Analyses/", basecall_group, "/BaseCalled_template")
  )
  move <- BaseCalled_template$Move #how the basecaller "moves" through the called sequence, and allows for a mapping from basecall to raw data
  move <- as.numeric(move) # change data type as moves are originally stored as raw

  #read parameters stored in channel_id group
  channel_id <- rhdf5::h5readAttributes(
    fast5_file_path,
    paste0(fast5_readname, "/channel_id")
  ) # parent dir for attributes (within fast5 file)
  digitisation <- channel_id$digitisation # number of quantisation levels in the analog to digital converter
  digitisation <- as.integer(digitisation)
  offset <- channel_id$offset # analog to digital signal error
  offset <- as.integer(offset)
  range <- channel_id$range # difference between the smallest and greatest values
  range <- as.integer(range)
  sampling_rate <- channel_id$sampling_rate # number of data points collected per second
  sampling_rate <- as.integer(sampling_rate)

  #read parameters (attrs) stored in basecall_1d_template
  basecall_1d_template <- rhdf5::h5readAttributes(
    fast5_file_path,
    paste0(
      fast5_readname,
      "/Analyses/",
      basecall_group,
      "/Summary/basecall_1d_template"
    )
  ) # parent dir for attributes (within fast5)
  stride <- basecall_1d_template$block_stride #  this parameter allows to sample data elements along a dimension
  called_events <- basecall_1d_template$called_events # number of events (nanopore translocations) recorded by device for given read

  # close all handled instances (groups, attrs) of fast5 file
  rhdf5::h5closeAll()

  #winsorize signal (remove cliffs) thanks to this, all signals are plotted without current jets
  signal_q <- stats::quantile(
    x = signal,
    probs = c(0.002, 0.998),
    na.rm = TRUE,
    type = 7
  )
  minimal_val <- signal_q[1L]
  maximal_val <- signal_q[2L]
  signal[signal < minimal_val] <- minimal_val
  signal[signal > maximal_val] <- maximal_val
  signal <- as.integer(signal)

  # number of events expanded for whole signal vec (this is estimation of signal length, however keep in mind that decimal values are ignored)
  number_of_events <- called_events * stride
  signal_length <- length(signal)

  #handling move data
  moves_sample_wise_vector <- c(
    rep(move, each = stride),
    rep(NA, signal_length - number_of_events)
  )

  #creating signal df
  signal_df <- data.frame(
    position = seq(1, signal_length, 1),
    signal = signal[1:signal_length],
    moves = moves_sample_wise_vector
  ) # this is signal converted to dframe

  # signal segmentation factor
  #adapter sequence
  signal_df$segment[signal_df$position < polya_start_position] <- "adapter"
  #transcript sequence
  signal_df$segment[
    signal_df$position >= transcript_start_position
  ] <- "transcript"
  #polya sequence
  signal_df$segment[
    signal_df$position >= polya_start_position &
      signal_df$position <= polya_end_position
  ] <- "poly(A)"
  signal_df$segment <- as.factor(signal_df$segment)

  # trim tail region:
  trim_position_upstream <- polya_start_position - 150
  trim_position_downstream <- transcript_start_position + 150

  signal_df <- signal_df[trim_position_upstream:trim_position_downstream, ]

  # plotting squiggle

  if (rescale == TRUE) {
    #additional df cols need for rescaling data
    signal_df <- dplyr::mutate(signal_df, time = position / sampling_rate)
    signal_df <- dplyr::mutate(
      signal_df,
      pA = ((signal + offset) * (range / digitisation))
    )
    #plotting signal rescaled to picoamps per second
    squiggle <- ggplot2::ggplot(data = signal_df, ggplot2::aes(x = time)) +
      ggplot2::geom_line(ggplot2::aes(y = pA, color = segment)) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = c("#089bcc", "#f56042", "#3a414d"))
    g.line <- ggplot2::geom_vline(
      xintercept = polya_start_position / sampling_rate,
      color = "#700f25"
    )
    g.line2 <- ggplot2::geom_vline(
      xintercept = polya_end_position / sampling_rate,
      color = "#0f3473"
    )
    g.labs <- ggplot2::labs(
      title = paste0("Read ", readname),
      x = "time [s]",
      y = "signal [pA]"
    )
    g.moves <- ggplot2::geom_line(ggplot2::aes(y = moves * 150), alpha = 0.3)
  } else if (rescale == FALSE) {
    #plotting raw signal
    squiggle <- ggplot2::ggplot(data = signal_df, ggplot2::aes(x = position)) +
      ggplot2::geom_line(ggplot2::aes(y = signal, color = segment)) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = c("#089bcc", "#f56042", "#3a414d"))
    g.line <- ggplot2::geom_vline(
      xintercept = polya_start_position,
      color = "#700f25"
    )
    g.line2 <- ggplot2::geom_vline(
      xintercept = polya_end_position,
      color = "#0f3473"
    )
    g.labs <- ggplot2::labs(
      title = paste0("Read ", readname),
      x = "position",
      y = "signal [raw]"
    )
    g.moves <- ggplot2::geom_line(ggplot2::aes(y = moves * 1000), alpha = 0.3)
  }

  if (moves == TRUE) {
    plot_squiggle <- squiggle + g.line + g.line2 + g.labs + g.moves
  } else if (moves == FALSE) {
    plot_squiggle <- squiggle + g.line + g.line2 + g.labs
  }

  return(plot_squiggle)
}

################################################################################
# INDEPENDENT OF FAST5 FORMAT
################################################################################

#' Draws a portion of poly(A) tail squiggle (chunk) for given read.
#'
#' This function allows to visualise a single fragment of the poly(A) tail area
#' defined by the fragment name (chunk_name). Intended for use on the output of
#' create_tail_chunk_list_moved function.
#'
#' Currently, draws only the raw signal, without the option to scale
#' to picoamperes [pA].
#'
#' @param chunk_name character string. Name of the given read segment (chunk)
#' within the given tail_chunk_list object.
#'
#' @param tail_chunk_list character string. The list object produced
#' by create_chunk_list function.
#'
#' @return ggplot2 object with squiggle plot depicting given fragment of
#' analyzed poly(A) tail.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#'example <- ninetails::plot_tail_chunk(
#'  chunk_name = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b_1",
#'  tail_chunk_list = tcl)
#'
#' print(example)
#'
#'}

plot_tail_chunk <- function(chunk_name, tail_chunk_list) {
  #assertions
  if (missing(tail_chunk_list)) {
    stop(
      "List of tail chunks is missing. Please provide a valid tail_chunk_list argument.",
      call. = FALSE
    )
  }

  if (missing(chunk_name)) {
    stop(
      "Chunk_name is missing. Please provide a valid chunk_name argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.list(tail_chunk_list),
    paste(
      "Given tail_chunk_list object is not a list. Please provide a valid argument."
    )
  )

  #TODO
  #add assertion checking whether given chunk is present in the tail chunk list or not

  #retrieve actual signal name from chiunk_name
  signal_name <- gsub("\\_.*", "", chunk_name)

  #retrieve the signal values from tail_chunk_list
  signal <- tail_chunk_list[[signal_name]][[chunk_name]][[1]]

  #create signal dataframe for plotting
  signal_df <- data.frame(
    position = seq(1, length(signal), 1),
    signal = signal[1:length(signal)]
  )

  #plot
  plot_squiggle <- ggplot2::ggplot(
    data = signal_df,
    ggplot2::aes(x = position)
  ) +
    ggplot2::geom_line(ggplot2::aes(y = signal), color = "#3271a8") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste0("Read ", chunk_name))

  return(plot_squiggle)
}

################################################################################
# VISUAL REPRESENTATION OF GAFS
################################################################################

#' Creates a visual representation of gramian angular field corresponding to the
#' given poly(A) tail fragment (chunk).
#'
#' The function uses a rainbow palette from grDevices, which makes the matrix
#' more visually pleasing than greyscale (however this is just a visual sugar
#' for human user, as the computer "sees" the data in greyscale
#' (each of 2 dimensions as single-channel matrix).
#'
#' IMPORTANT NOTE!
#' Please keep in mind that the matrices produced by ninetails pipeline are
#' originally two-dimensional arrays (GASF & GADF combined). Each of the
#' dimensions are plotted alltogether (collapsed) as a single depiction.
#' However, they can be splitted, but this requires additional processing steps.
#' For the purpose of ONT signal classification, this combined GASF + GADF
#' approach turned out to be the most suitable, thus the single-dimension
#' extraction is currently not implemented within ninetails plotting functions
#' (which does not mean it would not be).
#'
#' User can control the size of the plot by defining the dimensions within the
#' code chunk in R/RStudio. However, please keep in mind that the ninetails'
#' default built-in model was trained on 100x100 gasfs.
#'
#' @param gaf_name character string. Name of the given read segment (chunk)
#' for which the gaf is meant to be plotted. This is the name of the given gaf
#' within the gaf_list produced by the \code{\link{create_gaf_list}} function.
#'
#' @param gaf_list A list of gaf matrices organized by the read ID_index.
#' @param save_file logical [TRUE/FALSE]. If TRUE, the gaf plot 100x100 (pixels)
#' would be saved in the current workng directory. If FALSE, the plot would be
#' displayed only in the window. This parameter is set to FALSE by default.
#'
#' @return gramian angular field representing given fragment
#' of nanopore read signal.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#'example_gaf <- ninetails::plot_gaf(
#'  gaf_name = "5c2386e6-32e9-4e15-a5c7-2831f4750b2b_1",
#'  gaf_list = gl,
#'  save_file = TRUE)
#'
#'print(example_gaf)
#'
#'}
plot_gaf <- function(gaf_name,
                     gaf_list,
                     save_file = FALSE) {
  #assertions
  if (missing(gaf_name)) {
    stop(
      "Gaf_name is missing. Please provide a valid gaf_name argument.",
      call. = FALSE
    )
  }

  if (missing(gaf_list)) {
    stop(
      "List of GAFs is missing. Please provide a valid gaf_list argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.list(gaf_list),
    paste(
      "Given gaf_list object is not a list. Please provide a valid argument."
    )
  )
  assert_condition(
    gaf_name %in% names(gaf_list),
    "Given gaf_list does not contain provided gaf_name. Please provide a valid gaf_name argument."
  )

  #extract gaf of interest
  gaf <- gaf_list[[gaf_name]]
  #reshape the data so the ggplot will accept their format
  gaf <- reshape2::melt(gaf)

  #plot
  plt <- gaf %>%
    ggplot2::ggplot(ggplot2::aes(Var2, Var1, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colours = c("red", "green", "blue"),
      guide = "none"
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
    ) +
    ggplot2::theme_void()

  if (save_file == TRUE) {
    ggplot2::ggsave(
      paste0(gaf_name, ".png"),
      device = "png",
      width = 100,
      height = 100,
      units = "px",
      bg = "transparent"
    )
  }

  return(plt)
}

#' Creates a visual representation of multiple gramian angular fields
#' based on provided gaf_list (plots all gafs from the given list).
#'
#' The function saves the plots to files with predefined size 100x100 pixels.
#'
#' A feature useful for creating custom data sets for network training.
#' It is recommended to use this feature with caution. Producing multiple graphs
#' from extensive data sets may cause the system to crash.
#'
#' This function plots one-channel (100,100,1) as well as
#' multi-channel gafs (e.g. 100,100,2).
#'
#' IMPORTANT NOTE! In current version, this function plots multi-channel matrices
#' in collapsed manner. If one wants to separate color spaces to GASF/GADF
#' channels or to R,G,B space, this function would not be suitable. The data
#' would require additional processing steps!
#'
#' @param gaf_list A list of gaf matrices organized by the read ID_index.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return multiple png files containing gramian angular fields representing
#' given fragment of nanopore read signal.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' ninetails::plot_multiple_gaf(gaf_list = gl, num_cores = 10)
#'
#'}

plot_multiple_gaf <- function(gaf_list, num_cores) {

  #assertions
  if (missing(gaf_list)) {
    stop(
      "List of GAFs is missing. Please provide a valid gaf_list argument.",
      call. = FALSE
    )
  }

  if (missing(num_cores)) {
    stop(
      "Number of declared cores is missing. Please provide a valid num_cores argument.",
      call. = FALSE
    )
  }

  assert_condition(
    is.list(gaf_list),
    paste(
      "Given gaf_list object is not a list. Please provide a valid argument."
    )
  )
  assert_condition(
    is.numeric(num_cores),
    paste0(
      "Declared core number must be numeric. Please provide a valid argument."
    )
  )

  # creating cluster for parallel computing
  my_cluster <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(my_cluster))

  doSNOW::registerDoSNOW(my_cluster)
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  # this is list of indexes required for parallel computing; the main list is split for chunks
  index_list = split(
    1:length(names(gaf_list)),
    ceiling(1:length(names(gaf_list)) / 100)
  )

  # header for progress bar
  cat(paste('Drawing graphs...', '\n', sep = ''))

  # progress bar
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(index_list),
    style = 3,
    width = 50,
    char = "=",
    file = stderr()
  )

  #create empty list for extracted data
  plot_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)) {
    plot_list <- c(
      plot_list,
      foreach::foreach(nam = names(gaf_list)[index_list[[indx]]]) %dopar%
        ninetails::plot_gaf(nam, gaf_list, save_file = FALSE)
    )
    utils::setTxtProgressBar(pb, indx)
  }

  #label each signal according to corresponding read name to avoid confusion
  gaf_names <- names(gaf_list)
  names(plot_list) <- gaf_names

  #plot all the files
  mapply(
    ggplot2::ggsave,
    file = paste0(names(plot_list), ".png"),
    plot = plot_list,
    device = "png",
    width = 100,
    height = 100,
    units = "px",
    bg = "transparent"
  )
}


################################################################################
# DATA EXPLORATORY PLOTS
################################################################################

#' Plotting read classes data per category assigned to the analyzed reads.
#'
#' This function requires as input the dataframe with read_classes provided by
#' ninetails pipeline. Function works in 3 flavours, by plotting either: \itemize{
#' \item detailed classification (based on column "comments")
#' \item crude classification (based on column "class")
#' \item reads decorated with non-As exclusively
#' }
#'
#' Function  based on the Nanotail equivalent: https://github.com/LRB-IIMCB/nanotail
#' Many thanks to smeagol (Pawel Krawczyk) for advice & support!
#'
#' @param class_data A dataframe or tibble containing read_classes predictions
#' made by ninetails pipeline
#'
#' @param grouping_factor character string. A grouping variable (e.g. "sample_name")
#'
#' @param frequency logical [TRUE/FALSE]. If TRUE, the frequency will be plotted.
#' If FALSE, raw counts will be shown. This parameter is set to TRUE by default.
#'
#' @param type character string ["R"/"N"/"A"]. This variable controls the level
#' of detail of the resulting plot:
#' \itemize{
#' \item "R" - detailed classification (based on column "comments")
#' \item "N" - crude classification (based on column "class")
#' \item "A" - reads decorated with non-As exclusively
#' }
#' By default, the "R" option is set.
#'
#' @return ggplot object with read class prediction
#' @export
#'
#' @examples
#'\dontrun{
#'
#' ninetails::plot_class_counts(class_data=read_classes_dataframe,
#'                              grouping_factor="sample_name",
#'                              frequency=TRUE,
#'                              type="R")
#' }
#'
#'
plot_class_counts <- function(class_data,
                              grouping_factor = NA,
                              frequency = TRUE,
                              type = "R") {

  #assertions
  if (missing(class_data)) {
    stop(
      "Class_data is missing. Please provide a valid class_data argument",
      call. = FALSE
    )
  }

  #else: make sure that it is an object with rows
  if (!is.data.frame(class_data) || nrow(class_data) == 0) {
    stop(
      "Empty data frame provided as an input (class_data). Please provide valid input"
    )
  }

  assert_condition(
    is.logical(frequency),
    "Non-boolean value provided for option frequency"
  )
  assert_condition(
    is.character(type),
    "Non-character argument is not alowed for `type`. Please provide valid type (choose from: 'R', 'N', 'A')"
  )

  if (type == "R") {
    class_counts <- ninetails::count_class(
      class_data = class_data,
      grouping_factor = grouping_factor,
      detailed = TRUE
    )

    basic_colnames = c("comments", "n")
    assert_condition(
      basic_colnames[1] %in% colnames(class_counts),
      "comments column is missing in the input. Invalid output of count_class()."
    )
    assert_condition(
      basic_colnames[2] %in% colnames(class_counts),
      "n column is missing in the input. Invalid output of count_class()."
    )

    class_counts$comments <- factor(
      class_counts$comments,
      levels = c("NIN", "IRL", "QCF", "MPU", "MAU", "YAY")
    )
    if (ncol(class_counts) > 2) {
      grouping_colname = setdiff(colnames(class_counts), basic_colnames)

      class_plot <- ggplot2::ggplot(
        class_counts,
        ggplot2::aes(x = !!rlang::sym(grouping_colname), fill = comments, y = n)
      ) +
        ggplot2::ylab("count") +
        ggplot2::theme_bw() +
        ggplot2::labs(
          x = "\nsample",
          fill = "Class:",
          title = "Assigned read classes"
        )
      if (frequency) {
        class_plot <- class_plot +
          ggplot2::geom_bar(stat = "identity", position = "fill") +
          ggplot2::ylab("frequency")
      } else {
        class_plot <- class_plot +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::ylab("count")
      }
      class_plot <- class_plot +
        ggplot2::labs(
          x = "\nsample",
          fill = "Class:",
          title = "Assigned read classes"
        )
    } else {
      class_plot <- ggplot2::ggplot(
        class_counts,
        ggplot2::aes(x = comments, y = n)
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          x = "\nsample",
          y = "count\n",
          fill = "Class:",
          title = "Assigned read classes"
        )
    }

    class_plot <- class_plot +
      ggplot2::scale_fill_manual(
        values = c(
          "YAY" = "#ff6600",
          "QCF" = "#808080",
          "MAU" = "#00aad4",
          "IRL" = "#4d4d4d",
          "NIN" = "#1a1a1a",
          "MPU" = "#cccccc"
        ),
        labels = c(
          "YAY" = "move transition present, nonA residue detected",
          "QCF" = "nanopolish qc failed",
          "MAU" = "move transition absent, nonA residue undetected",
          "IRL" = "insufficient read length",
          "NIN" = "not included in the analysis (pass only = T)",
          "MPU" = "move transition present, nonA residue undetected"
        )
      )
  } else if (type == "N") {
    class_counts <- ninetails::count_class(
      class_data = class_data,
      grouping_factor = grouping_factor,
      detailed = FALSE
    )

    basic_colnames = c("class", "n")
    assert_condition(
      basic_colnames[1] %in% colnames(class_counts),
      "class column is missing in the input. Invalid output of count_class()."
    )
    assert_condition(
      basic_colnames[2] %in% colnames(class_counts),
      "n column is missing in the input. Invalid output of count_class()."
    )

    class_counts$class <- factor(
      class_counts$class,
      levels = c("unclassified", "blank", "decorated")
    )

    if (ncol(class_counts) > 2) {
      grouping_colname = setdiff(colnames(class_counts), basic_colnames)

      class_plot <- ggplot2::ggplot(
        class_counts,
        ggplot2::aes(x = !!rlang::sym(grouping_colname), fill = class, y = n)
      ) +
        ggplot2::ylab("count") +
        ggplot2::theme_bw() +
        ggplot2::labs(
          x = "\nsample",
          fill = "Class:",
          title = "Assigned read classes"
        )
      if (frequency) {
        class_plot <- class_plot +
          ggplot2::geom_bar(stat = "identity", position = "fill") +
          ggplot2::ylab("frequency")
      } else {
        class_plot <- class_plot +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::ylab("count")
      }
      class_plot <- class_plot +
        ggplot2::labs(
          x = "\nsample",
          fill = "Class:",
          title = "Assigned read classes"
        )
    } else {
      class_plot <- ggplot2::ggplot(
        class_counts,
        ggplot2::aes(x = class, y = n)
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          x = "\nsample",
          y = "count\n",
          fill = "Class:",
          title = "Assigned read classes"
        )
    }

    class_plot <- class_plot +
      ggplot2::scale_fill_manual(
        values = c(
          "decorated" = "#ff6600",
          "blank" = "#00aad4",
          "unclassified" = "#808080"
        )
      )
  } else if (type == "A") {
    class_counts <- ninetails::count_class(
      class_data = class_data,
      grouping_factor = grouping_factor,
      detailed = FALSE
    )

    basic_colnames = c("class", "n")
    assert_condition(
      basic_colnames[1] %in% colnames(class_counts),
      "class column is missing in the input. Invalid output of count_class()."
    )
    assert_condition(
      basic_colnames[2] %in% colnames(class_counts),
      "n column is missing in the input. Invalid output of count_class()."
    )

    class_counts$class <- factor(
      class_counts$class,
      levels = c("unclassified", "blank", "decorated")
    )

    if (ncol(class_counts) > 2) {
      grouping_colname = setdiff(colnames(class_counts), basic_colnames)

      if (frequency) {
        class_counts <- class_counts %>%
          dplyr::left_join(
            class_counts %>%
              dplyr::group_by(!!rlang::sym(grouping_colname)) %>%
              dplyr::summarize(total = sum(n))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::filter(class == "decorated") %>%
          dplyr::mutate(prop = n / total)

        class_plot <- ggplot2::ggplot(
          class_counts,
          ggplot2::aes(
            x = !!rlang::sym(grouping_colname),
            fill = class,
            y = prop
          )
        ) +
          ggplot2::ylab("count") +
          ggplot2::theme_bw() +
          ggplot2::labs(
            x = "\nsample",
            fill = "Class:",
            title = "Reads with non-As in poly(A)"
          ) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::ylab("frequency")
      } else {
        class_plot <- ggplot2::ggplot(
          class_counts %>% dplyr::filter(class == "decorated"),
          ggplot2::aes(x = !!rlang::sym(grouping_colname), fill = class, y = n)
        ) +
          ggplot2::ylab("count") +
          ggplot2::theme_bw() +
          ggplot2::labs(
            x = "\nsample",
            fill = "Class:",
            title = "Reads with non-As in poly(A)"
          ) +
          ggplot2::geom_bar(stat = "identity", position = "stack") +
          ggplot2::ylab("count")
      }
    } else {
      class_plot <- ggplot2::ggplot(
        class_counts %>% dplyr::filter(class == "decorated"),
        ggplot2::aes(x = class, y = n)
      ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          x = "\nsample",
          y = "count\n",
          fill = "Class:",
          title = "Reads with non-As in poly(A)"
        )
    }

    class_plot <- class_plot +
      ggplot2::scale_fill_manual(values = c("decorated" = "#ff6600"))
  } else {
    stop(
      "Wrong type of plot defined. Please provide a valid type (choose from: 'R', 'N', 'A')",
      call. = FALSE
    )
  }

  return(class_plot)
}


#' Plot counts of nonadenosine residues found in ninetails output data
#'
#' This function may report the non-A occurrences in two flavours: by read
#' or by residue. To illustrate this, let's have a look at the following example:
#'
#' read_1:
#' contains 3xU residues
#'
#' read_2:
#' contains 2xC residue, 1xG residue, 1xU residue
#'
#' function launched in by read mode (by_read==TRUE):
#' read_1 - taken into account once, because it non-As of the same (single) type
#' read_2 - taken into account thrice, because it has non-As of three types
#'
#' function launched in by residue mode (by_read==FALSE):
#' read_1 - taken into account thrice, because it has 3 non-A residues reported
#' read_2 - taken into account fourfold, because it has 4 nonAs reported
#'
#' The user can switch between those modes depending on the desired informations
#' to be reported.
#'
#' @param residue_data A dataframe or tibble containig non-A residue predictions
#' made by ninetails pipeline
#'
#' @param grouping_factor grouping variable (e.g. "sample_name")
#'
#' @param by_read logical [TRUE/FALSE]. If TRUE, the count/frequency data
#' per reads with given residues will be plotted (i.e. how many reads contain
#' given residue. If FALSE, the number of residues will be plotted).
#' Set to FALSE by default.
#'
#' @param frequency logical [TRUE/FALSE]. If TRUE, the frequency will be plotted.
#' If FALSE, raw counts will be shown. This parameter is set to TRUE by default.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'\dontrun{
#' ninetails::plot_residue_counts(residue_data=nonadenosine_residues_dataframe,
#'                                grouping_factor="sample_name",
#'                                frequency=TRUE)
#' }
plot_residue_counts <- function(residue_data,
                                grouping_factor = NA,
                                by_read = FALSE,
                                frequency = TRUE) {

  #assertions
  if (missing(residue_data)) {
    stop(
      "Residue_data dataframe is missing. Please provide a valid residue_data argument",
      call. = FALSE
    )
  }

  #else: make sure that it is an object with rows
  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  if (by_read == TRUE) {
    residue_data2 <- residue_data[
      !duplicated(residue_data[c("readname", "prediction")]),
    ]
    residue_counts <- ninetails::count_residues(
      residue_data = residue_data2,
      grouping_factor = grouping_factor
    )
    title_suffix <- "reported by read"
  } else {
    residue_counts <- ninetails::count_residues(
      residue_data = residue_data,
      grouping_factor = grouping_factor
    )
    title_suffix <- "reported by residue"
  }

  basic_colnames = c("prediction", "n")
  assert_condition(
    basic_colnames[1] %in% colnames(residue_counts),
    "prediction column is missing in the input. Invalid output of count_residues()."
  )
  assert_condition(
    basic_colnames[2] %in% colnames(residue_counts),
    "n column is missing in the input. Invalid output of count_residues()."
  )

  assert_condition(
    is.logical(by_read),
    "Non-boolean value provided for option by_read"
  )
  assert_condition(
    is.logical(frequency),
    "Non-boolean value provided for option frequency"
  )

  # if there were multiple samples compared
  if (ncol(residue_counts) > 2) {
    grouping_colname = setdiff(colnames(residue_counts), basic_colnames)

    residue_plot <- ggplot2::ggplot(
      residue_counts,
      ggplot2::aes(x = !!rlang::sym(grouping_colname), fill = prediction, y = n)
    ) +
      ggplot2::ylab("count") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "\nsample",
        fill = "Residue:",
        title = paste0("Non-A residue counts ", title_suffix)
      )

    if (frequency) {
      #residue_plot <- residue_plot + ggplot2::geom_bar(stat="identity",position="fill") + ggplot2::ylab("frequency")
      residue_counts <- residue_counts %>%
        tidyr::pivot_wider(names_from = prediction, values_from = n) %>%
        dplyr::mutate(
          total = C + G + U,
          pC = C / total,
          pG = G / total,
          pU = U / total
        ) %>%
        dplyr::select(group, pC, pG, pU) %>%
        dplyr::rename(C = pC, G = pG, U = pU) %>%
        tidyr::gather(key = prediction, value = n, C:G:U)

      #unify factor levels
      residue_counts$prediction <- factor(
        residue_counts$prediction,
        levels = c("C", "G", "U")
      )

      residue_plot <- ggplot2::ggplot(
        residue_counts,
        ggplot2::aes(
          x = !!rlang::sym(grouping_colname),
          fill = prediction,
          y = n
        )
      ) +
        ggplot2::ylab("frequency") +
        ggplot2::theme_bw() +
        ggplot2::labs(
          x = "\nsample",
          fill = "Residue:",
          title = paste0("Non-A residue frequency ", title_suffix)
        ) +
        ggplot2::geom_bar(stat = "identity", position = "dodge")
    } else {
      residue_plot <- residue_plot +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::ylab("count")
    }

    residue_plot <- residue_plot +
      ggplot2::labs(
        x = "\nsample",
        fill = "Residue:",
        title = paste0("Non-A residue counts ", title_suffix)
      )
  } else {
    residue_plot <- ggplot2::ggplot(
      residue_counts,
      ggplot2::aes(x = prediction, y = n)
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(
        x = "\nsample",
        y = "count\n",
        fill = "Residue:",
        title = paste0("Non-A residue counts ", title_suffix)
      )
  }

  #color scheme
  residue_plot <- residue_plot +
    ggplot2::scale_fill_manual(values = c("#3a424f", "#50a675", "#b0bdd4"))

  return(residue_plot)
}

#' Plots qc data (qc_tag) inherited from nanopolish polya function.
#'
#' Provides insight into the proportion of each quality tag category per sample
#' or experiment condition defined by the user.
#'
#' This is the ninetails' implementation of the
#' \code{\link[nanotail:plot_nanopolish_qc]{name}} function originally written
#' by P. Krawczyk (smeagol) and incorporated within the NanoTail package.
#'
#' For original source code, see:
#' https://github.com/LRB-IIMCB/nanotail/blob/master/R/polya_plots.R
#'
#' The variable names were adjusted according to the naming convention within
#' ninetails to avoid confusion.
#'
#' Many thanks to the author of original source code for help and advice.
#'
#' @param processing_info the output of nanopolish_qc function
#'
#' @param frequency logical [TRUE/FALSE].
#'
#' @return a ggplot object
#' @export
#'
#'
plot_nanopolish_qc <- function(processing_info, frequency = TRUE) {

  #assertions
  if (missing(processing_info)) {
    stop(
      "Nanopolish processing info is missing. Please provide a valid processing_info argument",
      call. = FALSE
    )
  }

  if (!is.data.frame(processing_info) || nrow(processing_info) == 0) {
    stop(
      "Empty data frame provided as an input (processing_info). Please provide valid input"
    )
  }

  basic_colnames = c("qc_tag", "n")
  assert_condition(
    basic_colnames[1] %in% colnames(processing_info),
    "qc_tag column is missing in the input. Is that valid output of nanopolish_qc()?"
  )
  assert_condition(
    basic_colnames[2] %in% colnames(processing_info),
    "n column is missing in the input. Is that valid output of nanopolish_qc()?"
  )
  assert_condition(
    is.logical(frequency),
    "Non-boolean value provided for option frequency"
  )

  # if there were multiple samples compared
  if (ncol(processing_info) > 2) {
    grouping_colname = setdiff(colnames(processing_info), basic_colnames)
    nanopolish_qc_plot <- ggplot2::ggplot(
      processing_info,
      ggplot2::aes(x = !!rlang::sym(grouping_colname), fill = qc_tag, y = n)
    )
    if (frequency) {
      nanopolish_qc_plot <- nanopolish_qc_plot +
        ggplot2::geom_bar(stat = "identity", position = "fill") +
        ggplot2::labs(
          x = "\nsample",
          y = "frequency\n",
          fill = "Nanopolish QC tag:"
        ) +
        ggplot2::ggtitle("Read frequency plot")
    } else {
      nanopolish_qc_plot <- nanopolish_qc_plot +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::labs(
          x = "\nsample",
          y = "count\n",
          fill = "Nanopolish QC tag:"
        ) +
        ggplot2::ggtitle("Read count plot")
    }
    nanopolish_qc_plot <- nanopolish_qc_plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.7)
      )
  } else {
    nanopolish_qc_plot <- ggplot2::ggplot(
      processing_info,
      ggplot2::aes(x = qc_tag, y = n)
    ) +
      ggplot2::geom_bar(stat = "identity")
  }

  #define color values & theme
  nanopolish_qc_plot <- nanopolish_qc_plot +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(
      values = c("#D7191C", "#FDAE61", "#FFFFBF", "#ABD9E9", "#2C7BB6")
    )

  return(nanopolish_qc_plot)
}


#' Plots poly(A) tail length (or estimated non-A position) distribution
#' in analyzed sample(s).
#'
#' This function plots distributions of either poly(A) tail lengths or estimated
#' non-A residues' positions across the dataset using the user-predefned grouping
#' variable (e.g. samples, conditions etc.). The grouping variable must be a
#' column within the input dataset passed to the function.
#'
#' User can specify this in samples_table if the ninetails pipeline output is
#' intended to be loaded into R session by \code{\link{read_residue_multiple}}
#' and \code{\link{read_class_multiple}} functions or added manually by
#' other means.
#'
#' The function takes as an input merged ninetails' output dataset, which can be
#' produced with the \code{\link{merge_nonA_tables}} function.
#'
#' The function allows to mark measures of central tendency, either mean, median
#' or mode (as vertical dashed line). One of these values may be shown on the plot
#' with additional caption in lower right corner. This is an optional feature.
#'
#' Mean and median are computed using base R functions.
#'
#' The mode value (if specified by value_to_show argument) is computed using
#' the \code{\link{get_mode}} function. This is a custom helper function written
#' based on this thread: https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
#'
#'
#' This returns the density mode for
#' normalized data to avoid unexpected results (messed-up plots) which may occur
#' if the distribution is bi- or multimodal.
#'
#' This is the ninetails' implementation of plot_polya_distribution()
#' from Nanotail backage by P. Krawczyk (smaegol) written under
#' author's permission.
#'
#' The original code is available here:
#' https://github.com/smaegol/nanotail/blob/master/R/polya_plots.R
#'
#' The function was simplified & adjusted to ninetails' naming convention
#' to avoid confusion & overwriting, if both packages are loaded in the same
#' R session.
#'
#' @param input_data the ninetails pipeline output data preprocessed with the
#' \code{\link{merge_nonA_tables}} function or residue_data or other suitable data
#'
#' @param variable_to_plot [character] string, the variable to be plotted defined
#' by the user. By default (if not provided by the user) the polya_length would
#' be plotted. Note that this variable has to be the column of input_data.
#'
#' @param grouping_factor [character] string, the grouping variable defined
#' by the user. Note that this variable has to be the column of input_data.
#'
#' @param max_length [numeric] maximum length of plotted tail data
#'
#' @param value_to_show [character] string; one of the measures of central
#' tendency: either the "mode", the "median" or the "mean" value, which user
#' wants to be displayed on the plot. By default, none is specified.
#'
#' @param ndensity logical [TRUE/FALSE]. If TRUE, the normalized density would
#' be shown. If false - the data would not be normalized. It is set to "TRUE"
#' by default.
#'
#' @param title logical [TRUE/FALSE]. If TRUE, the title + subtitles would be
#' displayed. If not - the title & subtitle would not be visible.
#'
#' @return a ggplot object
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' plt <- ninetails::plot_tail_distribution(input_data = merged_nonA_tables,
#'                                          variable_to_plot = "polya_length",
#'                                          grouping_factor = "group",
#'                                          max_length = 200,
#'                                          value_to_show = "median",
#'                                          ndensity=T,
#'                                          title=F)
#'
#'
#'}
plot_tail_distribution <- function(input_data,
                                   variable_to_plot = "polya_length",
                                   grouping_factor = NA,
                                   max_length = NA,
                                   value_to_show = NA,
                                   ndensity = T,
                                   title = F) {

  # Assertions
  if (missing(input_data)) {
    stop(
      "Ninetails data are missing. Please provide a valid input_data argument",
      call. = FALSE
    )
  }

  assert_condition(
    is.character(variable_to_plot),
    paste0(
      "Variable_to_plot must be a string. Please provide a valid argument."
    )
  )

  if (variable_to_plot == "polya_length") {
    x_caption <- ggplot2::xlab("poly(A) length")
    plot_title <- "Poly(A) length distribution"
  } else if (variable_to_plot == "est_nonA_pos") {
    x_caption <- ggplot2::xlab("estimated non-A position")
    plot_title <- "Non-A residues distribution"
  } else {
    x_caption <- ggplot2::xlab(variable_to_plot)
    plot_title <- paste0(variable_to_plot, "distribution")
  }

  if (!is.na(grouping_factor)) {
    assert_condition(
      grouping_factor %in% colnames(input_data),
      paste0(grouping_factor, " is not a column of input dataset")
    )

    plot_tails <- ggplot2::ggplot(
      input_data,
      ggplot2::aes(
        x = !!rlang::sym(variable_to_plot),
        color = !!rlang::sym(grouping_factor)
      )
    ) +
      x_caption
  } else {
    plot_tails <- ggplot2::ggplot(
      input_data,
      ggplot2::aes(x = !!rlang::sym(variable_to_plot))
    ) +
      x_caption
  }

  if (!is.na(max_length)) {
    assert_condition(
      is.numeric(max_length),
      "Please provide numeric value for max_length"
    )
    plot_tails <- plot_tails +
      ggplot2::scale_x_continuous(limits = c(0, max_length))
  }

  if (ndensity) {
    plot_tails <- plot_tails +
      ggplot2::geom_line(
        stat = "density",
        size = 1,
        ggplot2::aes(y = ..ndensity..)
      ) +
      ggplot2::ylab("normalized density")
  } else {
    plot_tails <- plot_tails +
      ggplot2::geom_line(
        stat = "density",
        size = 1,
        ggplot2::aes(y = ..density..)
      ) +
      ggplot2::ylab("density")
  }

  if (!is.na(value_to_show)) {
    if (value_to_show == "median") {
      center_value = "median_value"
      center_values <- input_data %>%
        dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
        dplyr::summarize(
          median_value = stats::median(polya_length, na.rm = TRUE)
        )
    } else if (value_to_show == "mean") {
      center_value = "mean_value"
      center_values <- input_data %>%
        dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
        dplyr::summarize(mean_value = mean(polya_length, na.rm = TRUE))
    } else if (value_to_show == "mode") {
      center_value = "mode_value"
      center_values <- input_data %>%
        dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
        dplyr::summarize(
          mode_value = get_mode(polya_length, na.rm = TRUE, method = "density")
        )
    } else {
      stop(
        "Wrong value_to_show specified (should be 'median', 'mode', 'mean' or none (NA)"
      )
    }

    plot_tails <- plot_tails +
      ggplot2::geom_vline(
        data = center_values,
        ggplot2::aes(
          xintercept = !!rlang::sym(center_value),
          color = !!rlang::sym(grouping_factor)
        ),
        linetype = "longdash",
        show.legend = FALSE
      )
  } else {
    plot_tails <- plot_tails
  }

  plot_tails <- plot_tails +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw()

  if (title) {
    plot_tails <- plot_tails +
      ggplot2::labs(
        title = plot_title,
        caption = paste0("Dashed lines - ", value_to_show, " values\n")
      )
  }

  return(plot_tails)
}

#' Plots panel characteristics of ninetails output.
#'
#' Creates a multipanel plot with comprehensive characteristics of input data
#' produced by ninetails software. Those panel charts provide the most comprehensive
#' characterization of a given pool of reads (representing particular transcript
#' or set of transcripts, respectively).
#'
#' Includes 5 panels A-E containing various subplots with the following content:\itemize{
#' \item A - read categories - result of read classification according to presence/absence
#' of non-As in their poly(A) tails
#' \item B - non-A residues - frequency of reads containing given residue among all
#' of the reads harboring non-A nucleotides
#' \item C - distribution of lengths of poly(A) tails - overall distribution of lengths
#' in comparison to the lengths of tails decorated with non-As
#' \item D - normalized distribution of non-A residues in poly(A) tails -
#' non-A nucleotide positions normalized to the length of reads in which given residue occurs
#' \item raw distribution of non-A residues in poly(A) tails - crude depiction
#' of positions along the tail range
#' }
#'
#' @param input_residue_data A dataframe or tibble containig non-A residue predictions
#' made by ninetails pipeline
#'
#' @param input_class_data A dataframe or tibble containing read_classes predictions
#' made by ninetails pipeline. Mutually exclusive with parameter input_merged_nonA_tables_data
#'
#' @param input_merged_nonA_tables_data A dataframe or tibble containing merged_nonA_tables data
#' produced by the \code{\link{merge_nonA_tables}} function. Mutually exclusive with parameter
#' input_class_data
#'
#' @param type [character] either "default" or "moderna" - controls the predefined
#' settings of the plot. In case of moderna adds the default UCUAG pentamer position
#' to the distribution subplots. By default, the "default" option is set.
#'
#' @param max_length [numeric] maximum length of plotted tail data
#'
#' @param direction_5_prime [logical] by default set to TRUE. controls the direction in which
#' the non-A positions are reported. Either from 5' end, according to the most common convention,
#' or from 3' end (as ONT direct RNA seq proceeds).
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#'\dontrun{
#'
#' ninetails::plot_panel_characteristics(input_residue_data=residue_data,
#'                                       input_class_data=class_data,
#'                                       input_merged_nonA_tables_data=NULL,
#'                                       type="default",
#'                                       max_length=100,
#'                                       direction_5_prime=T)
#'
#'}
#'
plot_panel_characteristics <- function(input_residue_data,
                                       input_class_data = NULL,
                                       input_merged_nonA_tables_data = NULL,
                                       type = "default",
                                       max_length = 300,
                                       direction_5_prime = TRUE) {

  # ASSERTIONS
  ##############################################################################

  if (!is.na(max_length)) {
    assert_condition(
      is.numeric(max_length),
      "Please provide numeric value for max_length"
    )
  }

  if (!is.data.frame(input_residue_data) || nrow(input_residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (input_residue_data). Please provide valid input"
    )
  }

  assert_condition(
    is.character(type),
    "Type must be a string (either 'default' or 'moderna'). Please provide a valid argument."
  )

  if (is.null(input_class_data) && is.null(input_merged_nonA_tables_data)) {
    stop(
      "At least one dataframe should be provided - either input_class_data or input_merged_nonA_tables_data"
    )
  }
  if (!is.null(input_class_data) && !is.null(input_merged_nonA_tables_data)) {
    stop(
      "Only one dataframe should be provided - either input_class_data or input_merged_nonA_tables_data"
    )
  }

  # PROCESSING INPUTS
  ##############################################################################

  if (!is.null(input_class_data)) {
    input_merged_nonA_tables_data <- ninetails::merge_nonA_tables(
      residue_data = input_residue_data,
      class_data = input_class_data,
      pass_only = F
    )
  }
  if (!is.null(input_merged_nonA_tables_data)) {
    input_merged_nonA_tables_data <- input_merged_nonA_tables_data
  }

  ## TAIL DISTRIBUTION DATA

  colnames_to_save <- c(
    "sample",
    "group",
    "readname",
    "prediction",
    "est_nonA_pos",
    "polya_length",
    "class",
    "comments",
    "transcript",
    "ensembl_transcript_id_full",
    "ensembl_transcript_id_short",
    "prediction_C",
    "prediction_G",
    "prediction_U",
    "nonA_residues"
  )

  tail_distribution_data <- input_merged_nonA_tables_data %>%
    dplyr::select(
      -dplyr::one_of(setdiff(
        names(input_merged_nonA_tables_data),
        colnames_to_save
      ))
    )

  tail_distribution_data <- dplyr::bind_rows(
    tail_distribution_data %>%
      dplyr::filter(is.na(nonA_residues)) %>%
      dplyr::mutate(type = "blank"),
    tail_distribution_data %>%
      dplyr::filter(!is.na(nonA_residues)) %>%
      dplyr::mutate(type = "nonA"),
    tail_distribution_data %>% dplyr::mutate(type = "total")
  )

  # factor level order! hotfix
  tail_distribution_data$type <- factor(
    tail_distribution_data$type,
    levels = c("nonA", "total", "blank"),
    ordered = TRUE
  )

  ## RESIDUE DATA
  #clean the table from unneeded cols
  colnames_to_save <- c(
    "readname",
    "prediction",
    "est_nonA_pos",
    "polya_length"
  )
  input_residue_data <- input_residue_data %>%
    dplyr::select(
      -dplyr::one_of(setdiff(names(input_residue_data), colnames_to_save))
    ) %>%
    dplyr::group_by(prediction) %>%
    dplyr::add_count() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = paste0("(n=", n, ")"))

  if (direction_5_prime == TRUE) {
    # getting rid of decimals
    input_residue_data$round_pos <- round(input_residue_data$est_nonA_pos)
    input_residue_data$round_length <- round(input_residue_data$polya_length)

    # bin columns
    upper_limit <- 5 * ceiling(max(input_residue_data$round_pos) / 5)
    input_residue_data$binned_positions <- cut(
      input_residue_data$round_pos,
      breaks = seq(0, upper_limit, by = 5),
      right = TRUE
    )
    upper_limit <- 5 * ceiling(max(input_residue_data$round_length) / 5)
    input_residue_data$binned_lengths <- cut(
      input_residue_data$round_length,
      breaks = seq(0, upper_limit, by = 5),
      right = TRUE
    )

    #convert intervals into numeric cols with upper value
    input_residue_data$binned_lengths <- gsub(
      ".*\\,",
      "",
      input_residue_data$binned_lengths
    )
    input_residue_data$binned_lengths <- gsub(
      '.{1}$',
      '',
      input_residue_data$binned_lengths
    )
    input_residue_data$binned_lengths <- as.numeric(
      input_residue_data$binned_lengths
    )
    input_residue_data$binned_positions <- gsub(
      ".*\\,",
      "",
      input_residue_data$binned_positions
    )
    input_residue_data$binned_positions <- gsub(
      '.{1}$',
      '',
      input_residue_data$binned_positions
    )
    input_residue_data$binned_positions <- as.numeric(
      input_residue_data$binned_positions
    )

    #adding columns for cumulative counts per positions
    input_residue_data <- input_residue_data %>%
      dplyr::mutate(
        counts_of_reads_equal_or_longer_than_est_position = sapply(
          binned_positions,
          function(x) sum(x <= binned_lengths)
        )
      ) %>%
      dplyr::group_by(binned_positions) %>%
      dplyr::mutate(
        counts_of_reads_with_nonA_in_given_position = dplyr::n()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        ygreki = counts_of_reads_with_nonA_in_given_position /
          counts_of_reads_equal_or_longer_than_est_position
      )

    group_count <- length(unique(input_residue_data$binned_lengths))

    subtitle_info <- "Non-A positions reported from 5' end"
  } else if (direction_5_prime == FALSE) {
    # replace est_nonA_pos
    input_residue_data <- input_residue_data %>%
      dplyr::mutate(est_nonA_pos_2 = polya_length - est_nonA_pos) %>%
      dplyr::select(-est_nonA_pos) %>%
      dplyr::rename(est_nonA_pos = est_nonA_pos_2)

    # getting rid of decimals
    input_residue_data$round_pos <- round(input_residue_data$est_nonA_pos)
    input_residue_data$round_length <- round(input_residue_data$polya_length)

    # bin columns
    upper_limit <- 5 * ceiling(max(input_residue_data$round_pos) / 5)
    input_residue_data$binned_positions <- cut(
      input_residue_data$round_pos,
      breaks = seq(0, upper_limit, by = 5),
      right = TRUE
    )
    upper_limit <- 5 * ceiling(max(input_residue_data$round_length) / 5)
    input_residue_data$binned_lengths <- cut(
      input_residue_data$round_length,
      breaks = seq(0, upper_limit, by = 5),
      right = TRUE
    )

    #convert intervals into numeric cols with upper value
    input_residue_data$binned_lengths <- gsub(
      ".*\\,",
      "",
      input_residue_data$binned_lengths
    )
    input_residue_data$binned_lengths <- gsub(
      '.{1}$',
      '',
      input_residue_data$binned_lengths
    )
    input_residue_data$binned_lengths <- as.numeric(
      input_residue_data$binned_lengths
    )
    input_residue_data$binned_positions <- gsub(
      ".*\\,",
      "",
      input_residue_data$binned_positions
    )
    input_residue_data$binned_positions <- gsub(
      '.{1}$',
      '',
      input_residue_data$binned_positions
    )
    input_residue_data$binned_positions <- as.numeric(
      input_residue_data$binned_positions
    )

    #adding columns for cumulative counts per positions
    input_residue_data <- input_residue_data %>%
      dplyr::mutate(
        counts_of_reads_equal_or_longer_than_est_position = sapply(
          binned_positions,
          function(x) sum(x <= binned_lengths)
        )
      ) %>%
      dplyr::group_by(binned_positions) %>%
      dplyr::mutate(
        counts_of_reads_with_nonA_in_given_position = dplyr::n()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        ygreki = counts_of_reads_with_nonA_in_given_position /
          counts_of_reads_equal_or_longer_than_est_position
      )

    group_count <- length(unique(input_residue_data$binned_lengths))

    subtitle_info <- "Non-A positions reported from 3' end"
  } else {
    stop(
      "Unknown argument defined. Please provide correct parameter",
      .call = FALSE
    )
  }

  #drop unnecessary columns
  input_residue_data$round_pos <- NULL
  input_residue_data$round_length <- NULL

  ## READ COUNTS FOR METRIC PLOT
  summarized_nonA <- ninetails::summarize_nonA(
    input_merged_nonA_tables_data,
    summary_factors = c("group"),
    transcript_id_column = NULL
  )

  ### split for 2 separate dframes (count, hits)
  summarized_nonA <- lapply(
    split.default(
      summarized_nonA[-(1:3)],
      sub("_.*", "", names(summarized_nonA)[-(1:3)])
    ),
    function(x) cbind(summarized_nonA[1:3], x)
  )
  ###extract only counts
  summarized_nonA <- summarized_nonA[[1]]
  summarized_nonA <- summarized_nonA %>%
    dplyr::select(-polya_median, -polya_mean) %>%
    dplyr::mutate_if(is.numeric, ~ round(., 3)) %>%
    tidyr::pivot_longer(
      cols = c(
        counts_total,
        counts_blank,
        counts_nonA,
        counts_C,
        counts_G,
        counts_U
      ),
      names_to = "source",
      values_to = "counts"
    )
  summarized_nonA$cat <- c(rep("main", 3), rep("mods", 3))

  ### divide into 2 separate subplots for general categories & c,g,u preds
  main_metrics <- summarized_nonA %>%
    dplyr::filter(cat == "main") %>%
    dplyr::mutate(label = paste0("(n=", counts, ")"))
  main_metrics$source <- factor(
    main_metrics$source,
    levels = c("counts_nonA", "counts_blank", "counts_total")
  )

  cgu_metrics <- summarized_nonA %>%
    dplyr::filter(cat == "mods") %>%
    dplyr::mutate(label = paste0("(n=", counts, ")"))
  cgu_metrics$source <- factor(
    cgu_metrics$source,
    levels = c("counts_C", "counts_G", "counts_U")
  )

  # PLOTTING
  ##############################################################################

  # general read classification plots
  ## general read categories:
  general_read_categories <- ggplot2::ggplot(
    data = main_metrics,
    ggplot2::aes(x = source, y = counts, fill = source)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(
      values = c("#ff6600", "#00aad4", "#174e73"),
      labels = c("with non-As", "blank", "total"),
      guide = ggplot2::guide_legend(reverse = TRUE)
    ) +
    ggplot2::scale_x_discrete(
      breaks = main_metrics$source,
      labels = main_metrics$label
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = "Read categories", tag = "A")

  ## residue freq per read:
  residue_counts <- ggplot2::ggplot(
    data = cgu_metrics,
    ggplot2::aes(x = forcats::fct_rev(source), y = counts, fill = source)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(
      values = c("#3a424f", "#50a675", "#b0bdd4"),
      labels = c("reads with C", "reads with G", "reads with U")
    ) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(
      breaks = cgu_metrics$source,
      labels = cgu_metrics$label
    ) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = "Non-A residues", tag = "B") +
    ggplot2::coord_flip()

  # distribution plot:
  distrib_plot <- ninetails::plot_tail_distribution(
    input_data = tail_distribution_data,
    variable_to_plot = "polya_length",
    max_length = max_length,
    ndensity = F,
    title = F,
    grouping_factor = "type"
  ) +
    ggplot2::scale_color_manual(
      values = c("#ff6600", "#174e73", "#00aad4"),
      labels = c("with non-As", "total", "blank")
    ) +
    ggplot2::labs(
      title = "Distribution of lengths of poly(A) tails",
      color = "read type",
      tag = "C"
    ) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )

  # NORMALIZED PLOT
  binned_length_pos <- ggplot2::ggplot(
    input_residue_data,
    ggplot2::aes(
      y = ygreki / group_count,
      x = binned_lengths,
      fill = prediction
    )
  ) +
    ggplot2::geom_col(stat = "identity", position = "stack") +
    ggplot2::facet_wrap(~prediction, ncol = 1) +
    ggplot2::scale_fill_manual(values = c("#3a424f", "#50a675", "#b0bdd4")) +
    ggplot2::scale_x_continuous(limits = c(0, max_length)) +
    ggplot2::labs(
      title = "Distribution of non-A residues in poly(A) tails (normalized)",
      subtitle = paste0(
        "Positions were binned every 5 nucleotides (",
        subtitle_info,
        ")"
      ),
      y = "normalized non-A frequency",
      x = "poly(A) length",
      tag = "D",
      fill = "non-A residue"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous()

  rel_density_plot <- ggplot2::ggplot(
    data = input_residue_data,
    ggplot2::aes(
      x = forcats::fct_rev(prediction),
      y = est_nonA_pos,
      fill = prediction
    )
  )

  # relative density plot (fisheye + rug):
  rel_density_plot <- rel_density_plot +
    ggdist::stat_halfeye(
      adjust = .5,
      width = .6,
      .width = 0,
      justification = -.1,
      point_colour = NA
    ) +
    gghalves::geom_half_point(
      side = "l",
      shape = 124,
      range_scale = 0,
      size = 6,
      alpha = .3,
      show.legend = FALSE,
      color = "#1d2e3b"
    ) +
    ggplot2::coord_cartesian(xlim = c(1.2, NA), clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("#3a424f", "#50a675", "#b0bdd4")) +
    ggplot2::scale_y_continuous(limits = c(0, max_length)) +
    ggplot2::scale_x_discrete(
      breaks = input_residue_data$prediction,
      labels = input_residue_data$label
    ) +

    ggplot2::labs(
      title = "Distribution of non-A residues in poly(A) tails (raw positions)",
      subtitle = paste0(subtitle_info),
      y = "estimated position in poly(A) tail",
      x = "residue",
      tag = "E",
      fill = "non-A residue"
    ) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )

  # ADD OPTIONS TO PLOTS
  ##############################################################################

  if (type == "default") {
    #distribution plot
    distrib_plot <- distrib_plot
  } else if (type == "moderna") {
    #distribution plot
    distrib_plot <- distrib_plot +
      ggplot2::geom_vline(
        xintercept = 100,
        color = "red",
        size = 5,
        alpha = 0.2
      ) +
      ggplot2::labs(caption = "default Moderna pentamer position marked in red")

    rel_density_plot <- rel_density_plot +
      ggplot2::geom_hline(
        yintercept = 100,
        color = "red",
        size = 5,
        alpha = 0.2
      ) +
      ggplot2::labs(caption = "default Moderna pentamer position marked in red")
  } else {
    stop(
      "Unknown type of data defined. Please provide correct parameter",
      .call = FALSE
    )
  }

  # FINAL ASSEMBLY OF PLOT PANELS
  ##############################################################################

  # layout for patchwork
  design = "
  AB
  CC
  DD
  EE
  "

  # final plot organisation:
  final <- general_read_categories +
    residue_counts +
    distrib_plot +
    binned_length_pos +
    rel_density_plot +
    patchwork::plot_layout(design = design, heights = c(2, 4, 7, 7))

  #return(rel_density_plot)
  return(final)
}


#' Scatterplot of nonA residue positions within poly(A) tail
#'
#' This function allows to produce the scatterplot of raw non-A
#' residue predictions (y-axis) along the user-defined tail span
#' (x-axis). The resulting plot also contains the rugs and ridges
#' along the axes, which are a visual aid designed to better
#' show the distributions of given residue.
#'
#' @param residue_data A dataframe or tibble containig non-A residue predictions
#' made by ninetails pipeline.
#'
#' @param base character. One of the following ["C"/"G","U"]. This parameter
#' defines which nonadenosine is to be plotted. It is obligatory to prevent
#' the overplotting.
#'
#' @param max_length numeric [1]. This parameter controls maximum length
#' of the poly(A) tail to be taken into consideration for plotting.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' ninetails::plot_rug_density(residue_data=residue_data, base="C", max_length=100)
#' }
#'
plot_rug_density <- function(residue_data,
                             base,
                             max_length) {
  #assertions
  if (base == "C") {
    colorscale <- ggplot2::scale_fill_manual(values = c("#3a424f"))
    colorfill <- ggplot2::scale_color_manual(values = c("#3a424f"))
  } else if (base == "G") {
    colorscale <- ggplot2::scale_fill_manual(values = c("#50a675"))
    colorfill <- ggplot2::scale_color_manual(values = c("#50a675"))
  } else if (base == "U") {
    colorscale <- ggplot2::scale_fill_manual(values = c("#b0bdd4"))
    colorfill <- ggplot2::scale_color_manual(values = c("#b0bdd4"))
  } else {
    stop("Wrong base. Base must be either C, G or U")
  }

  data <- residue_data %>% dplyr::filter(prediction == base)

  #main plot
  pmain <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = polya_length,
      y = polya_length - est_nonA_pos,
      color = prediction
    )
  ) +
    ggplot2::geom_point(show.legend = FALSE) +
    ggplot2::geom_rug(
      sides = "tr",
      alpha = 0.2,
      size = 1.5,
      col = "#3a424f",
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(expand = c(0.1, 0.1)) +
    ggplot2::scale_x_continuous(expand = c(0.1, 0.1)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(size = 0.45, color = "#3a424f"),
      axis.line.y = ggplot2::element_line(size = 0.45, color = "#3a424f")
    ) +
    ggplot2::xlim(0, max_length) +
    ggplot2::ylim(0, max_length) +
    colorscale +
    colorfill +
    ggplot2::labs(
      x = "poly(A) length [nt]",
      y = "estimated non-A positions [nt]"
    )

  # Marginal densities along x axis
  xdens <- cowplot::axis_canvas(pmain, axis = "x") +
    ggplot2::geom_density(
      data = data,
      ggplot2::aes(x = polya_length, fill = prediction),
      alpha = 0.5,
      size = 0.2,
      show.legend = FALSE
    ) +
    colorscale +
    colorfill
  # Marginal densities along y axis
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  ydens <- cowplot::axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
    ggplot2::geom_density(
      data = data,
      ggplot2::aes(x = polya_length - est_nonA_pos, fill = prediction),
      alpha = 0.5,
      size = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::coord_flip() +
    colorscale +
    colorfill

  p1 <- cowplot::insert_xaxis_grob(
    pmain,
    xdens,
    grid::unit(.2, "null"),
    position = "top"
  )
  p2 <- cowplot::insert_yaxis_grob(
    p1,
    ydens,
    grid::unit(.2, "null"),
    position = "right"
  )
  p3 <- cowplot::ggdraw(p2)

  return(p3)
}

#' Plot abundances of reads with given amount of non-A residues per read
#'
#' This function plots frequencies of reads containing one, two or more separate
#' instances (occurrences) of non-As reported by ninetails. The frequency
#' is computed with respect to the total amount of decorated reads in the
#' analyzed dataset.
#'
#' @param residue_data A dataframe or tibble containig non-A residue predictions
#' made by ninetails pipeline
#'
#' @param grouping_factor character string. A grouping variable (e.g. "sample_name")
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' plt <- ninetails::plot_nonA_abundance(residue_data = residue_data,
#'                                       grouping_factor = "sample_name")
#' plt
#' }
plot_nonA_abundance <- function(residue_data, grouping_factor = NA) {

  #assertions
  if (missing(residue_data)) {
    stop(
      "Residue_data is missing. Please provide a valid residue_data argument",
      call. = FALSE
    )
  }

  if (!is.data.frame(residue_data) || nrow(residue_data) == 0) {
    stop(
      "Empty data frame provided as an input (residue_data). Please provide valid input"
    )
  }

  if (!is.na(grouping_factor)) {
    assert_condition(
      grouping_factor %in% colnames(residue_data),
      paste0(grouping_factor, " is not a column of input dataset")
    )
  }

  nonA_counts <- ninetails::count_nonA_abundance(
    residue_data = residue_data,
    grouping_factor = grouping_factor
  )

  nonA_counts <- nonA_counts %>%
    tidyr::pivot_wider(names_from = instances, values_from = count) %>%
    dplyr::mutate(total = single + two + more) %>%
    dplyr::mutate_all(~ replace_na(., 0)) %>%
    dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
    dplyr::mutate(
      psingle = single / total,
      ptwo = two / total,
      pmore = more / total
    ) %>%
    dplyr::select(!!rlang::sym(grouping_factor), psingle, ptwo, pmore) %>%
    dplyr::rename(single = psingle, two = ptwo, more = pmore) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(
      cols = c(single, two, more),
      names_to = "instances",
      values_to = "value"
    ) %>%
    dplyr::distinct()

  nonA_counts$instances <- factor(
    nonA_counts$instances,
    levels = c("single", "two", "more")
  )

  tp <- ggplot2::ggplot(
    nonA_counts,
    ggplot2::aes(x = !!rlang::sym(grouping_factor), y = value, fill = instances)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(y = "frequency", fill = "number of occurrences per read") +
    ggplot2::scale_fill_manual(
      values = c("single" = "#91b7db", "two" = "#0978e3", "more" = "#0f304f")
    ) +
    ggplot2::theme_bw()

  return(tp)
}
