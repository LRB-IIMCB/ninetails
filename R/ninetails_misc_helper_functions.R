################################################################################
# ASSERTION HELPERS (base R replacements for checkmate/assertthat)
################################################################################
#' Assert condition is TRUE, stop with message if FALSE
#'
#' Base R replacement for assertthat::assert_that
#'
#' @param cond Logical condition to check
#' @param msg Error message if condition is FALSE
#'
#' @return Invisible TRUE if assertion passes
#' @keywords internal
#'
assert_condition <- function(cond, msg = "Assertion failed") {
  if (!isTRUE(cond)) {
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

#' Test if x is a single non-empty character string
#'
#' Base R replacement for checkmate::test_string
#'
#' @param x Object to test
#' @param min.chars Minimum number of characters (default 1)
#' @param null.ok Allow NULL? (default FALSE)
#'
#' @return TRUE if x is a valid string, FALSE otherwise
#' @keywords internal
#'
is_string <- function(x, min.chars = 1, null.ok = FALSE) {
  if (is.null(x)) return(null.ok)
  is.character(x) && length(x) == 1 && !is.na(x) && nchar(x) >= min.chars
}

#' Assert file exists with informative error
#'
#' Base R replacement for checkmate::assert_file_exists
#'
#' @param x File path to check
#' @param arg_name Optional argument name for clearer error message
#'
#' @return Invisible TRUE if file exists
#' @keywords internal
#'
assert_file_exists <- function(x, arg_name = NULL) {
  if (!file.exists(x)) {
    if (is.null(arg_name)) {
      stop("File does not exist: ", x, call. = FALSE)
    } else {
      stop(arg_name, " file does not exist: ", x, call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Assert directory exists with informative error
#'
#' @param x Directory path to check
#' @param arg_name Optional argument name for clearer error message
#'
#' @return Invisible TRUE if directory exists
#' @keywords internal
#'
assert_dir_exists <- function(x, arg_name = NULL) {
  if (!dir.exists(x)) {
    if (is.null(arg_name)) {
      stop("Directory does not exist: ", x, call. = FALSE)
    } else {
      stop(arg_name, " directory does not exist: ", x, call. = FALSE)
    }
  }
  invisible(TRUE)
}


#' Check for no NA values
#'
#' Base R replacement for assertthat::noNA
#'
#' @param x Vector to check
#'
#' @return TRUE if no NA values present, FALSE otherwise
#' @keywords internal
#'
no_na <- function(x) {
  !anyNA(x)
}


################################################################################
# SIGNAL PROCESSING HELPERS
################################################################################
#' Winsorizes nanopore signal.
#'
#' Removes high cliffs (extending above & below
#' signal range). Slightly affects signal extremes.
#'
#' Based on the A. Signorell's DescTools
#' https://github.com/AndriSignorell/DescTools/blob/master/R/DescTools.r
#'
#'
#' @param signal numeric vector. A vector corresponding to given ONT signal.
#'
#' @return winsorized signal (numeric vector).
#' @export
#'
#' @examples
#' \dontrun{
#'
#' winsorized_signal <- winsorize_signal(signal= sample(200:300))
#'
#' }
#'
winsorize_signal <- function(signal){

  #assertions

  if (missing(signal)) {
    stop("Signal vector is missing. Please provide a valid signal argument.", .call = FALSE)
  }

  assert_condition(is.numeric(signal),
                   "Signal vector must be numeric. Please provide a valid argument.")
  assert_condition(is.atomic(signal),
                   "Signal vector must be atomic. Please provide a valid argument.")
  assert_condition(no_na(signal),
                   "Signal vector must not contain missing values. Please provide a valid argument.")

  signal_q <- stats::quantile(x=signal, probs=c(0.005, 0.995), na.rm=TRUE, type=7)
  minimal_val <- signal_q[1L]
  maximal_val <- signal_q[2L]

  signal[signal<minimal_val] <- minimal_val
  signal[signal>maximal_val] <- maximal_val

  winsorized_signal <- as.integer(signal)

  return(winsorized_signal)

}


#' Loads keras model for multiclass signal prediction.
#'
#' @param keras_model_path either missing or character string. Full path of the
#' .h5 file with keras model used to predict signal classes. If function is
#' called without this argument(argument is missing) the default pretrained
#' model will be loaded. Otherwise, the dir with custom model shall be provided.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' load_keras_model(keras_model_path = "/path/to/the/model/in/hdf5_format")
#' }
#'
load_keras_model <- function(keras_model_path){
  if (rlang::is_missing(keras_model_path)) {
    path_to_default_model <- system.file("extdata", "cnn_model", "gasf_gadf_combined_model_20220808.h5", package="ninetails")
    keras_model <- keras::load_model_hdf5(path_to_default_model)
  } else {
    keras_model <- keras::load_model_hdf5(keras_model_path)
  }
}

#' Check if fast5 file is multi-read format
#'
#' Checks whether the provided fast5 file structure corresponds to a multi-read
#' fast5 file format. Multi-read fast5 files contain multiple reads stored
#' in groups named with "read_" prefix.
#'
#' @param fast5_file_structure A data frame returned by rhdf5::h5ls() containing
#'   the HDF5 structure of a fast5 file. Must contain a 'name' column.
#'
#' @return Logical. TRUE if the file is multi-read fast5 format, FALSE otherwise.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Check if a fast5 file is multi-read format
#' fast5_path <- "path/to/file.fast5"
#' file_structure <- rhdf5::h5ls(fast5_path, recursive = FALSE)
#' is_multifast5(file_structure)
#' }
#'
is_multifast5 <- function(fast5_file_structure) {

  assert_condition(is.data.frame(fast5_file_structure),
                   "fast5_file_structure must be a data frame (output of rhdf5::h5ls())")
  assert_condition("name" %in% colnames(fast5_file_structure),
                   "fast5_file_structure must contain 'name' column")

  sum(grepl('read_', fast5_file_structure$name)) > 0
}


#' Check if fast5 file contains RNA reads
#'
#' Checks whether the provided fast5 file contains RNA sequencing reads
#' by examining the experiment_type attribute in the context_tags group.
#'
#' @param fast5_file Character string. Path to the fast5 file.
#' @param read_id Character string. The read identifier within the fast5 file
#'   (e.g., "read_abc123"). Can be obtained from rhdf5::h5ls() output.
#'
#' @return Logical. TRUE if the file contains RNA reads, FALSE otherwise.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Check if a fast5 file contains RNA reads
#' fast5_path <- "path/to/file.fast5"
#' file_structure <- rhdf5::h5ls(fast5_path, recursive = FALSE)
#' read_id <- file_structure$name[1]
#' is_RNA(fast5_path, read_id)
#' }
#'
is_RNA <- function(fast5_file, read_id) {

  assert_condition(is_string(fast5_file),
                   "fast5_file must be a character string path")
  assert_file_exists(fast5_file, "fast5_file")
  assert_condition(is_string(read_id),
                   "read_id must be a character string")

  read_context_tags <- rhdf5::h5readAttributes(fast5_file, paste0(read_id, "/context_tags"))
  read_context_tags$experiment_type == "rna"
}


#' Checks if the provided directory contains fast5 files in the correct format.
#'
#' This function analyses the structure of the first fast5 file in the given
#' directory and checks whether it fulfills the analysis requirements (if the
#' file is multifast5, basecalled by Guppy basecaller and containing provided
#' basecall_group). Otherwise the function throws an error (with description).
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#'
#' @return outputs the text info with basic characteristics of the data.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' check_fast5_filetype <- function(workspace = '/path/to/guppy/workspace',
#'                                  basecalled_group = 'Basecall_1D_000')
#'
#' }
#'
#'
#'
# This lookup function is inspired by adnaniazi's explore-basecaller-and-fast5type.R from tailfindr
# https://github.com/adnaniazi/tailfindr/blob/master/R/explore-basecaller-and-fast5type.R

check_fast5_filetype <- function(workspace,
                                 basecall_group){

  #Assertions
  if (missing(workspace)) {
    stop("Directory with basecalled fast5s is missing. Please provide a valid workspace argument.", call. =FALSE)
  }

  if (missing(basecall_group)) {
    stop("Basecall group is missing. Please provide a valid basecall_group argument.", call. =FALSE)
  }

  assert_condition(is.character(workspace),
                   "Path to fast5 files is not a character string. Please provide valid path to basecalled fast5 files.")


  #list fast5 files in given dir
  fast5_files_list <- list.files(path = workspace, pattern = "\\.fast5$", recursive = TRUE, full.names = TRUE)

  #count fast5 files
  num_fast5_files <- length(fast5_files_list)

  cat(paste0('[', as.character(Sys.time()), '] ','Found ', num_fast5_files, ' fast5 file(s) in provided directory.\n'))

  #closer look into the first file on the list
  selected_fast5_file <-fast5_files_list[1]
  selected_fast5_file_structure <- rhdf5::h5ls(file.path(selected_fast5_file), recursive = FALSE)

  selected_fast5_read <- selected_fast5_file_structure$name[1]

  cat(paste0('[', as.character(Sys.time()), '] ','Analyzing one of the given fast5 files to check', '\n','if the data are in required format... \n'))


  # checking whether fast5 file is single or multi
  if (!is_multifast5(selected_fast5_file_structure)) {
    stop("The provided fast5 is single fast5 file. Please provide multifast5 file(s).", call. = FALSE)
  }

  # check whether file is basecalled or not
  tryCatch(selected_basecall_group <- rhdf5::h5read(selected_fast5_file,paste0(selected_fast5_read,"/Analyses/", basecall_group)), error = function(e) { cat("The previewed fast5 file does not contain defined basecall_group. Check basecall_group - is it valid? Ninetails requires fast5 files basecalled by Guppy.") })

  if (exists('selected_basecall_group')) {
    # checking whether the fast5 file contains RNA ONT reads
    if (!is_RNA(selected_fast5_file, selected_fast5_read)) {
      stop("The provided fast5 does not contain RNA reads. Please provide multifast5 file(s) with RNA reads.", call. = FALSE)
    }


    # retrieve basecaller & basecalling model (read attributes)
    selected_basecall_group <- rhdf5::h5readAttributes(selected_fast5_file,paste0(selected_fast5_read,"/Analyses/", basecall_group))
    basecaller_used <- selected_basecall_group$name
    model_used <- selected_basecall_group$model_type

    # retrieve guppy basecaller version (read attributes)
    path_to_guppy_version <- rhdf5::h5readAttributes(selected_fast5_file,paste0(selected_fast5_read,"/tracking_id"))
    guppy_version <- path_to_guppy_version$guppy_version

    # close all handled instances (groups, attrs) of fast5 file
    rhdf5::h5closeAll()

    cat('  Previewed fast5 file parameters:\n')
    cat('    data type: RNA \n')
    cat('    fast5 file type: multifast5 \n')
    cat('    basecaller used:',basecaller_used,' \n')
    cat('    basecaller version:',guppy_version,' \n')
    cat('    basecalling model:',model_used,' \n',' \n')

  } else {
    # close all handled instances (groups, attrs) of fast5 file
    rhdf5::h5closeAll()
    stop(paste0('[', as.character(Sys.time()), '] ','Ninetails encountered an error. Please provide fast5 files basecalled by Guppy software. '), call. =FALSE)
  }

}




#' Substitutes 0s surrounded by adjacent nonzeros in pseudomove vector
#' to facilitate position-centering function.
#'
#' This function helps to avoid redundancy introduced by z-score thesholding
#' algo (this happens when signal is jagged), so one segment would be reported
#' instead of multiple.
#'
#' @param pseudomoves numeric vector produced by z-score
#' filtering algo (\code{\link{filter_signal_by_threshold}}) corresponding
#' to the tail region of the read of interest as delimited by
#' nanopolish polya function.
#'
#' @return a numeric vector of adjusted pseudomoves (smoothened)
#' @export
#'
#' @examples
#' \dontrun{
#'
#' substitute_gaps(pseudomoves = pseudomoves_vector)
#'
#'}
#'
substitute_gaps <- function(pseudomoves){
  rle_pseudomoves <- rle(pseudomoves)
  indx <- rle_pseudomoves$lengths < 2 & rle_pseudomoves$values == 0 & c(Inf, rle_pseudomoves$values[-length(rle_pseudomoves$values)]) == c(rle_pseudomoves$values[-1], Inf)
  if (any(indx)) rle_pseudomoves$values[indx] <- rle_pseudomoves$values[which(indx)-1]
  adjusted_pseudomoves <- inverse.rle(rle_pseudomoves)

  return(adjusted_pseudomoves)
}

################################################################################
# ADDITIONAL DATA PROCESSING HELPERS
################################################################################
#' Calculates statistical mode of given vector.
#'
#' This function operates in 2 flavours. Either of them can be defined
#' by the "method" parameter. There are 2 options available: "density"
#' and "value".
#'
#' If the "density" (default) option was chosen, then the
#' statistical mode (most frequent value) would be computed from the
#' normalized data density distribution. In this mode, the functon returns
#' a single value.
#'
#' If the "value" option was selected, then the function return the mode based
#' on actual value. This means that, if the dataset was bi- or multimodal, the
#' returned vector would contain more than one most frequent values.
#'
#' This function was written based on the following SO thread:
#' https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
#' Special thanks to Chris and hugovdberg!
#'
#' @param x vector with values [numeric] to compute the mode for
#'
#' @param method [character] string; which method of computing statistical mode
#' is meant to be used by the function. 2 options available: "value"/"density".
#' The latter is set by default.
#'
#' @param na.rm logical [TRUE/FALSE] whether or not to remove NA values from the
#' calculation. Set to FALSE by default.
#'
#' @return statistical mode of given vector of values [numeric]
#' @export
#'
#' @examples
#'\dontrun{
#'
#'test1 <- c(rep(2,5), rep(3,4), rep(1,4), rep(8,2), rep(7,3), rep(5,3)) # 1 most freq val
#'test2 <- c(rep(2,5), rep(3,4), rep(1,4), rep(8,5), rep(7,3), rep(5,3)) # 2 most freq val
#'test3 <- c(rep(2,5), rep(3,4), rep(1,4), rep(8,5), rep(7,5), rep(5,3)) # 3 most freq val
#'test4 <- c("Lorem", "ipsum", "dolor", "sit", "amet")
#'
#'result <- ninetails::get_mode(x=test1, method= "value")
#'print(result)
#'class(result)
#'
#'}
#'
get_mode <- function(x, method ="density", na.rm = FALSE) {

  # assertions
  assert_condition(is.numeric(x),
                   "Provided vector must be numeric. Please provide a valid argument.")

  x <- unlist(x)
  if (na.rm) {
    x <- x[!is.na(x)]
  }

  if (method %in% c("value", "density", "") | is.na(method)) {
    # Return actual mode (from the real values in dataset)
    if (method %in% c("density", "")) {

      # Return density mode for normalized data - only for numeric!)
      d <- stats::density(x)
      return(d$x[d$y==max(d$y)])
      #return(modeest::mlv(x,na.rm=TRUE,method="parzen", abc=T)) #in some cases this method produces "weird" output

    } else if (method %in% c("value")) {

      uniqx <- unique(x)
      n <- length(uniqx)
      freqx <- tabulate(match(x, uniqx))
      modex <- freqx == max(freqx)
      return(uniqx[which(modex)])
    }
  } else {
    stop("Wrong mode provided. Please provide either 'density' or 'value'", call. =FALSE)
  }
}

#'  Correcting class recognition in ninetails data
#'
#'  This helper fnction solves issues with backcompatibility of ninetails'
#'  read class naming convention.
#'
#'  Previously, the read classes were named as "modified", "unmodified" and
#'  "unclassified" (pre-release versions). However, it has been changed
#'  to more precise terms: "decorated", "blank" and "unclassified", since the
#'  presence of non-A within poly(A) tail does not necessary have to be
#'  the result of post-transcriptional, cytoplasmic modification. It might be
#'  caused by the semi-templated tails or polymerase slippage during RNA
#'  synthesis/canonical tailing or caused by other potential yet undefined
#'  mechanisms.
#'
#' @param df data frame with ninetails read classification results (class_data)
#' or the output of \code{\link{merge_nonA_tables}} function.
#'
#' @return corrected dataframe (with changed class labels)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' class_data <- ninetails::correct_labels(class_data)
#' }
#'
correct_labels <- function(df) {

  assert_condition(is.data.frame(df),
                   "Provided input must be dataframe. Please provide a valid argument.")

  if("class" %in% colnames(df)) {
    df$class <- ifelse(df$class %in% c("unmodified", "unclassified"),
                       "blank",
                       ifelse(df$class == "blank",
                              "unclassified",
                              "decorated"))
  }

  return(df)
}

################################################################################
# REQUIRED BY DORADO-DEPENDENT PIPELINES
################################################################################
#' Filter Dorado summary file for reads fulfilling ninetails quality criteria
#'
#' This function takes a Dorado summary file (from ONT Dorado basecaller)
#' or a data frame containing equivalent summary information,
#' and filters out reads that do not meet ninetails quality and alignment criteria.
#' Specifically, it removes reads with missing alignments, low mapping quality,
#' poly(A) start positions indicating most likely artifacts,
#' or poly(A) tails shorter than 10 bases.
#'
#' @param dorado_summary Character path to a Dorado summary file (tab-delimited),
#' or a data frame containing summary information.
#'
#' @returns A filtered tibble or data frame containing only reads that meet
#' the filtering criteria.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From file
#' filtered <- filter_dorado_summary("dorado_summary.txt")
#'
#' # From data frame
#' df <- data.frame(
#'   read_id = c("read1", "read2"),
#'   alignment_direction = c("+", "*"),
#'   alignment_mapq = c(60, 0),
#'   poly_tail_start = c(100, 0),
#'   poly_tail_length = c(20, 5)
#' )
#' filtered <- filter_dorado_summary(df)
#' }
filter_dorado_summary <- function(dorado_summary){

  # Variable binding (suppressing R CMD check from throwing an error)
  alignment_direction <- alignment_mapq <- poly_tail_start <- poly_tail_length <- NULL

  # If input is a file path, read with vroom
  if (is.character(dorado_summary)) {
    dorado_summary <- vroom::vroom(dorado_summary, delim = "\t")
  }
  # Ensure it's a data.frame/tibble
  if (!is.data.frame(dorado_summary)) {
    stop("Input must be a data.frame, tibble, or a valid file path.")
  }
  required_cols <- c("read_id", "poly_tail_length")
  missing_cols <- setdiff(required_cols, names(dorado_summary))
  if (length(missing_cols) > 0) {
    stop(sprintf("Required columns missing: %s", paste(missing_cols, collapse = ", ")))
  }

  #These are the reads that fulfill the quality criteria for classification by Ninetails:
  #(I) they are mapped to the reference (the alignment has a direction, either + or –, not a placeholder *),
  #(II) the mapping quality is greater than 0,
  #(III) poly(A) tail coordinates are correctly indicated (the poly(A) start
  #cannot be 0 in the signal because, in DRS, the adapter region passes through the pore first),
  #(IV) the poly(A) tail must be at least 10 nt long (shorter tails cannot be reliably processed by the CNN).

  dorado_summary_filtered <- dorado_summary %>%
    dplyr::filter(alignment_direction!="*" &
                  alignment_mapq!=0 &
                  poly_tail_start!=0 &
                  poly_tail_length>=10)


  return(dorado_summary_filtered)

}



#' Check and handle existing output directory for ninetails analysis
#'
#' This function checks if the specified output directory already exists and contains
#' files that might be overwritten by the ninetails analysis. If the directory exists
#' and is not empty, it prompts the user for action and logs the decision.
#'
#' This function is not intended to be used outside the pipeline wrapper.
#'
#' @section User Interaction:
#' When an existing non-empty directory is detected, the function will:
#' \itemize{
#'   \item Display the directory path and file count
#'   \item Prompt the user to choose: abort analysis or overwrite existing files
#'   \item Wait for user input (a/A for abort, o/O for overwrite)
#'   \item Log the user's decision and proceed accordingly
#' }
#'
#' @section Directory States:
#' The function handles several directory states:
#' \itemize{
#'   \item \strong{Non-existent}: Creates the directory and logs creation
#'   \item \strong{Empty}: Uses existing directory and logs confirmation
#'   \item \strong{Non-empty}: Prompts user for overwrite decision
#' }
#'
#' @param save_dir Character string. Full path to the output directory where
#'   ninetails results will be saved.
#' @param log_message Function for logging messages to both console and log file.
#'   Should accept parameters: message, type, section.
#'
#' @returns Logical. Returns TRUE if the analysis should proceed, FALSE if the
#'   user chose to abort the analysis.
#'
#' @section Implementation Notes:
#' \itemize{
#'   \item Uses \code{readline()} for interactive user input
#'   \item Logs all decisions and actions for audit trail
#'   \item Handles edge cases like permission errors
#'   \item Validates user input with retry mechanism
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with logging function
#' log_func <- function(msg, type = "INFO", section = NULL) {
#'   cat(sprintf("[%s] %s\n", type, msg))
#' }
#'
#' # Check directory and proceed if allowed
#' should_proceed <- check_output_directory("/path/to/output", log_func)
#' if (should_proceed) {
#'   # Continue with analysis
#' } else {
#'   # User chose to abort
#' }
#' }
check_output_directory <- function(save_dir, log_message) {

  # Check if directory exists
  if (!dir.exists(save_dir)) {
    # Directory doesn't exist, create it
    tryCatch({
      dir.create(save_dir, recursive = TRUE)
      log_message(sprintf("Created output directory: %s", save_dir), "INFO", "Directory Setup")
      return(TRUE)
    }, error = function(e) {
      log_message(sprintf("Failed to create output directory: %s", e$message), "ERROR", "Directory Setup")
      stop(sprintf("Cannot create output directory: %s", e$message), call. = FALSE)
    })
  }

  # Directory exists, check if it's empty
  existing_files <- list.files(save_dir, recursive = TRUE, include.dirs = FALSE)

  if (length(existing_files) == 0) {
    # Directory is empty, safe to proceed
    log_message(sprintf("Using existing empty directory: %s", save_dir), "INFO", "Directory Setup")
    return(TRUE)
  }

  # Directory contains files, prompt user
  log_message(sprintf("Output directory already exists and contains %d files", length(existing_files)),
              "WARNING", "Directory Conflict")
  log_message(sprintf("Directory path: %s", save_dir), "WARNING")

  # Display warning to user
  cat("\n")
  cli::cli_alert_warning("Output directory already exists and is not empty!")
  cli::cli_text(sprintf("Directory: {.path %s}", save_dir))
  cli::cli_text(sprintf("Contains: {.val %d} files", length(existing_files)))
  cat("\n")

  # Show some example files if there are many
  if (length(existing_files) <= 5) {
    cli::cli_text("Existing files:")
    for (file in existing_files) {
      cli::cli_li(file)
    }
  } else {
    cli::cli_text("Example existing files:")
    for (file in existing_files[1:3]) {
      cli::cli_li(file)
    }
    cli::cli_text(sprintf("... and %d more files", length(existing_files) - 3))
  }

  cat("\n")
  cli::cli_text("Choose an action:")
  cli::cli_li("{.strong a} = Abort analysis (recommended if unsure)")
  cli::cli_li("{.strong c} = Continue and overwrite some of existing files")
  cat("\n")

  # Get user input with validation
  while (TRUE) {
    user_input <- readline(prompt = "Enter your choice (a/c): ")
    user_input <- tolower(trimws(user_input))

    if (user_input %in% c("a", "abort")) {
      log_message("User chose to abort analysis due to existing files", "INFO", "User Decision")
      cli::cli_alert_info("Analysis aborted by user.")
      return(FALSE)

    } else if (user_input %in% c("c", "continue")) {
      log_message("User chose to continue and overwrite existing files", "WARNING", "User Decision")
      log_message(sprintf("Proceeding with analysis in directory: %s", save_dir), "WARNING")
      cli::cli_alert_warning("Proceeding with analysis. Existing files may be overwritten.")
      return(TRUE)

    } else {
      cli::cli_alert_danger("Invalid input. Please enter 'a' to abort or 'c' to continue.")
    }
  }
}

################################################################################
# cDNA PIPELINE HELPERS
################################################################################
#' Count trailing occurrences of a character in a string
#'
#' This helper function counts how many times a specific character appears
#' consecutively at the end of a string, working backwards from the last character.
#'
#' @param string Character string to analyze
#' @param char Single character to count at the end of the string
#'
#' @return Integer count of trailing character occurrences
#' @export
#'
#' @examples
#' \dontrun{
#' count_trailing_chars("ACGTTTTT", "T")  # Returns 4
#' count_trailing_chars("ACGT", "A")      # Returns 0
#' }
count_trailing_chars <- function(string, char) {

  # assertions
  assert_condition(is_string(string),
                   "string must be a character string")
  assert_condition(is_string(char) && nchar(char) == 1,
                   "char must be a single character")

  count <- 0
  len <- nchar(string)
  for (i in len:1) {
    if (substr(string, i, i) == char) {
      count <- count + 1
    } else {
      break
    }
  }
  return(count)
}

#' Generate reverse complement of a DNA sequence
#'
#' This function creates the reverse complement of a DNA sequence string,
#' handling standard nucleotide bases (A, T, G, C) and ambiguous bases (N).
#'
#' @param sequence Character string containing DNA sequence
#'
#' @return Character string with reverse complement sequence
#' @export
#'
#' @examples
#' \dontrun{
#' reverse_complement("ATCG")    # Returns "CGAT"
#' reverse_complement("AAATTT")  # Returns "AAATTT"
#' }
reverse_complement <- function(sequence) {

  # assertions
  assert_condition(is_string(sequence),
                   "sequence must be a character string")


  complement_map <- c(A="T", T="A", G="C", C="G", N="N")
  chars <- strsplit(sequence, "")[[1]]
  comp <- sapply(chars, function(c) complement_map[c], USE.NAMES = FALSE)
  paste(rev(comp), collapse = "")
}

#' Calculate edit distance with sliding window matching
#'
#' This function implements edit distance calculation with HW mode, where
#' the query sequence can align anywhere within the target sequence using
#' a sliding window approach. It finds the minimum edit distance across
#' all possible alignments.
#'
#' @param query Character string with the primer/query sequence
#' @param target Character string with the target sequence window
#'
#' @return Integer representing the minimum edit distance found
#' @export
#'
#' @examples
#' \dontrun{
#' edit_distance_hw("ATCG", "XXATCGXX")  # Returns 0 (perfect match inside)
#' edit_distance_hw("ATCG", "ATTT")      # Returns minimum distance
#' }
edit_distance_hw <- function(query, target) {

  # assertions
  assert_condition(is_string(query),
                   "query must be a character string")
  assert_condition(is_string(target),
                   "target must be a character string")

  query_len <- nchar(query)
  target_len <- nchar(target)

  if (query_len == 0 || target_len == 0) {
    return(query_len)
  }

  # Simulate HW mode by sliding query across target
  min_dist <- .Machine$integer.max

  # Slide query across target
  if (target_len >= query_len) {
    for (i in 1:(target_len - query_len + 1)) {
      window <- substr(target, i, i + query_len - 1)
      dist <- utils::adist(query, window)[1, 1]
      min_dist <- min(min_dist, dist)
    }
  }

  # Also try full alignment
  full_dist <- adist(query, target)[1, 1]
  min_dist <- min(min_dist, full_dist)

  return(min_dist)
}
