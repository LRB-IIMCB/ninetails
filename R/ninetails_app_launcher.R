#' Launch the Ninetails Analysis Dashboard
#'
#' Opens an interactive Shiny application for exploring ninetails results.
#' Supports two modes:
#'
#' \strong{Multi-sample mode} (\code{config}): Loads a YAML configuration file
#' pointing to multiple samples with class and residue data.
#'
#' \strong{Single-sample mode}: Provide individual file paths. All tabs are
#' available when \code{class_file} and \code{residue_file} are supplied.
#' The Signal Viewer tab requires \code{summary_file} and \code{pod5_dir}.
#'
#' @param config Character string (optional). Path to a YAML configuration
#'   file. When provided, single-sample arguments are ignored.
#'
#' @param summary_file Character string (optional). Path to a dorado summary
#'   file. Required for Signal Viewer tab in single-sample mode.
#'
#' @param pod5_dir Character string (optional). Path to the directory containing
#'   POD5 files. Required for Signal Viewer tab in single-sample mode.
#'
#' @param class_file Character string (optional). Path to the read_classes
#'   output file from ninetails. Enables analysis tabs in single-sample mode.
#'
#' @param residue_file Character string (optional). Path to a nonadenosine_residues
#'   output file from ninetails.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @details
#' The YAML configuration file should have the following structure:
#'
#' \preformatted{
#' samples:
#'   sample_label_1:
#'     sample_name: KO_rep1
#'     group: KO
#'     class_path: /path/to/read_classes.txt
#'     residue_path: /path/to/nonadenosine_residues.txt
#'     dorado_summary: /path/to/dorado_summary.txt   # optional
#'     pod5_dir: /path/to/pod5/                      # optional
#' }
#'
#' @return Launches a Shiny application (does not return a value).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Multi-sample mode
#' ninetails::launch_signal_browser(config = "/path/to/config.yml")
#'
#' # Single-sample: all tabs
#' ninetails::launch_signal_browser(
#'   summary_file = "/path/to/dorado_summary.txt",
#'   pod5_dir     = "/path/to/pod5/",
#'   class_file   = "/path/to/read_classes.txt",
#'   residue_file = "/path/to/nonadenosine_residues.txt"
#' )
#'
#' # Single-sample: signal viewer only
#' ninetails::launch_signal_browser(
#'   summary_file = "/path/to/dorado_summary.txt",
#'   pod5_dir     = "/path/to/pod5/"
#' )
#'
#' }
#'
launch_signal_browser <- function(config = NULL,
                                  summary_file = NULL,
                                  pod5_dir = NULL,
                                  class_file = NULL,
                                  residue_file = NULL,
                                  ...) {

  for (pkg in c("shiny", "plotly", "htmltools")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' required. Install with: install.packages('", pkg, "')",
           call. = FALSE)
    }
  }

  # Initialize data containers
  class_data    <- NULL
  residue_data  <- NULL
  merged_data   <- NULL
  signal_config <- list()

  # ---- Multi-sample mode (YAML config) ----
  if (!is.null(config)) {

    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' required. Install with: install.packages('yaml')",
           call. = FALSE)
    }

    assert_condition(is_string(config), "config must be a character string")
    assert_file_exists(config, "config")

    cfg <- yaml::read_yaml(config)
    if (is.null(cfg$samples) || length(cfg$samples) == 0) {
      stop("No samples found in config file.", call. = FALSE)
    }

    class_list   <- list()
    residue_list <- list()

    cat(paste0("[", Sys.time(), "] Loading ", length(cfg$samples), " samples...\n"))

    for (sid in names(cfg$samples)) {
      s <- cfg$samples[[sid]]
      sname <- s$sample_name %||% sid
      grp   <- s$group %||% "ungrouped"

      if (!is.null(s$class_path) && file.exists(s$class_path)) {
        tryCatch({
          cd <- vroom::vroom(s$class_path, show_col_types = FALSE)
          cd$sample_name <- sname; cd$group <- grp
          class_list[[sid]] <- cd
        }, error = function(e) {
          warning("class_path for '", sid, "': ", e$message, call. = FALSE)
        })
      }

      if (!is.null(s$residue_path) && file.exists(s$residue_path)) {
        tryCatch({
          rd <- vroom::vroom(s$residue_path, show_col_types = FALSE)
          rd$sample_name <- sname; rd$group <- grp
          residue_list[[sid]] <- rd
        }, error = function(e) {
          warning("residue_path for '", sid, "': ", e$message, call. = FALSE)
        })
      }

      if (!is.null(s$dorado_summary) && !is.null(s$pod5_dir)) {
        signal_config[[sname]] <- list(
          dorado_summary = s$dorado_summary,
          pod5_dir = s$pod5_dir
        )
      }
    } # for each sample

    if (length(class_list) == 0) {
      stop("No class data loaded. Check class_path entries.", call. = FALSE)
    }

    class_data <- dplyr::bind_rows(class_list)
    class_data$sample_name <- as.factor(class_data$sample_name)
    class_data$group <- as.factor(class_data$group)

    if (length(residue_list) > 0) {
      residue_data <- dplyr::bind_rows(residue_list)
      residue_data$sample_name <- as.factor(residue_data$sample_name)
      residue_data$group <- as.factor(residue_data$group)
    }

    cat(paste0("[", Sys.time(), "] Loaded ",
               format(nrow(class_data), big.mark = ","), " reads from ",
               length(unique(class_data$sample_name)), " samples.\n"))

  } else {
    # ---- Single-sample mode ----

    # Load class data if provided
    if (!is.null(class_file)) {
      assert_condition(is_string(class_file), "class_file must be a character string")
      if (file.exists(class_file)) {
        tryCatch({
          class_data <- vroom::vroom(class_file, show_col_types = FALSE)
          cat(paste0("[", Sys.time(), "] Loaded ",
                     format(nrow(class_data), big.mark = ","),
                     " reads from class file.\n"))
        }, error = function(e) {
          warning("Error loading class_file: ", e$message, call. = FALSE)
        })
      } else {
        warning("class_file not found: ", class_file, call. = FALSE)
      }
    }

    # Load residue data if provided
    if (!is.null(residue_file)) {
      assert_condition(is_string(residue_file), "residue_file must be a character string")
      if (file.exists(residue_file)) {
        tryCatch({
          residue_data <- vroom::vroom(residue_file, show_col_types = FALSE)
        }, error = function(e) {
          warning("Error loading residue_file: ", e$message, call. = FALSE)
        })
      } else {
        warning("residue_file not found: ", residue_file, call. = FALSE)
      }
    }

    # Signal config for single sample
    if (!is.null(summary_file) && !is.null(pod5_dir)) {
      signal_config[["single"]] <- list(
        dorado_summary = summary_file,
        pod5_dir = pod5_dir
      )
    }
  }

  # Attempt merged table
  if (!is.null(class_data) && !is.null(residue_data) && nrow(residue_data) > 0) {
    tryCatch({
      merged_data <- ninetails::merge_nonA_tables(
        class_data, residue_data, pass_only = FALSE
      )
    }, error = function(e) {
      warning("Could not create merged table: ", e$message, call. = FALSE)
    })
  }

  # Pass everything via shinyOptions
  shiny::shinyOptions(
    ninetails.class_data    = class_data,
    ninetails.residue_data  = residue_data,
    ninetails.merged_data   = merged_data,
    ninetails.signal_config = signal_config,
    ninetails.summary_file  = summary_file %||% "",
    ninetails.pod5_dir      = pod5_dir %||% "",
    ninetails.residue_file  = residue_file %||% ""
  )

  app_dir <- system.file("app", package = "ninetails")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("App not found. Expected at: inst/app/", call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
