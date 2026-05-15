#' Launch the Ninetails Analysis Dashboard (Guppy Legacy)
#'
#' Opens an interactive Shiny application for exploring ninetails results
#' from the Guppy legacy pipeline (\code{check_tails_guppy}). This is the
#' fast5-based counterpart of \code{\link{launch_signal_browser}}, designed
#' for data basecalled with Guppy \eqn{\le} 6.0.0 using Nanopolish poly(A)
#' coordinates.
#'
#' The dashboard provides the same six tabs as the Dorado version
#' (Classification, Residues, Poly(A) length, Signal Viewer, Download,
#' About), with two Guppy-specific additions:
#'
#' \describe{
#'   \item{\strong{Nanopolish QC}}{Additional plot in the Classification
#'     tab showing the distribution of Nanopolish QC tags (PASS, ADAPTER,
#'     NOREGION, SUFFCLIP, etc.) via \code{plot_nanopolish_qc()}.}
#'   \item{\strong{Fast5 Signal Viewer}}{Uses \code{plot_squiggle_fast5()}
#'     and \code{plot_tail_range_fast5()} instead of POD5-based functions.
#'     Includes a \code{moves} toggle to show/hide basecaller move
#'     transitions. No Python dependency required.}
#' }
#'
#' @param config Character string (optional). Path to a YAML configuration
#'   file defining multiple samples. When provided, single-sample arguments
#'   are ignored.
#'
#' @param nanopolish_file Character string (optional). Path to the Nanopolish
#'   polya output file. Required for the Signal Viewer tab in single-sample
#'   mode.
#'
#' @param sequencing_summary_file Character string (optional). Path to the
#'   Guppy sequencing summary file. Required for the Signal Viewer tab in
#'   single-sample mode.
#'
#' @param workspace Character string (optional). Path to the directory
#'   containing multi-fast5 files. Required for the Signal Viewer tab in
#'   single-sample mode.
#'
#' @param class_file Character string (optional). Path to the
#'   \code{read_classes} output file from \code{check_tails_guppy()}.
#'
#' @param residue_file Character string (optional). Path to a
#'   \code{nonadenosine_residues} output file from ninetails.
#'
#' @param basecall_group Character string. Fast5 hierarchy level for
#'   basecall data extraction. Default: \code{"Basecall_1D_000"}.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @details
#' The YAML configuration file should have the following structure:
#'
#' \preformatted{
#' samples:
#'   sample_label:
#'     sample_name: WT_rep1
#'     group: WT
#'     class_path: /path/to/read_classes.txt
#'     residue_path: /path/to/nonadenosine_residues.txt
#'     nanopolish: /path/to/nanopolish_output.tsv   # or polya_path
#'     sequencing_summary: /path/to/sequencing_summary.txt  # or seq_summary
#'     workspace: /path/to/fast5/                        # optional
#' }
#'
#' Signal viewer fields accept alternative names: \code{nanopolish} or
#' \code{polya_path} for the Nanopolish polya output, and
#' \code{sequencing_summary} or \code{seq_summary} for the Guppy
#' sequencing summary.
#'
#' @return Launches a Shiny application (does not return a value).
#'
#' @seealso \code{\link{launch_signal_browser}} for the Dorado DRS version,
#'   \code{\link{check_tails_guppy}} for the Guppy pipeline.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Multi-sample mode
#' ninetails::launch_signal_browser_guppy(config = "config_guppy.yml")
#'
#' # Single-sample: all tabs
#' ninetails::launch_signal_browser_guppy(
#'   nanopolish_file = "/path/to/nanopolish_output.tsv",
#'   sequencing_summary_file  = "/path/to/sequencing_summary.txt",
#'   workspace = "/path/to/fast5/",
#'   class_file = "/path/to/read_classes.txt",
#'   residue_file = "/path/to/nonadenosine_residues.txt"
#' )
#'
#' }
#'
launch_signal_browser_guppy <- function(config = NULL,
                                        nanopolish_file = NULL,
                                        sequencing_summary_file = NULL,
                                        workspace = NULL,
                                        class_file = NULL,
                                        residue_file = NULL,
                                        basecall_group = "Basecall_1D_000",
                                        ...) {

  for (pkg in c("shiny", "plotly", "htmltools")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' required. Install with: install.packages('", pkg, "')",
           call. = FALSE)
    }
  }

  # Initialize data containers
  class_data <- NULL
  residue_data <- NULL
  merged_data <- NULL
  signal_config <- list()

  #  Multi-sample mode (YAML config)
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
      sname <- if (!is.null(s$sample_name)) s$sample_name else sid
      grp   <- if (!is.null(s$group)) s$group else "ungrouped"

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

      # Guppy signal config: accept nanopolish/polya_path and sequencing_summary/seq_summary
      np_val <- if (!is.null(s$nanopolish)) s$nanopolish else if (!is.null(s$polya_path)) s$polya_path else NULL
      ss_val <- if (!is.null(s$sequencing_summary)) s$sequencing_summary else if (!is.null(s$seq_summary)) s$seq_summary else NULL
      ws_val <- s$workspace
      if (!is.null(np_val) && !is.null(ss_val) && !is.null(ws_val)) {
        signal_config[[sname]] <- list(
          nanopolish_path = np_val,
          sequencing_summary_path = ss_val,
          workspace = ws_val
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
    #  Single-sample mode

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
    if (!is.null(nanopolish_file) && !is.null(sequencing_summary_file) &&
        !is.null(workspace)) {
      signal_config[["single"]] <- list(
        nanopolish_path = nanopolish_file,
        sequencing_summary_path = sequencing_summary_file,
        workspace = workspace
      )
    }

    # Add default sample_name and group if missing
    if (!is.null(class_data)) {
      if (!"sample_name" %in% names(class_data))
        class_data$sample_name <- as.factor("sample_1")
      if (!"group" %in% names(class_data))
        class_data$group <- as.factor("ungrouped")
    }
    if (!is.null(residue_data)) {
      if (!"sample_name" %in% names(residue_data))
        residue_data$sample_name <- as.factor("sample_1")
      if (!"group" %in% names(residue_data))
        residue_data$group <- as.factor("ungrouped")
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
    ninetails.class_data = class_data,
    ninetails.residue_data  = residue_data,
    ninetails.merged_data = merged_data,
    ninetails.signal_config = signal_config,
    ninetails.basecall_group = basecall_group
  )

  app_dir <- system.file("app_guppy", package = "ninetails")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Guppy app not found. Expected at: inst/app_guppy/", call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
