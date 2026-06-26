#' Launch the Ninetails Analysis Dashboard
#'
#' Opens an interactive Shiny application for exploring and visualizing
#' ninetails poly(A) tail composition analysis results. The dashboard
#' provides a comprehensive overview of read classification, non-adenosine
#' residue composition, poly(A) tail length distributions, and raw nanopore
#' signal inspection - all in a single interface with interactive filters,
#' per-sample breakdowns, and a configurable report generator.
#'
#' The dashboard is organized into six tabs:
#'
#' \describe{
#'   \item{\strong{Classification}}{Summary value boxes (samples,
#'     transcripts, total/blank/decorated reads) and bar charts showing
#'     read classification and non-A abundance per sample or condition.
#'     Three classification views are available: summary, detailed (by
#'     comment code), and decorated-only.}
#'   \item{\strong{Residues}}{Distribution of non-adenosine residue types
#'     (C, G, U) with per-sample rug density plots showing positional
#'     distribution along poly(A) tails (subsampled to 1,000 points per
#'     residue type per sample), and a searchable summary table when
#'     merged data is available.}
#'   \item{\strong{Poly(A) length}}{Density distributions of poly(A)
#'     tail lengths with condition filtering, central tendency overlays,
#'     8 color palettes, and a summary statistics table (n, mean, median,
#'     SD, SEM) that updates with all filters.}
#'   \item{\strong{Signal Viewer}}{Raw nanopore signal visualization from
#'     POD5 files with two sub-tabs: Static Viewer (ggplot2 full signal
#'     and zoomed poly(A) region) and Dynamic Explorer (interactive Plotly
#'     with zoom/pan). Includes non-A residue overlay highlighting,
#'     filterable read lists, and Previous/Next navigation.}
#'   \item{\strong{Download}}{Configurable report generator producing a
#'     self-contained HTML file. Users select which sections to include
#'     (classification, abundance, residue frequency, rug density, poly(A)
#'     distribution, example signal plots) and configure plot settings.
#'     Supports per-transcript sub-reports for up to 3 selected
#'     transcripts. All plots include descriptive annotations.}
#'   \item{\strong{About}}{Package version, citation (Nat Commun 2025),
#'     links to GitHub, Wiki, pkgdown site, Zenodo DOI, laboratory and
#'     developer contact information.}
#' }
#'
#' The dashboard supports two usage modes:
#'
#' \strong{Multi-sample mode} (\code{config}): Loads a YAML configuration
#' file pointing to multiple samples with class and residue data. Enables
#' comparative analysis across samples and experimental groups. All
#' plotting functions use \code{sample_name} and \code{group} columns for
#' grouping.
#'
#' \strong{Single-sample mode}: Provide individual file paths. All
#' analysis tabs are available when \code{class_file} and
#' \code{residue_file} are supplied. The Signal Viewer tab requires
#' \code{summary_file} and \code{pod5_dir}. Default \code{sample_name}
#' and \code{group} columns are added automatically if missing.
#'
#' @param config Character string (optional). Path to a YAML configuration
#'   file defining multiple samples with file paths and experimental
#'   groups. When provided, single-sample arguments are ignored.
#'   See Details for the expected format.
#'
#' @param summary_file Character string (optional). Path to a Dorado
#'   summary file (tab-separated, with columns \code{read_id},
#'   \code{filename}, \code{poly_tail_start}, \code{poly_tail_end}).
#'   Required for the Signal Viewer tab in single-sample mode.
#'
#' @param pod5_dir Character string (optional). Path to the directory
#'   containing POD5 files from the sequencing run. Required for the
#'   Signal Viewer tab in single-sample mode.
#'
#' @param class_file Character string (optional). Path to the
#'   \code{read_classes} output file from a ninetails pipeline
#'   (\code{check_tails_dorado_DRS}, \code{check_tails_dorado_cDNA}, or
#'   \code{check_tails_guppy}). Enables the Classification, Residues,
#'   Poly(A) length, and Download tabs.
#'
#' @param residue_file Character string (optional). Path to a
#'   \code{nonadenosine_residues} output file from ninetails. Enables
#'   residue-specific plots (abundance, residue counts, rug density)
#'   and non-A overlay on signal plots.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}},
#'   such as \code{port}, \code{host}, or \code{launch.browser}.
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
#' Required fields per sample: \code{sample_name}, \code{group},
#' \code{class_path}, \code{residue_path}. The \code{dorado_summary} and
#' \code{pod5_dir} fields are optional and only needed for the Signal
#' Viewer tab.
#'
#' A template configuration file is included in the package:
#' \code{system.file("extdata", "config_template.yml", package = "ninetails")}
#'
#' The dashboard requires the following packages (listed in Suggests):
#' \code{shiny}, \code{plotly}, \code{htmltools}, \code{DT},
#' \code{base64enc}, \code{cowplot}. For multi-sample mode: \code{yaml}.
#' For the Signal Viewer: a Python environment with the \code{pod5}
#' package, accessible via \code{reticulate}.
#'
#' Static assets (\code{logo.png}, \code{favicon.ico},
#' \code{IIMCB_logo.png}) should be placed in \code{inst/app/www/}.
#'
#' For deployment to Shiny Server, a wrapper script is included at
#' \code{system.file("deployment", "app.R", package = "ninetails")}.
#' The wrapper has two configuration variables at the top:
#' \code{config_path} (path to YAML config) and \code{PYTHON_PATH}
#' (path to Python binary with \code{pod5} module, required for the
#' Signal Viewer tab). See \code{vignette("shiny_app")} for full
#' deployment instructions.
#'
#' @return Launches a Shiny application (does not return a value).
#'
#' @seealso
#' \code{\link{check_tails_dorado_DRS}} for the main analysis pipeline,
#' \code{\link{merge_nonA_tables}} for merging output tables,
#' \code{\link{annotate_with_biomart}} for adding gene symbols,
#' \code{\link{plot_class_counts}} and other plotting functions for
#' static equivalents of dashboard plots.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Multi-sample mode (recommended for comparative analysis)
#' ninetails::launch_signal_browser(config = "/path/to/config.yml")
#'
#' # Single-sample: all tabs (classification + residues + signal viewer)
#' ninetails::launch_signal_browser(
#'   summary_file = "/path/to/dorado_summary.txt",
#'   pod5_dir     = "/path/to/pod5/",
#'   class_file   = "/path/to/read_classes.txt",
#'   residue_file = "/path/to/nonadenosine_residues.txt"
#' )
#'
#' # Single-sample: signal viewer only (no analysis tabs)
#' ninetails::launch_signal_browser(
#'   summary_file = "/path/to/dorado_summary.txt",
#'   pod5_dir     = "/path/to/pod5/"
#' )
#'
#' # Custom port and host for remote access
#' ninetails::launch_signal_browser(
#'   config = "config.yml",
#'   port = 8080,
#'   host = "0.0.0.0"
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

  ####### Multi-sample mode (YAML config) ######
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
    ####### Single-sample mode ######

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

    # Add default sample_name and group columns if missing
    # (required by merge_nonA_tables and downstream plotting functions)
    if (!is.null(class_data)) {
      if (!"sample_name" %in% names(class_data)) {
        class_data$sample_name <- as.factor("sample_1")
      }
      if (!"group" %in% names(class_data)) {
        class_data$group <- as.factor("ungrouped")
      }
    }
    if (!is.null(residue_data)) {
      if (!"sample_name" %in% names(residue_data)) {
        residue_data$sample_name <- as.factor("sample_1")
      }
      if (!"group" %in% names(residue_data)) {
        residue_data$group <- as.factor("ungrouped")
      }
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


#' Launch the cDNA orientation validation browser
#'
#' Minimum-useful companion to \code{\link{launch_signal_browser}}
#' specifically for verifying cDNA orientation classifications against
#' the raw signal layout. Shows one signal at a time with the algorithm's
#' \code{tail_type} call, an independent size-ratio-inferred layout
#' (see \code{\link{infer_cdna_layout}}), and an agreement marker so
#' disagreements stand out visually.
#'
#' Intended as a validation tool while
#' \code{\link{check_tails_dorado_cDNA}} is still under construction;
#' will be superseded by a fuller cDNA dashboard once cDNA-native CNN
#' submodels are in place.
#'
#' @param dorado_summary Character string or data frame. Path to a
#'   Dorado summary file or a data frame with at least \code{read_id},
#'   \code{filename}, \code{poly_tail_start}, \code{poly_tail_end}
#'   columns. The Dorado >=1.4.0 \code{input_filename} alias is also
#'   accepted.
#' @param pod5_dir Character string. Path to the directory containing
#'   POD5 files for the run.
#' @param orientation_data Character string or data frame. Path to a TSV
#'   produced by \code{\link{detect_orientation_multiple}}, or that
#'   function's in-memory return value. Must contain \code{read_id} and
#'   \code{tail_type}; \code{reference} and \code{sequence} are used if
#'   present.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @return Launches a Shiny application (does not return a value).
#'
#' @seealso \code{\link{launch_signal_browser}},
#'   \code{\link{detect_orientation_multiple}},
#'   \code{\link{infer_cdna_layout}},
#'   \code{\link{cdna_layout_agreement_marker}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using files on disk
#' ninetails::launch_cdna_signal_browser(
#'   dorado_summary   = "/path/to/dorado_summary.txt",
#'   pod5_dir         = "/path/to/pod5/",
#'   orientation_data = "/path/to/sequence_with_tail_type.tsv"
#' )
#'
#' # Using an in-memory orientation tibble (e.g. straight out of the pipeline)
#' classified <- ninetails::detect_orientation_multiple(
#'   sequence_files = my_seq_files, num_cores = 4
#' )
#' ninetails::launch_cdna_signal_browser(
#'   dorado_summary   = "/path/to/dorado_summary.txt",
#'   pod5_dir         = "/path/to/pod5/",
#'   orientation_data = classified
#' )
#' }
launch_cdna_signal_browser <- function(dorado_summary,
                                       pod5_dir,
                                       orientation_data,
                                       ...) {

  for (pkg in c("shiny", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' required. Install with: install.packages('", pkg, "')",
           call. = FALSE)
    }
  }

  # #### pod5_dir ####
  if (missing(pod5_dir)) {
    stop("pod5_dir is required.", call. = FALSE)
  }
  assert_condition(is_string(pod5_dir),
                   "pod5_dir must be a character string")
  assert_dir_exists(pod5_dir, "pod5_dir")

  # #### dorado_summary ####
  if (missing(dorado_summary)) {
    stop("dorado_summary is required.", call. = FALSE)
  }
  if (is.character(dorado_summary) && length(dorado_summary) == 1) {
    assert_file_exists(dorado_summary, "dorado_summary")
    summary_df <- vroom::vroom(dorado_summary, show_col_types = FALSE)
  } else if (is.data.frame(dorado_summary)) {
    summary_df <- as.data.frame(dorado_summary)
  } else {
    stop("dorado_summary must be a file path or data frame.", call. = FALSE)
  }

  # Handle Dorado >=1.4.0 column rename: input_filename -> filename
  if (!"filename" %in% colnames(summary_df) &&
      "input_filename" %in% colnames(summary_df)) {
    summary_df$filename <- summary_df$input_filename
  }

  required_summary_cols <- c("read_id", "filename",
                             "poly_tail_start", "poly_tail_end")
  missing_summary_cols <- setdiff(required_summary_cols, colnames(summary_df))
  if (length(missing_summary_cols) > 0) {
    stop("dorado_summary is missing required columns: ",
         paste(missing_summary_cols, collapse = ", "), call. = FALSE)
  }

  # #### orientation_data ####
  if (missing(orientation_data)) {
    stop(
      "orientation_data is required. Pass a file path or the in-memory ",
      "output of detect_orientation_multiple().",
      call. = FALSE)
  }
  if (is.character(orientation_data) && length(orientation_data) == 1) {
    assert_file_exists(orientation_data, "orientation_data")
    orient_df <- vroom::vroom(orientation_data, show_col_types = FALSE)
  } else if (is.data.frame(orientation_data)) {
    orient_df <- as.data.frame(orientation_data)
  } else {
    stop("orientation_data must be a file path or data frame.", call. = FALSE)
  }

  required_orient_cols <- c("read_id", "tail_type")
  missing_orient_cols <- setdiff(required_orient_cols, colnames(orient_df))
  if (length(missing_orient_cols) > 0) {
    stop("orientation_data is missing required columns: ",
         paste(missing_orient_cols, collapse = ", "), call. = FALSE)
  }

  optional_cols <- intersect(c("reference", "sequence"), colnames(orient_df))
  keep_cols <- c("read_id", "tail_type", optional_cols)
  orient_subset <- orient_df[, keep_cols, drop = FALSE]

  # #### Join ####
  joined_df <- dplyr::left_join(summary_df, orient_subset, by = "read_id")

  unmatched <- sum(is.na(joined_df$tail_type))
  if (unmatched > 0) {
    warning(sprintf(
      paste0("%d / %d reads in dorado_summary had no matching tail_type ",
             "in orientation_data; they will be displayed as 'unidentified'."),
      unmatched, nrow(joined_df)), call. = FALSE)
    joined_df$tail_type[is.na(joined_df$tail_type)] <- "unidentified"
  }

  cat(sprintf(
    "[%s] cDNA browser loaded %s reads (%d polyA, %d polyT, %d unidentified)\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    format(nrow(joined_df), big.mark = ","),
    sum(joined_df$tail_type == "polyA",        na.rm = TRUE),
    sum(joined_df$tail_type == "polyT",        na.rm = TRUE),
    sum(joined_df$tail_type == "unidentified", na.rm = TRUE)))

  shiny::shinyOptions(
    ninetails.cdna_data     = joined_df,
    ninetails.cdna_pod5_dir = pod5_dir
  )

  app_dir <- system.file("app_cdna", package = "ninetails")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("App not found. Expected at: inst/app_cdna/", call. = FALSE)
  }

  shiny::runApp(app_dir, ...)
}
