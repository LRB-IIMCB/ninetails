################################################################################
# Ninetails Dashboard — Universal Shiny Server deployment wrapper
#
# Auto-detects whether the YAML config describes Dorado DRS or Guppy legacy
# data based on the fields present in the first sample entry:
#   - pod5_dir or dorado_summary  →  Dorado DRS app
#   - workspace or nanopolish     →  Guppy legacy app
#
# Place this file as app.R in the Shiny Server app directory alongside
# config.yml and a www/ folder with static assets.
#
# Directory structure:
#   /srv/shiny-server/my_analysis/
#   ├── app.R           (this file)
#   ├── config.yml      (YAML configuration)
#   └── www/            (symlink or copy from inst/app/www/ or inst/app_guppy/www/)
#
# Setup:
#   1. Install ninetails and all Suggests on the server
#   2. Copy this file as app.R to the Shiny Server app directory
#   3. Place config.yml alongside it
#   4. Copy or symlink www/ assets
#   5. Edit PYTHON_PATH below if needed (Dorado Signal Viewer only)
#   6. Restart Shiny Server
#
################################################################################


################################################################################
# DEPLOYMENT CONFIGURATION — edit these before launching
################################################################################

# Path to YAML config (relative to this app.R directory)
config_path <- file.path(getwd(), "config.yml")

# Python environment for signal extraction (Dorado Signal Viewer only).
# Set to the Python binary with the 'pod5' module installed.
# Not needed for Guppy pipeline or if Signal Viewer is not used.
PYTHON_PATH <- NULL

# Basecall group for Guppy fast5 access (ignored for Dorado)
BASECALL_GROUP <- "Basecall_1D_000"

################################################################################

library(ninetails)
library(shiny)
library(ggplot2)
library(dplyr)
library(vroom)

# Apply Python path if specified
if (!is.null(PYTHON_PATH) && nzchar(PYTHON_PATH)) {
  Sys.setenv(RETICULATE_PYTHON = PYTHON_PATH)
  cat(paste0("[", Sys.time(), "] Python path set to: ", PYTHON_PATH, "\n"))
}

# ---- Validate config ----

if (!file.exists(config_path)) {
  stop("config.yml not found in app directory: ", getwd(),
       "\nPlace your YAML configuration file alongside this app.R.",
       call. = FALSE)
}

cfg <- yaml::read_yaml(config_path)
if (is.null(cfg$samples) || length(cfg$samples) == 0) {
  stop("No samples found in config.yml", call. = FALSE)
}

# ---- Auto-detect pipeline ----

first_sample <- cfg$samples[[1]]
is_dorado <- !is.null(first_sample$pod5_dir) || !is.null(first_sample$dorado_summary)
is_guppy  <- !is.null(first_sample$workspace) || !is.null(first_sample$nanopolish)

if (is_dorado) {
  pipeline <- "dorado"
  app_name <- "app"
  cat(paste0("[", Sys.time(), "] Detected Dorado DRS pipeline\n"))
} else if (is_guppy) {
  pipeline <- "guppy"
  app_name <- "app_guppy"
  cat(paste0("[", Sys.time(), "] Detected Guppy legacy pipeline\n"))
} else {
  # No signal fields — detect from class/residue data only
  pipeline <- "dorado"
  app_name <- "app"
  cat(paste0("[", Sys.time(), "] No signal config detected, defaulting to Dorado app\n"))
}

# ---- Load data ----

class_list   <- list()
residue_list <- list()
signal_config <- list()

cat(paste0("[", Sys.time(), "] Loading ", length(cfg$samples), " samples...\n"))

for (sid in names(cfg$samples)) {
  s <- cfg$samples[[sid]]
  sname <- if (!is.null(s$sample_name)) s$sample_name else sid
  grp   <- if (!is.null(s$group)) s$group else "ungrouped"

  if (!is.null(s$class_path) && file.exists(s$class_path)) {
    tryCatch({
      cd <- vroom::vroom(s$class_path, show_col_types = FALSE)
      cd$sample_name <- as.factor(sname)
      cd$group <- as.factor(grp)
      class_list[[sid]] <- cd
    }, error = function(e) {
      warning("class_path for '", sid, "': ", e$message, call. = FALSE)
    })
  }

  if (!is.null(s$residue_path) && file.exists(s$residue_path)) {
    tryCatch({
      rd <- vroom::vroom(s$residue_path, show_col_types = FALSE)
      rd$sample_name <- as.factor(sname)
      rd$group <- as.factor(grp)
      residue_list[[sid]] <- rd
    }, error = function(e) {
      warning("residue_path for '", sid, "': ", e$message, call. = FALSE)
    })
  }

  # Build signal config based on pipeline
  if (pipeline == "dorado") {
    if (!is.null(s$dorado_summary) && !is.null(s$pod5_dir)) {
      signal_config[[sname]] <- list(
        dorado_summary = s$dorado_summary,
        pod5_dir = s$pod5_dir
      )
    }
  } else {
    if (!is.null(s$nanopolish) && !is.null(s$sequencing_summary) &&
        !is.null(s$workspace)) {
      signal_config[[sname]] <- list(
        nanopolish_path = s$nanopolish,
        sequencing_summary_path = s$sequencing_summary,
        workspace = s$workspace
      )
    }
  }
}

if (length(class_list) == 0) {
  stop("No class data loaded. Check class_path entries in config.yml.",
       call. = FALSE)
}

class_data <- dplyr::bind_rows(class_list)
class_data$sample_name <- as.factor(class_data$sample_name)
class_data$group <- as.factor(class_data$group)

residue_data <- if (length(residue_list) > 0) {
  rd <- dplyr::bind_rows(residue_list)
  rd$sample_name <- as.factor(rd$sample_name)
  rd$group <- as.factor(rd$group)
  rd
} else { NULL }

merged_data <- NULL
if (!is.null(class_data) && !is.null(residue_data) && nrow(residue_data) > 0) {
  tryCatch({
    merged_data <- ninetails::merge_nonA_tables(
      class_data, residue_data, pass_only = FALSE)
  }, error = function(e) {
    warning("Could not create merged table: ", e$message, call. = FALSE)
  })
}

cat(paste0("[", Sys.time(), "] Loaded ",
           format(nrow(class_data), big.mark = ","), " reads from ",
           length(unique(class_data$sample_name)), " samples.\n"))

# ---- Pass data to app via shinyOptions ----

shiny::shinyOptions(
  ninetails.class_data    = class_data,
  ninetails.residue_data  = residue_data,
  ninetails.merged_data   = merged_data,
  ninetails.signal_config = signal_config,
  ninetails.summary_file  = "",
  ninetails.pod5_dir      = "",
  ninetails.residue_file  = "",
  ninetails.basecall_group = BASECALL_GROUP
)

# ---- Source the appropriate app ----

app_dir <- system.file(app_name, package = "ninetails")
if (!nzchar(app_dir) || !dir.exists(app_dir)) {
  stop("ninetails ", pipeline, " app not found at inst/", app_name,
       "/. Is the package installed?", call. = FALSE)
}

cat(paste0("[", Sys.time(), "] Launching ", pipeline, " dashboard from: ",
           app_dir, "\n"))

app_env <- new.env(parent = globalenv())
source(file.path(app_dir, "app.R"), local = app_env)

shiny::shinyApp(ui = app_env$ui, server = app_env$server)
