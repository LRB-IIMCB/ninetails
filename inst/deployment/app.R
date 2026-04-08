################################################################################
# Ninetails Dashboard — Shiny Server deployment wrapper
#
# Place this file in the Shiny Server app directory alongside config.yml.
# Shiny Server will run this file directly. It loads data from the YAML
# config and delegates to the app bundled in the installed ninetails package.
#
# Directory structure:
#   /srv/shiny-server/ninetails/
#   ├── app.R           (this file)
#   ├── config.yml      (your YAML configuration)
#   └── www/            (static assets — symlink or copy from inst/app/www/)
#       ├── logo.png
#       ├── favicon.ico
#       └── IIMCB_logo.png
#
# Setup:
#   1. Install ninetails and all Suggests on the server
#   2. Copy this file + config.yml to the app directory
#   3. Edit PYTHON_PATH below (required for Signal Viewer tab)
#   4. Symlink or copy www/ assets:
#        ln -s $(Rscript -e "cat(system.file('app/www', package='ninetails'))") www
#   5. Restart Shiny Server or let it auto-detect
#
################################################################################


################################################################################
# DEPLOYMENT CONFIGURATION — edit these before launching
################################################################################

# Path to YAML config (relative to this app.R directory)
config_path <- file.path(getwd(), "config.yml")

# Python environment for signal extraction (Signal Viewer tab).
# Set this to the Python binary that has the 'pod5' module installed.
# If Signal Viewer is not needed, leave as NULL.
#
# Examples:
#   PYTHON_PATH <- "/usr/bin/python3"
#   PYTHON_PATH <- "/home/user/miniconda3/envs/r-reticulate/bin/python"
#   PYTHON_PATH <- "/home/user/.virtualenvs/ninetails/bin/python"
#
PYTHON_PATH <- NULL

################################################################################

library(ninetails)
library(shiny)
library(ggplot2)
library(dplyr)
library(vroom)
library(reticulate)

# Apply Python path if specified (must happen before any reticulate calls)
if (!is.null(PYTHON_PATH) && nzchar(PYTHON_PATH)) {
  Sys.setenv(RETICULATE_PYTHON = PYTHON_PATH)
  cat(paste0("[", Sys.time(), "] Python path set to: ", PYTHON_PATH, "\n"))
}

####### Validate config ######

if (!file.exists(config_path)) {
  stop("config.yml not found in app directory: ", getwd(),
       "\nPlace your YAML configuration file alongside this app.R.",
       call. = FALSE)
}

####### Load data (mirrors launch_signal_browser() logic)######

cfg <- yaml::read_yaml(config_path)
if (is.null(cfg$samples) || length(cfg$samples) == 0) {
  stop("No samples found in config.yml", call. = FALSE)
}

class_list   <- list()
residue_list <- list()
signal_config <- list()

cat(paste0("[", Sys.time(), "] Loading ", length(cfg$samples), " samples...\n"))

for (sid in names(cfg$samples)) {
  s <- cfg$samples[[sid]]
  sname <- s$sample_name %||% sid
  grp   <- s$group %||% "ungrouped"

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

  if (!is.null(s$dorado_summary) && !is.null(s$pod5_dir)) {
    signal_config[[sname]] <- list(
      dorado_summary = s$dorado_summary,
      pod5_dir = s$pod5_dir
    )
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
      class_data, residue_data, pass_only = FALSE
    )
  }, error = function(e) {
    warning("Could not create merged table: ", e$message, call. = FALSE)
  })
}

cat(paste0("[", Sys.time(), "] Loaded ",
           format(nrow(class_data), big.mark = ","), " reads from ",
           length(unique(class_data$sample_name)), " samples.\n"))

####### Pass data to app via shinyOptions ######

shiny::shinyOptions(
  ninetails.class_data    = class_data,
  ninetails.residue_data  = residue_data,
  ninetails.merged_data   = merged_data,
  ninetails.signal_config = signal_config,
  ninetails.summary_file  = "",
  ninetails.pod5_dir      = "",
  ninetails.residue_file  = ""
)

####### Source the app from the installed package ######

app_dir <- system.file("app", package = "ninetails")
if (!nzchar(app_dir) || !dir.exists(app_dir)) {
  stop("ninetails app not found. Is the package installed?", call. = FALSE)
}

# Source into a local environment, extract ui and server
app_env <- new.env(parent = globalenv())
source(file.path(app_dir, "app.R"), local = app_env)

# Shiny Server expects shinyApp() as the last expression
shiny::shinyApp(ui = app_env$ui, server = app_env$server)
