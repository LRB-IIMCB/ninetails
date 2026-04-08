################################################################################
# Ninetails Analysis Dashboard
#
# Interactive Shiny application for exploring poly(A) tail composition
# analysis results produced by the ninetails package. Provides a unified
# interface for inspecting read classification, non-adenosine residue
# composition, poly(A) length distributions, and raw nanopore signals.
#
# The dashboard is designed for post-analysis exploration of ninetails
# output. It accepts either a YAML configuration file (multi-sample mode)
# or individual file paths (single-sample mode) and adapts its active
# tabs based on which data is available. All plots are generated using
# the same plotting functions available in the ninetails package, with
# additional interactive features (filters, condition selectors,
# togglable descriptions, per-sample rug density plots with subsampling,
# and a configurable HTML report generator).
#
# Location: inst/app/app.R
# Launch:   ninetails::launch_signal_browser()
# Docs:     vignette("shiny_app", package = "ninetails")
#
# Tabs:
#   1. Classification — value boxes (samples, transcripts, reads,
#      blank, decorated) + plot_class_counts + plot_nonA_abundance
#      in responsive flex layout with togglable descriptions
#   2. Residues — plot_residue_counts + per-sample rug density plots
#      (C/G/U, max 1000 points each) + summary DT table
#   3. Poly(A) length — plot_tail_distribution with condition filter,
#      palette picker, central tendency, reset button, and summary
#      statistics table (n, mean, median, SD, SEM)
#   4. Signal Viewer — sub-tabs: Static Viewer (ggplot2) / Dynamic
#      Explorer (plotly) with non-A residue overlay, read navigation,
#      and filterable read list
#   5. Download — configurable report generator with section checkboxes,
#      plot settings, per-transcript sub-reports (up to 3), and signal
#      plots (5 random reads per category per sample)
#   6. About — package info, citation, IIMCB logo, links, credits
#
# Static assets (inst/app/www/):
#   - logo.png      (copy from man/figures/logo.png)
#   - favicon.ico   (copy from pkgdown/favicon/)
#   - IIMCB_logo.png (IIMCB institute logo)
#
# Dependencies (Suggests):
#   shiny, plotly, htmltools, DT, base64enc, cowplot, yaml, reticulate
#
################################################################################

library(shiny)
library(ggplot2)
library(dplyr)
library(vroom)
library(reticulate)

# Disable scientific notation on plot axes (e.g. 2,000 instead of 2e+03)
options(scipen = 999)


################################################################################
# DATA FROM LAUNCHER
################################################################################

class_data    <- shiny::getShinyOption("ninetails.class_data")
residue_data  <- shiny::getShinyOption("ninetails.residue_data")
merged_data   <- shiny::getShinyOption("ninetails.merged_data")
signal_config <- shiny::getShinyOption("ninetails.signal_config", list())

default_summary <- shiny::getShinyOption("ninetails.summary_file", default = "")
default_pod5    <- shiny::getShinyOption("ninetails.pod5_dir", default = "")
default_residue <- shiny::getShinyOption("ninetails.residue_file", default = "")

has_class   <- !is.null(class_data) && nrow(class_data) > 0
has_residue <- !is.null(residue_data) && nrow(residue_data) > 0
has_merged  <- !is.null(merged_data) && nrow(merged_data) > 0
has_signal  <- length(signal_config) > 0 ||
  nzchar(default_summary) || nzchar(default_pod5)

# Startup diagnostics
cat("--- Ninetails Dashboard Startup ---\n")
if (has_class)   cat("  class_data:   ", nrow(class_data), "rows\n")
if (has_residue) cat("  residue_data: ", nrow(residue_data), "rows\n")
if (has_merged)  cat("  merged_data:  ", nrow(merged_data), "rows\n")
cat("  signal_config:", length(signal_config), "sample(s)\n")

# Dynamic grouping variables
available_groups <- character(0)
if (has_class) {
  candidates <- c("sample_name", "group")
  available_groups <- candidates[candidates %in% names(class_data)]
  if (length(available_groups) == 0) available_groups <- "contig"
}

# Transcript label column: prefer symbol > contig
transcript_col <- "contig"
if (has_class && "symbol" %in% names(class_data)) transcript_col <- "symbol"

# Pre-compute summary stats for value boxes
n_samples_val     <- if (has_class && "sample_name" %in% names(class_data))
  length(unique(class_data$sample_name)) else 1
n_transcripts_val <- if (has_class) length(unique(class_data[[transcript_col]])) else 0
n_reads_val       <- if (has_class) nrow(class_data) else 0
n_blank_val       <- if (has_class && "class" %in% names(class_data))
  sum(class_data$class == "blank", na.rm = TRUE) else 0
n_decorated_val   <- if (has_class && "class" %in% names(class_data))
  sum(class_data$class == "decorated", na.rm = TRUE) else 0


################################################################################
# COLOR PALETTES
################################################################################

palettes <- list(
  kanto  = c("#D62728FF", "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#9467BDFF",
             "#8C564BFF", "#79AF97FF", "#EFC000FF", "#E377C2FF", "#7F7F7FFF",
             "#bdff17", "#b5004b", "#17BECFFF", "#87bdf5", "#3a424f"),
  johto  = c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#EFC000FF",
             "#925E9FFF", "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF",
             "#8F7700FF", "#EE4C97FF", "#6F99ADFF", "#B09C85FF", "#CD534CFF"),
  hoenn  = c("#5A9599FF", "#FF6F00FF", "#49a33b", "#8A4198FF", "#84D7E1FF",
             "#FF6348FF", "#3C5488FF", "#FF95A8FF", "#3D3B25FF", "#ADE2D0FF",
             "#CD534CFF", "#FFDC91FF", "#50a675", "#C71000FF", "#008EA0FF"),
  sinnoh = c("#ff1726", "#196abd", "#1fed29", "#ff6200", "#159485",
             "#00aad4", "#db63e6", "#374e55ff", "#87bdf5", "#bdff17",
             "#b5004b", "#fdff70", "#00dbd4", "#ffb300", "#6b43b0"),
  hisui  = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF",
             "#42B540FF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF",
             "#FDAE61", "#79AF97FF", "#6A6599FF", "#EE4C97FF", "#80796BFF"),
  unova  = c("#D7191C", "#2C7BB6", "#FDAE61", "#ABD9E9", "#3a424f",
             "#50a675", "#b0bdd4", "#fdff70", "#808080", "#00aad4",
             "#1a1a1a", "#cccccc", "#ff6600", "#216778", "#b5004b"),
  kalos  = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
             "#8491B4FF", "#374E55FF", "#00A1D5FF", "#79AF97FF", "#c7081b",
             "#0559ff", "#3d0101", "#4a0869", "#b03b0c", "#0f05a1"),
  alola  = c("#ff1726", "#00aad4", "#db63e6", "#ffb300", "#1fed29",
             "#87bdf5", "#196abd", "#b5004b", "#FDAE61", "#3a424f",
             "#216778", "#50a675", "#808080", "#00dbd4", "#6b43b0")
)


################################################################################
# HELPER: custom value box HTML
################################################################################

.value_box_ui <- function(value, label, color = "#2C7BB6",
                          icon_url = "logo.png") {
  shiny::div(
    class = "ninetails-vbox",
    style = paste0("border-left: 5px solid ", color, ";"),
    # Watermark icon
    shiny::tags$img(
      src = icon_url, class = "vbox-watermark",
      style = "opacity: 0.06;"
    ),
    # Content
    shiny::div(class = "vbox-content",
               shiny::div(class = "vbox-value", style = paste0("color: ", color, ";"),
                          value),
               shiny::div(class = "vbox-label", style = paste0("color: ", color, ";"),
                          label)
    ) # div content
  ) # div vbox
}


################################################################################
# UI
################################################################################

ui <- shiny::fluidPage(

  shiny::tags$head(
    shiny::tags$link(rel = "stylesheet",
                     href = "https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;600;700&display=swap"),
    shiny::tags$link(rel = "icon", type = "image/x-icon", href = "favicon.ico"),
    shiny::tags$style(shiny::HTML("
      body { font-family: 'Open Sans', sans-serif; background-color: #f9fafb; }

      /* ---- Header ---- */
      .header {
        background-color: #dce2e8; color: #333; padding: 16px 20px;
        margin: -15px -15px 20px -15px;
        display: flex; align-items: center; gap: 14px;
      }
      .header h2 { margin: 0; font-weight: bold; font-size: 22px; }
      .header p { margin: 4px 0 0 0; color: #555; font-size: 14px; }

      /* ---- Custom value boxes ---- */
      .vbox-row {
        display: flex; gap: 12px; margin-bottom: 15px; flex-wrap: wrap;
      }
      .ninetails-vbox {
        background: #fff; border-radius: 6px; padding: 18px 20px;
        box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        position: relative; overflow: hidden; min-height: 80px;
        flex: 1 1 0; min-width: 140px;
      }
      .vbox-watermark {
        position: absolute; right: 10px; top: 50%;
        transform: translateY(-50%);
        height: 60px; width: 60px; pointer-events: none;
      }
      .vbox-content { position: relative; z-index: 1; }
      .vbox-value { font-size: 26px; font-weight: 700; line-height: 1.1; }
      .vbox-label { font-size: 14px; font-weight: 700; margin-top: 3px; }

      /* ---- Cards ---- */
      .card {
        background: #fff; border-radius: 6px; padding: 15px;
        margin-bottom: 15px; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
      }
      .card h4 {
        margin-top: 0; color: #333; font-weight: 600;
        border-bottom: 2px solid #2C7BB6; padding-bottom: 8px;
      }
      .filter-section {
        background: #f5f7f9; padding: 15px; border-radius: 6px; margin-bottom: 15px;
      }
      .placeholder-msg { text-align: center; padding: 60px 20px; color: #888; }
      .placeholder-msg h4 { color: #aaa; }
      .read-info {
        background: #f0f5fa; padding: 12px 15px; border-radius: 6px;
        border-left: 4px solid #2C7BB6;
      }
      .read-info p { margin: 5px 0; }
      .classification-info {
        background: #e4eef7; padding: 8px 12px; border-radius: 4px;
        margin-top: 10px; font-size: 12px; color: #1a4971;
      }
      .residue-info {
        background: #e8f5e9; padding: 8px 12px; border-radius: 4px;
        margin-top: 8px; font-size: 12px; color: #2e7d32;
        border-left: 3px solid #50a675;
      }
      .paths-section {
        background: #f0f5fa; padding: 12px 15px; border-radius: 6px;
        margin-bottom: 12px; border-left: 3px solid #2C7BB6;
      }
      .paths-section .form-group { margin-bottom: 8px; }
      .btn-primary { background-color: #2C7BB6; border-color: #2569a0; }
      .btn-primary:hover { background-color: #2569a0; }
      .selectize-input, .selectize-dropdown { font-size: 12px; }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        color: #2C7BB6; border-bottom: 2px solid #2C7BB6; font-weight: 600;
      }
      .nav-tabs > li > a { color: #555; }

      /* ---- Flex plot row (side-by-side when viewport allows) ---- */
      .plots-flex {
        display: flex; flex-wrap: wrap; gap: 15px;
      }
      .plots-flex > .card {
        flex: 1 1 400px; max-width: 100%;
      }
      @media (min-width: 1200px) {
        .plots-flex > .card { max-width: calc(50% - 8px); }
      }

      /* ---- About tab ---- */
      .about-section { max-width: 800px; margin: 0 auto; padding: 20px 0; }
      .about-section h3 { color: #2C7BB6; border-bottom: 2px solid #dce2e8; padding-bottom: 6px; }
      .about-section a { color: #2C7BB6; }
      .about-logo { display: block; margin: 0 auto 20px auto; }

      /* ---- Plot descriptions ---- */
      .plot-desc {
        background: #f5f7f9; padding: 10px 14px; border-radius: 4px;
        font-size: 13px; color: #555; line-height: 1.5;
        border-left: 3px solid #2C7BB6; margin-top: 8px;
      }
      .plot-desc b { color: #333; }
    ")) # tags$style
  ), # tags$head

  # ---- Header bar ----
  shiny::div(class = "header",
             htmltools::img(src = "logo.png", height = 45, width = 45, alt = "ninetails"),
             shiny::div(
               shiny::h2(htmltools::strong(paste0(
                 "Ninetails analysis dashboard ",
                 as.character(utils::packageVersion("ninetails"))
               ))),
               shiny::p("Poly(A) tail composition analysis and signal visualization")
             ) # div title
  ), # div header

  # ---- Main tabs ----
  shiny::tabsetPanel(
    id = "main_tabs", type = "tabs",


    ############################################################################
    # TAB 1: Classification (with value boxes on top)
    ############################################################################
    shiny::tabPanel("Classification", shiny::br(),
                    if (has_class) {
                      shiny::sidebarLayout(
                        shiny::sidebarPanel(width = 2,
                                            shiny::div(class = "filter-section",
                                                       shiny::selectInput("class_grouping", "Grouping variable",
                                                                          choices = available_groups, selected = available_groups[1]),
                                                       shiny::selectizeInput("class_contig_filter", "Filter by transcript",
                                                                             choices = NULL,
                                                                             options = list(placeholder = "All transcripts", maxOptions = 200)),
                                                       shiny::checkboxInput("class_frequency", "Show as frequency", FALSE),
                                                       shiny::selectInput("class_plot_type", "Classification view",
                                                                          choices = c("Summary (N)" = "N", "Detailed (R)" = "R",
                                                                                      "Decorated (A)" = "A"),
                                                                          selected = "N"),
                                                       shiny::hr(),
                                                       shiny::checkboxInput("show_class_desc", "Show descriptions", FALSE)
                                            ) # div filter-section
                        ), # sidebarPanel
                        shiny::mainPanel(width = 10,
                                         # Value boxes row (inside mainPanel, above plots)
                                         shiny::div(class = "vbox-row",
                                                    .value_box_ui(
                                                      format(n_samples_val, big.mark = ","),
                                                      "Samples analyzed", "#2C7BB6"),
                                                    .value_box_ui(
                                                      format(n_transcripts_val, big.mark = ","),
                                                      "Transcripts found", "#2C7BB6"),
                                                    .value_box_ui(
                                                      format(n_reads_val, big.mark = ","),
                                                      "Total reads", "#555555"),
                                                    .value_box_ui(
                                                      format(n_blank_val, big.mark = ","),
                                                      "Blank reads", "#089bcc"),
                                                    .value_box_ui(
                                                      format(n_decorated_val, big.mark = ","),
                                                      "Decorated reads", "#ff6600")
                                         ), # div vbox-row
                                         shiny::div(class = "plots-flex",
                                                    shiny::div(class = "card",
                                                               shiny::h4("Read Classification"),
                                                               shiny::plotOutput("class_plot", height = "450px"),
                                                               shiny::conditionalPanel("input.show_class_desc",
                                                                                       shiny::div(class = "plot-desc",
                                                                                                  shiny::HTML(paste0(
                                                                                                    "Read classification summary. Reads are categorized as ",
                                                                                                    "<b>decorated</b> (containing non-adenosine residues), ",
                                                                                                    "<b>blank</b> (no modifications detected), or ",
                                                                                                    "<b>unclassified</b> (insufficient data quality). ",
                                                                                                    "Classification codes: ",
                                                                                                    "<b>YAY</b> \u2014 move transition present, non-A detected; ",
                                                                                                    "<b>MPU</b> \u2014 move transition present, only A detected; ",
                                                                                                    "<b>MAU</b> \u2014 no move transition, unmodified tail; ",
                                                                                                    "<b>IRL</b> \u2014 insufficient read length (<10 nt); ",
                                                                                                    "<b>BAC</b> \u2014 bad coordinates; ",
                                                                                                    "<b>UNM</b> \u2014 unmapped read."
                                                                                                  ))
                                                                                       ) # div plot-desc
                                                               ) # conditionalPanel
                                                    ), # div card
                                                    if (has_residue) {
                                                      shiny::div(class = "card",
                                                                 shiny::h4("Non-A Abundance (reads with 1, 2, 3+ non-As)"),
                                                                 shiny::plotOutput("nonA_abundance_plot", height = "400px"),
                                                                 shiny::conditionalPanel("input.show_class_desc",
                                                                                         shiny::div(class = "plot-desc",
                                                                                                    shiny::HTML(paste0(
                                                                                                      "Frequency of reads containing one, two, or three or ",
                                                                                                      "more separate non-adenosine residues per read. ",
                                                                                                      "Computed relative to the total number of decorated ",
                                                                                                      "reads in each sample."
                                                                                                    ))
                                                                                         ) # div plot-desc
                                                                 ) # conditionalPanel
                                                      ) # div card
                                                    }
                                         ) # div plots-flex
                        ) # mainPanel
                      ) # sidebarLayout
                    } else {
                      shiny::div(class = "placeholder-msg",
                                 shiny::h4("No classification data loaded"),
                                 shiny::p("Provide class_file or a YAML config."))
                    }
    ), # tabPanel Classification


    ############################################################################
    # TAB 2: Residues (with rug density plots)
    ############################################################################
    shiny::tabPanel("Residues", shiny::br(),
                    if (has_residue) {
                      shiny::sidebarLayout(
                        shiny::sidebarPanel(width = 2,
                                            shiny::div(class = "filter-section",
                                                       shiny::selectInput("residue_grouping", "Grouping variable",
                                                                          choices = available_groups, selected = available_groups[1]),
                                                       shiny::selectizeInput("residue_contig_filter", "Filter by transcript",
                                                                             choices = NULL,
                                                                             options = list(placeholder = "All transcripts", maxOptions = 200)),
                                                       shiny::checkboxInput("residue_frequency", "Show as frequency", FALSE),
                                                       shiny::checkboxInput("residue_by_read", "Count by read", FALSE),
                                                       shiny::sliderInput("rug_max_length", "Max tail length (rug plot)",
                                                                          min = 10, max = 500, value = 200, step = 10),
                                                       shiny::hr(),
                                                       shiny::checkboxInput("show_residue_desc", "Show descriptions", FALSE)
                                            ) # div filter-section
                        ), # sidebarPanel
                        shiny::mainPanel(width = 10,
                                         shiny::div(class = "card",
                                                    shiny::h4("Residue Counts"),
                                                    shiny::plotOutput("residue_plot", height = "400px"),
                                                    shiny::conditionalPanel("input.show_residue_desc",
                                                                            shiny::div(class = "plot-desc",
                                                                                       shiny::HTML(paste0(
                                                                                         "Distribution of non-adenosine residue types ",
                                                                                         "(cytidine, guanosine, uridine) detected across samples. ",
                                                                                         "Toggle <b>Count by read</b> to count each read once per ",
                                                                                         "residue type, or count total individual residue occurrences."
                                                                                       ))
                                                                            ) # div plot-desc
                                                    ) # conditionalPanel
                                         ), # div card
                                         shiny::div(class = "card",
                                                    shiny::h4("Non-A Position Distribution (max 1000 points per residue per sample)"),
                                                    shiny::uiOutput("rug_plots_ui"),
                                                    shiny::conditionalPanel("input.show_residue_desc",
                                                                            shiny::div(class = "plot-desc",
                                                                                       shiny::HTML(paste0(
                                                                                         "Positional distribution of non-adenosine residues along ",
                                                                                         "poly(A) tails. Each point represents the estimated position ",
                                                                                         "(from the 3\u2019 end) of a detected modification. Marginal ",
                                                                                         "density curves show the overall distribution shape. Points ",
                                                                                         "are subsampled to max 1,000 per residue type per sample."
                                                                                       ))
                                                                            ) # div plot-desc
                                                    ) # conditionalPanel
                                         ), # div card rug
                                         if (has_merged) {
                                           shiny::div(class = "card",
                                                      shiny::h4("Summary Table"),
                                                      DT::DTOutput("residue_summary_table"))
                                         }
                        ) # mainPanel
                      ) # sidebarLayout
                    } else {
                      shiny::div(class = "placeholder-msg",
                                 shiny::h4("No residue data loaded"),
                                 shiny::p("Provide residue_file or a YAML config."))
                    }
    ), # tabPanel Residues


    ############################################################################
    # TAB 3: Poly(A) Distribution
    ############################################################################
    shiny::tabPanel("Poly(A) length", shiny::br(),
                    if (has_class) {
                      shiny::sidebarLayout(
                        shiny::sidebarPanel(width = 2,
                                            shiny::div(class = "filter-section",
                                                       shiny::selectInput("polya_grouping", "Grouping variable",
                                                                          choices = available_groups, selected = available_groups[1]),
                                                       shiny::selectizeInput("polya_contig_filter", "Filter by transcript",
                                                                             choices = NULL,
                                                                             options = list(placeholder = "All transcripts", maxOptions = 200)),
                                                       shiny::uiOutput("polya_condition_filter_ui"),
                                                       shiny::selectInput("polya_center", "Central tendency",
                                                                          choices = c("mean", "median", "mode", "none"), selected = "none"),
                                                       shiny::sliderInput("polya_max_length", "Max tail length",
                                                                          min = 0, max = 500, value = 200),
                                                       shiny::checkboxInput("polya_ndensity", "Normalized density", TRUE),
                                                       shiny::selectInput("polya_palette", "Color palette",
                                                                          choices = names(palettes), selected = "unova"),
                                                       shiny::hr(),
                                                       shiny::actionButton("polya_reset", "Reset all filters",
                                                                           class = "btn-default btn-sm", icon = shiny::icon("refresh"),
                                                                           style = "width: 100%;"),
                                                       shiny::hr(),
                                                       shiny::checkboxInput("show_polya_desc", "Show descriptions", FALSE)
                                            ) # div filter-section
                        ), # sidebarPanel
                        shiny::mainPanel(width = 10,
                                         shiny::div(class = "card",
                                                    shiny::h4("Poly(A) Tail Length Distribution"),
                                                    shiny::plotOutput("polya_dist_plot", height = "500px"),
                                                    shiny::conditionalPanel("input.show_polya_desc",
                                                                            shiny::div(class = "plot-desc",
                                                                                       shiny::HTML(paste0(
                                                                                         "Density distribution of poly(A) tail lengths across ",
                                                                                         "samples or experimental conditions. Use the condition ",
                                                                                         "filter to compare selected groups. Central tendency ",
                                                                                         "lines (mean, median, or mode) can be overlaid to ",
                                                                                         "highlight distribution centers."
                                                                                       ))
                                                                            ) # div plot-desc
                                                    ) # conditionalPanel
                                         ), # div card
                                         shiny::div(class = "card",
                                                    shiny::h4("Poly(A) Length Summary"),
                                                    DT::DTOutput("polya_summary_table"))
                        ) # mainPanel
                      ) # sidebarLayout
                    } else {
                      shiny::div(class = "placeholder-msg",
                                 shiny::h4("No data loaded"),
                                 shiny::p("Provide class_file or a YAML config."))
                    }
    ), # tabPanel Poly(A)


    ############################################################################
    # TAB 4: Signal Viewer
    ############################################################################
    shiny::tabPanel("Signal Viewer", shiny::br(),
                    shiny::fluidRow(

                      # ---- Left panel ----
                      shiny::column(3,
                                    shiny::div(class = "card",
                                               shiny::h4("Data"),
                                               if (length(signal_config) > 0 && names(signal_config)[1] != "single") {
                                                 shiny::div(class = "paths-section",
                                                            shiny::selectInput("signal_sample", "Select sample",
                                                                               choices = names(signal_config),
                                                                               selected = names(signal_config)[1]))
                                               } else if (length(signal_config) == 0) {
                                                 shiny::div(class = "paths-section",
                                                            shiny::textInput("signal_summary_file", "Summary file",
                                                                             value = "", placeholder = "/path/to/dorado_summary.txt"),
                                                            shiny::textInput("signal_pod5_dir", "POD5 directory",
                                                                             value = "", placeholder = "/path/to/pod5/"),
                                                            shiny::textInput("signal_residue_file",
                                                                             "Non-A residues file (optional)",
                                                                             value = "", placeholder = "/path/to/nonadenosine_residues.txt"),
                                                            shiny::actionButton("load_signal_data", "Load data",
                                                                                class = "btn-primary btn-sm",
                                                                                icon = shiny::icon("folder-open")))
                                               },
                                               shiny::uiOutput("signal_data_status")
                                    ), # div card data

                                    shiny::div(class = "card",
                                               shiny::h4("Filters"),
                                               shiny::div(class = "filter-section",
                                                          shiny::sliderInput("polya_length_range",
                                                                             "Poly(A) tail length range",
                                                                             min = 0, max = 500, value = c(10, 500), step = 1),
                                                          shiny::selectInput("comments_filter", "Decoration status",
                                                                             choices = c("All", "YAY", "MPU", "MAU", "UNM", "IRL", "BAC"),
                                                                             selected = "All"),
                                                          shiny::selectInput("sig_residue_filter", "Non-A residue type",
                                                                             choices = c("All", "C", "G", "U"), selected = "All"),
                                                          shiny::uiOutput("genome_filter_ui"),
                                                          shiny::uiOutput("mapq_filter_ui")
                                               ), # div filter-section
                                               shiny::hr(),
                                               shiny::uiOutput("read_selection_ui"),
                                               shiny::uiOutput("filter_summary")
                                    ) # div card filters
                      ), # column 3

                      # ---- Right panel with sub-tabs ----
                      shiny::column(9,
                                    shiny::tabsetPanel(id = "signal_sub_tabs", type = "tabs",

                                                       shiny::tabPanel("Static Viewer", shiny::br(),
                                                                       shiny::conditionalPanel("output.signal_loaded",
                                                                                               shiny::div(class = "card",
                                                                                                          shiny::h4("Selected Read Information"),
                                                                                                          shiny::uiOutput("read_info")),
                                                                                               shiny::div(class = "card",
                                                                                                          shiny::h4("Entire Signal"),
                                                                                                          shiny::plotOutput("full_signal_plot", height = "350px")),
                                                                                               shiny::div(class = "card",
                                                                                                          shiny::h4("Poly(A) Region (+/- 250 positions)"),
                                                                                                          shiny::plotOutput("polya_zoom_plot", height = "350px"))
                                                                       ), # conditionalPanel
                                                                       shiny::conditionalPanel("!output.signal_loaded",
                                                                                               shiny::div(class = "card",
                                                                                                          shiny::div(class = "placeholder-msg",
                                                                                                                     shiny::h4("No signal loaded"),
                                                                                                                     shiny::p("Load data and select a read.")))
                                                                       ) # conditionalPanel
                                                       ), # tabPanel Static Viewer

                                                       shiny::tabPanel("Dynamic Explorer", shiny::br(),
                                                                       shiny::conditionalPanel("output.signal_loaded",
                                                                                               shiny::div(class = "card",
                                                                                                          shiny::h4("Selected Read Information"),
                                                                                                          shiny::uiOutput("read_info_explorer")),
                                                                                               shiny::div(class = "card",
                                                                                                          shiny::h4("Interactive Signal Explorer"),
                                                                                                          shiny::p(style = "color: #666; font-size: 13px;",
                                                                                                                   "Click and drag to zoom. Double-click to reset."),
                                                                                                          plotly::plotlyOutput("interactive_signal_plot", height = "500px"))
                                                                       ), # conditionalPanel
                                                                       shiny::conditionalPanel("!output.signal_loaded",
                                                                                               shiny::div(class = "card",
                                                                                                          shiny::div(class = "placeholder-msg",
                                                                                                                     shiny::h4("No signal loaded"),
                                                                                                                     shiny::p("Load data and select a read.")))
                                                                       ) # conditionalPanel
                                                       ) # tabPanel Dynamic Explorer
                                    ) # tabsetPanel signal_sub_tabs
                      ) # column 9
                    ) # fluidRow
    ), # tabPanel Signal Viewer


    ############################################################################
    # TAB 5: Download (penultimate)
    ############################################################################
    shiny::tabPanel("Download", shiny::br(),
                    if (has_class) {
                      shiny::div(style = "max-width: 800px; margin: 0 auto;",
                                 shiny::div(class = "card",
                                            shiny::h4("Report Configuration"),
                                            shiny::p(style = "color: #666; font-size: 13px; margin-bottom: 15px;",
                                                     "Select which sections to include in the report.",
                                                     "Plots are embedded as images in a self-contained HTML file."),

                                            shiny::tags$h5(style = "color: #2C7BB6; font-weight: 600;",
                                                           "Global sections"),
                                            shiny::checkboxInput("rpt_classification", "Classification plot",
                                                                 value = TRUE),
                                            if (has_residue) shiny::checkboxInput("rpt_abundance",
                                                                                  "Non-A abundance (1, 2, 3+ per read)", value = TRUE),
                                            if (has_residue) shiny::checkboxInput("rpt_residue",
                                                                                  "Residue frequency", value = TRUE),
                                            if (has_residue) shiny::checkboxInput("rpt_rug",
                                                                                  "Rug density plots (per sample, with downsampling)", value = TRUE),
                                            shiny::checkboxInput("rpt_polya",
                                                                 "Poly(A) length distribution", value = TRUE),
                                            if (has_signal) shiny::checkboxInput("rpt_signals",
                                                                                 "Example signal plots (5 per category per sample)", value = FALSE),

                                            shiny::hr(),
                                            shiny::tags$h5(style = "color: #2C7BB6; font-weight: 600;",
                                                           "Plot settings"),
                                            shiny::selectInput("rpt_center", "Central tendency (poly(A) plot)",
                                                               choices = c("mean", "median", "mode", "none"),
                                                               selected = "none", width = "250px"),
                                            shiny::sliderInput("rpt_max_length", "Max poly(A) length",
                                                               min = 0, max = 500, value = 200, width = "250px"),
                                            shiny::selectInput("rpt_palette", "Color palette (poly(A) plot)",
                                                               choices = names(palettes), selected = "unova", width = "250px"),

                                            shiny::hr(),
                                            shiny::tags$h5(style = "color: #2C7BB6; font-weight: 600;",
                                                           "Per-transcript sections (optional)"),
                                            shiny::p(style = "color: #666; font-size: 12px;",
                                                     "Select up to 3 transcripts for detailed sub-reports.",
                                                     "Each will include classification, residue, and poly(A) plots."),
                                            shiny::selectizeInput("rpt_transcripts",
                                                                  "Select transcripts",
                                                                  choices = NULL, multiple = TRUE,
                                                                  options = list(placeholder = "None selected",
                                                                                 maxItems = 3, maxOptions = 200)),

                                            shiny::hr(),
                                            shiny::downloadButton("download_report", "Download report (.html)",
                                                                  class = "btn-primary",
                                                                  style = "width: 100%; margin-top: 10px;")
                                 ) # div card
                      ) # div wrapper
                    } else {
                      shiny::div(class = "placeholder-msg",
                                 shiny::h4("No data loaded"),
                                 shiny::p("Load data to enable report generation."))
                    }
    ), # tabPanel Download


    ############################################################################
    # TAB 6: About (last)
    ############################################################################
    shiny::tabPanel("About", shiny::br(),
                    shiny::div(class = "about-section",
                               htmltools::img(src = "logo.png", height = 120, width = 120,
                                              class = "about-logo", alt = "ninetails logo"),
                               shiny::h3("Ninetails"),
                               shiny::p(htmltools::strong(paste0(
                                 "Version ", as.character(utils::packageVersion("ninetails"))
                               ))),
                               shiny::p(
                                 "An R package for finding non-adenosine residues in poly(A) tails",
                                 "of Oxford Nanopore sequencing reads using convolutional neural",
                                 "networks based on raw current signal."
                               ),

                               shiny::h3("Citation"),
                               shiny::p(
                                 htmltools::HTML(paste0(
                                   "Gumi\u0144ska, N., Matylla-Kuli\u0144ska, K., Krawczyk, P.S. ",
                                   "<i>et al.</i> Direct profiling of non-adenosines in poly(A) tails ",
                                   "of endogenous and therapeutic mRNAs with Ninetails. ",
                                   "<i>Nat Commun</i> <b>16</b>, 2664 (2025). ",
                                   "<a href='https://doi.org/10.1038/s41467-025-57787-6' target='_blank'>",
                                   "https://doi.org/10.1038/s41467-025-57787-6</a>"
                                 ))
                               ),

                               shiny::h3("Links"),
                               shiny::tags$ul(
                                 shiny::tags$li(shiny::tags$a(
                                   href = "https://github.com/LRB-IIMCB/ninetails",
                                   target = "_blank", "GitHub repository")),
                                 shiny::tags$li(shiny::tags$a(
                                   href = "https://github.com/LRB-IIMCB/ninetails/wiki",
                                   target = "_blank", "Documentation (Wiki)")),
                                 shiny::tags$li(shiny::tags$a(
                                   href = "https://LRB-IIMCB.github.io/ninetails/",
                                   target = "_blank", "Package website (pkgdown)")),
                                 shiny::tags$li(shiny::tags$a(
                                   href = "https://doi.org/10.5281/zenodo.13309819",
                                   target = "_blank", "Zenodo DOI"))
                               ), # tags$ul

                               shiny::h3("Laboratory"),
                               htmltools::img(src = "IIMCB_logo.png", height = 90,
                                              class = "about-logo", alt = "IIMCB logo"),
                               shiny::p(htmltools::HTML(paste0(
                                 "<a href='https://www.iimcb.gov.pl/en/research/",
                                 "41-laboratory-of-rna-biology-era-chairs-group' target='_blank'>",
                                 "Laboratory of RNA Biology</a> (ERA Chairs Group), ",
                                 "International Institute of Molecular and Cell Biology in Warsaw (IIMCB)"
                               ))),

                               shiny::h3("Developer"),
                               shiny::p(htmltools::HTML(paste0(
                                 "Natalia Gumi\u0144ska ",
                                 "(<a href='mailto:nguminska@iimcb.gov.pl'>nguminska@iimcb.gov.pl</a>)"
                               ))),

                               shiny::hr(),
                               shiny::p(style = "color: #999; font-size: 12px; text-align: center;",
                                        paste0("Dashboard built with Shiny ",
                                               utils::packageVersion("shiny"),
                                               " | R ", R.version$major, ".", R.version$minor))
                    ) # div about-section
    ) # tabPanel About

  ) # tabsetPanel main_tabs
) # fluidPage


################################################################################
# SERVER
################################################################################

server <- function(input, output, session) {

  ##############################################################################
  # SHARED: contig filters (server-side selectize)
  ##############################################################################

  if (has_class) {
    tc <- sort(unique(as.character(class_data[[transcript_col]])))
    shiny::updateSelectizeInput(session, "class_contig_filter",
                                choices = c("All" = "", tc), selected = "", server = TRUE)
    shiny::updateSelectizeInput(session, "polya_contig_filter",
                                choices = c("All" = "", tc), selected = "", server = TRUE)
    # Report transcript selector
    shiny::updateSelectizeInput(session, "rpt_transcripts",
                                choices = tc, selected = character(0), server = TRUE)
  }
  if (has_residue) {
    tc_r <- sort(unique(as.character(residue_data[[transcript_col]])))
    shiny::updateSelectizeInput(session, "residue_contig_filter",
                                choices = c("All" = "", tc_r), selected = "", server = TRUE)
  }

  .filter_by_contig <- function(data, contig_value) {
    if (is.null(contig_value) || contig_value == "") return(data)
    data[data[[transcript_col]] == contig_value, , drop = FALSE]
  }


  ##############################################################################
  # CLASSIFICATION TAB
  ##############################################################################

  if (has_class) {

    class_data_filtered <- shiny::reactive({
      .filter_by_contig(class_data, input$class_contig_filter)
    })

    output$class_plot <- shiny::renderPlot({
      cd <- class_data_filtered(); shiny::req(nrow(cd) > 0)
      shiny::req(input$class_grouping, input$class_plot_type,
                 !is.null(input$class_frequency))
      gf <- if (input$class_grouping %in% names(cd)) input$class_grouping else NA
      p <- ninetails::plot_class_counts(class_data = cd,
                                        grouping_factor = gf, frequency = input$class_frequency,
                                        type = input$class_plot_type)
      # Remove subtitle to align with adjacent plots in flex row
      p + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
    })

    if (has_residue) {
      output$nonA_abundance_plot <- shiny::renderPlot({
        rd <- .filter_by_contig(residue_data, input$class_contig_filter)
        shiny::req(nrow(rd) > 0, input$class_grouping)
        gf <- if (input$class_grouping %in% names(rd)) input$class_grouping else NA
        tryCatch(
          ninetails::plot_nonA_abundance(residue_data = rd,
                                         grouping_factor = gf),
          error = function(e) {
            ggplot2::ggplot() +
              ggplot2::annotate("text", x = 0.5, y = 0.5,
                                label = paste("plot_nonA_abundance:", e$message),
                                size = 4, color = "#cc0000") +
              ggplot2::theme_void()
          }
        )
      })
    }
  }


  ##############################################################################
  # RESIDUES TAB (with rug density plots)
  ##############################################################################

  if (has_residue) {

    residue_data_filtered <- shiny::reactive({
      .filter_by_contig(residue_data, input$residue_contig_filter)
    })

    output$residue_plot <- shiny::renderPlot({
      rd <- residue_data_filtered(); shiny::req(nrow(rd) > 0)
      shiny::req(input$residue_grouping, !is.null(input$residue_frequency),
                 !is.null(input$residue_by_read))
      gf <- if (input$residue_grouping %in% names(rd)) input$residue_grouping else NA
      ninetails::plot_residue_counts(residue_data = rd,
                                     grouping_factor = gf, frequency = input$residue_frequency,
                                     by_read = input$residue_by_read)
    })

    # Rug density plots: one row of 3 plots (C, G, U) per sample
    # Each base type is subsampled to max 1000 points per sample
    # to keep plots readable and rendering fast.

    # Sanitize sample name for use as Shiny output ID
    .safe_id <- function(x) gsub("[^A-Za-z0-9]", "_", x)

    # Build dynamic UI: one card per sample with 3 rug plots
    output$rug_plots_ui <- shiny::renderUI({
      rd <- residue_data_filtered(); shiny::req(nrow(rd) > 0)

      # Determine samples present in filtered data
      if ("sample_name" %in% names(rd)) {
        sample_names <- sort(unique(as.character(rd$sample_name)))
      } else {
        sample_names <- "all"
      }

      # Create UI elements for each sample
      ui_list <- lapply(sample_names, function(sn) {
        sid <- .safe_id(sn)
        shiny::tagList(
          shiny::tags$h5(
            style = "color: #2C7BB6; font-weight: 600; margin-top: 15px; margin-bottom: 8px;",
            if (sn == "all") "All data" else sn
          ),
          shiny::fluidRow(
            shiny::column(4, shiny::plotOutput(
              paste0("rug_", sid, "_C"), height = "350px")),
            shiny::column(4, shiny::plotOutput(
              paste0("rug_", sid, "_G"), height = "350px")),
            shiny::column(4, shiny::plotOutput(
              paste0("rug_", sid, "_U"), height = "350px"))
          ) # fluidRow
        ) # tagList
      }) # lapply sample_names

      # Register render functions for each sample x base combination
      for (sn in sample_names) {
        local({
          local_sn <- sn
          sid <- .safe_id(local_sn)

          for (base_type in c("C", "G", "U")) {
            local({
              local_base <- base_type
              output_id <- paste0("rug_", sid, "_", local_base)

              output[[output_id]] <- shiny::renderPlot({
                rd <- residue_data_filtered(); shiny::req(nrow(rd) > 0)
                shiny::req(input$rug_max_length)

                # Filter by sample
                if (local_sn != "all" && "sample_name" %in% names(rd)) {
                  rd <- rd[rd$sample_name == local_sn, , drop = FALSE]
                }

                # Check if this base type exists for this sample
                if (!local_base %in% rd$prediction || nrow(rd) == 0) {
                  return(
                    ggplot2::ggplot() +
                      ggplot2::annotate("text", x = 0.5, y = 0.5,
                                        label = paste("No", local_base, "residues"),
                                        size = 4, color = "#999999") +
                      ggplot2::theme_void() +
                      ggplot2::ggtitle(local_base)
                  )
                }

                # Subsample to max 1000 points for this base type
                rd_base <- rd[rd$prediction == local_base, , drop = FALSE]
                n_avail <- nrow(rd_base)
                if (n_avail > 1000) {
                  set.seed(42)
                  keep_idx <- sample.int(n_avail, size = 1000, replace = FALSE)
                  # Rebuild rd: keep all non-target rows + subsampled target rows
                  rd_other <- rd[rd$prediction != local_base, , drop = FALSE]
                  rd <- rbind(rd_other, rd_base[keep_idx, ])
                }

                tryCatch(
                  ninetails::plot_rug_density(
                    residue_data = rd,
                    base = local_base,
                    max_length = input$rug_max_length
                  ),
                  error = function(e) {
                    ggplot2::ggplot() +
                      ggplot2::annotate("text", x = 0.5, y = 0.5,
                                        label = paste(local_base, "error:", e$message),
                                        size = 3.5, color = "#cc0000") +
                      ggplot2::theme_void()
                  }
                )
              }) # renderPlot
            }) # local base
          } # for base_type
        }) # local sample
      } # for sample_names

      do.call(shiny::tagList, ui_list)
    }) # renderUI rug_plots_ui

    if (has_merged) {
      output$residue_summary_table <- DT::renderDT({
        shiny::req(input$residue_grouping)
        md <- .filter_by_contig(merged_data, input$residue_contig_filter)
        shiny::req(nrow(md) > 0)
        tid_col <- if ("ensembl_transcript_id_short" %in% names(md)) {
          "ensembl_transcript_id_short"
        } else { NULL }
        gf <- if (input$residue_grouping %in% names(md)) {
          input$residue_grouping
        } else { "contig" }
        tryCatch({
          tbl <- ninetails::summarize_nonA(merged_nonA_tables = md,
                                           summary_factors = gf, transcript_id_column = tid_col)
          tbl <- dplyr::mutate(tbl, dplyr::across(
            tidyselect::vars_select_helpers$where(is.numeric), ~ round(., 3)))
          DT::datatable(tbl, options = list(scrollX = TRUE, pageLength = 15),
                        rownames = FALSE)
        }, error = function(e) {
          DT::datatable(data.frame(Message = paste("Error:", e$message)))
        })
      })
    }
  }


  ##############################################################################
  # POLY(A) DISTRIBUTION TAB
  ##############################################################################

  if (has_class) {

    # Dynamic UI: condition selector based on grouping variable
    output$polya_condition_filter_ui <- shiny::renderUI({
      shiny::req(input$polya_grouping)
      base <- if (has_merged) merged_data else class_data
      if (!input$polya_grouping %in% names(base)) return(NULL)
      conditions <- sort(unique(as.character(base[[input$polya_grouping]])))
      shiny::selectizeInput("polya_condition_filter",
                            "Select condition(s)",
                            choices = conditions,
                            selected = conditions,
                            multiple = TRUE,
                            options = list(placeholder = "All conditions",
                                           plugins = list("remove_button")))
    })

    # Reset button: restore all filters to defaults
    shiny::observeEvent(input$polya_reset, {
      shiny::updateSelectInput(session, "polya_grouping",
                               selected = available_groups[1])
      shiny::updateSelectizeInput(session, "polya_contig_filter",
                                  selected = "")
      shiny::updateSelectInput(session, "polya_center",
                               selected = "none")
      shiny::updateSliderInput(session, "polya_max_length",
                               value = 200)
      shiny::updateCheckboxInput(session, "polya_ndensity",
                                 value = TRUE)
      shiny::updateSelectInput(session, "polya_palette",
                               selected = "unova")
      # Condition filter is rebuilt by renderUI when grouping changes
    })

    polya_data_filtered <- shiny::reactive({
      base <- if (has_merged) merged_data else class_data
      # Apply contig filter
      base <- .filter_by_contig(base, input$polya_contig_filter)
      # Apply condition filter (only selected levels of grouping variable)
      if (!is.null(input$polya_condition_filter) &&
          !is.null(input$polya_grouping) &&
          input$polya_grouping %in% names(base)) {
        base <- base[base[[input$polya_grouping]] %in%
                       input$polya_condition_filter, , drop = FALSE]
      }
      base
    })

    output$polya_dist_plot <- shiny::renderPlot({
      pd <- polya_data_filtered(); shiny::req(nrow(pd) > 0)
      shiny::req(input$polya_grouping, input$polya_center, input$polya_palette)
      center_val <- if (input$polya_center == "none") NA else input$polya_center
      gf <- if (input$polya_grouping %in% names(pd)) input$polya_grouping else NA
      p <- ninetails::plot_tail_distribution(input_data = pd,
                                             grouping_factor = gf, max_length = input$polya_max_length,
                                             value_to_show = center_val, ndensity = input$polya_ndensity)
      pal_colors <- palettes[[input$polya_palette]]
      p + ggplot2::scale_color_manual(values = pal_colors)
    })

    output$polya_summary_table <- DT::renderDT({
      pd <- polya_data_filtered(); shiny::req(nrow(pd) > 0)
      shiny::req(input$polya_grouping)
      gf <- if (input$polya_grouping %in% names(pd)) input$polya_grouping else NULL

      # Detect poly(A) length column
      len_col <- if ("polya_length" %in% names(pd)) {
        "polya_length"
      } else if ("poly_tail_length" %in% names(pd)) {
        "poly_tail_length"
      } else { NULL }
      shiny::req(len_col)

      # Compute per-condition summary
      if (!is.null(gf)) {
        tbl <- pd %>%
          dplyr::group_by(!!rlang::sym(gf)) %>%
          dplyr::summarise(
            n_reads = dplyr::n(),
            mean_length = mean(!!rlang::sym(len_col), na.rm = TRUE),
            median_length = stats::median(!!rlang::sym(len_col), na.rm = TRUE),
            sd = stats::sd(!!rlang::sym(len_col), na.rm = TRUE),
            sem = stats::sd(!!rlang::sym(len_col), na.rm = TRUE) /
              sqrt(dplyr::n()),
            .groups = "drop"
          )
      } else {
        tbl <- pd %>%
          dplyr::summarise(
            n_reads = dplyr::n(),
            mean_length = mean(!!rlang::sym(len_col), na.rm = TRUE),
            median_length = stats::median(!!rlang::sym(len_col), na.rm = TRUE),
            sd = stats::sd(!!rlang::sym(len_col), na.rm = TRUE),
            sem = stats::sd(!!rlang::sym(len_col), na.rm = TRUE) /
              sqrt(dplyr::n())
          )
      }

      tbl <- dplyr::mutate(tbl, dplyr::across(
        tidyselect::vars_select_helpers$where(is.numeric), ~ round(., 2)))

      DT::datatable(tbl,
                    options = list(dom = "t", paging = FALSE, searching = FALSE,
                                   ordering = FALSE),
                    rownames = FALSE,
                    colnames = c(if (!is.null(gf)) gf else NULL,
                                 "Reads", "Mean length", "Median length", "SD", "SEM"))
    })
  }


  ##############################################################################
  # DOWNLOAD TAB
  ##############################################################################

  if (has_class) {
    output$download_report <- shiny::downloadHandler(
      filename = function() {
        paste0("ninetails_report_", format(Sys.Date(), "%Y%m%d"), ".html")
      },
      content = function(file) {

        shiny::withProgress(message = "Generating report...", value = 0, {

          tmp_dir <- tempdir()
          gf <- if ("sample_name" %in% names(class_data)) "sample_name" else NA

          # Helper: save a ggplot as base64 PNG string
          .plot_to_base64 <- function(p, w = 10, h = 6) {
            png_file <- tempfile(tmpdir = tmp_dir, fileext = ".png")
            tryCatch({
              ggplot2::ggsave(png_file, p, width = w, height = h, dpi = 150)
              b64 <- base64enc::dataURI(file = png_file, mime = "image/png")
              unlink(png_file)
              b64
            }, error = function(e) { unlink(png_file); NULL })
          }

          # Helper: add a section with a ggplot
          .add_plot_section <- function(html, title, plot_expr, w = 10, h = 6) {
            p <- tryCatch(plot_expr, error = function(e) NULL)
            if (!is.null(p)) {
              b64 <- .plot_to_base64(p, w, h)
              if (!is.null(b64)) {
                html <- c(html,
                          paste0("<h2>", title, "</h2>"),
                          paste0("<img src='", b64, "'>"))
              }
            }
            html
          }

          # Read checkbox states (FALSE if checkbox is absent from UI)
          want_class     <- isTRUE(input$rpt_classification)
          want_abundance <- isTRUE(input$rpt_abundance) && has_residue
          want_residue   <- isTRUE(input$rpt_residue) && has_residue
          want_rug       <- isTRUE(input$rpt_rug) && has_residue
          want_polya     <- isTRUE(input$rpt_polya)
          want_signals   <- isTRUE(input$rpt_signals) && length(signal_config) > 0

          # Report-specific settings
          rpt_center_val <- if (!is.null(input$rpt_center) &&
                                input$rpt_center != "none") {
            input$rpt_center
          } else { NA }
          rpt_max_len <- input$rpt_max_length %||% 200
          rpt_palette <- input$rpt_palette %||% "unova"
          rpt_transcripts <- input$rpt_transcripts  # NULL or character vector

          # ---- HTML header ----
          html <- c(
            "<!DOCTYPE html><html><head>",
            "<meta charset='UTF-8'>",
            paste0("<title>Ninetails Report - ", Sys.Date(), "</title>"),
            "<style>",
            "body { font-family: 'Helvetica Neue', sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; }",
            "h1 { color: #2C7BB6; } h2 { color: #333; border-bottom: 2px solid #dce2e8; padding-bottom: 6px; }",
            "h3 { color: #2C7BB6; margin-top: 30px; }",
            "h4 { color: #555; margin-top: 20px; }",
            "img { max-width: 100%; margin: 10px 0; display: block; }",
            ".stat { display: inline-block; padding: 10px 20px; margin: 5px; background: #f0f5fa; border-radius: 6px; border-left: 4px solid #2C7BB6; }",
            ".stat b { font-size: 20px; color: #2C7BB6; }",
            ".note { background: #f5f7f9; padding: 8px 12px; border-radius: 4px; font-size: 13px; color: #666; margin: 5px 0 15px 0; }",
            ".absent { background: #f5f5f5; padding: 20px; text-align: center; color: #999; border-radius: 6px; margin: 5px 0; font-size: 13px; }",
            ".setting { font-size: 12px; color: #888; font-style: italic; }",
            "</style></head><body>",
            "<h1>Ninetails Analysis Report</h1>",
            paste0("<p>Generated: ", Sys.time(), " | Package version: ",
                   utils::packageVersion("ninetails"), "</p>"),
            "<h2>Summary</h2>",
            paste0("<div class='stat'><b>", n_samples_val, "</b><br>Samples</div>"),
            paste0("<div class='stat'><b>", format(n_reads_val, big.mark = ","),
                   "</b><br>Total reads</div>"),
            paste0("<div class='stat'><b>", format(n_blank_val, big.mark = ","),
                   "</b><br>Blank</div>"),
            paste0("<div class='stat'><b>", format(n_decorated_val, big.mark = ","),
                   "</b><br>Decorated</div>"),
            paste0("<div class='stat'><b>", n_transcripts_val,
                   "</b><br>Transcripts</div>")
          )

          shiny::incProgress(0.05, detail = "Global plots...")

          # ---- Report descriptions (reusable) ----
          desc_class <- paste0(
            "Read classification summary. Reads are categorized as ",
            "<b>decorated</b> (containing non-adenosine residues), ",
            "<b>blank</b> (no modifications detected), or ",
            "<b>unclassified</b> (insufficient data quality). ",
            "Classification codes: ",
            "<b>YAY</b> \u2014 move transition present, non-A detected; ",
            "<b>MPU</b> \u2014 move transition present, only A detected; ",
            "<b>MAU</b> \u2014 no move transition, unmodified tail; ",
            "<b>IRL</b> \u2014 insufficient read length (<10 nt); ",
            "<b>BAC</b> \u2014 bad coordinates; ",
            "<b>UNM</b> \u2014 unmapped read.")
          desc_abundance <- paste0(
            "Frequency of reads containing one, two, or three or more ",
            "separate non-adenosine residues per read. Computed relative ",
            "to the total number of decorated reads in each sample.")
          desc_residue <- paste0(
            "Distribution of non-adenosine residue types (cytidine, ",
            "guanosine, uridine) detected across samples.")
          desc_rug <- paste0(
            "Positional distribution of non-adenosine residues along ",
            "poly(A) tails. Each point represents the estimated position ",
            "(from the 3\u2019 end) of a detected modification. Marginal ",
            "density curves show the overall distribution shape.")
          desc_polya <- paste0(
            "Density distribution of poly(A) tail lengths across samples ",
            "or experimental conditions.")

          # ---- Global classification ----
          if (want_class) {
            html <- .add_plot_section(html, "Read Classification",
                                      ninetails::plot_class_counts(class_data = class_data,
                                                                   grouping_factor = gf, frequency = FALSE, type = "N") +
                                        ggplot2::theme(plot.subtitle = ggplot2::element_blank()))
            html <- c(html, paste0("<div class='note'>", desc_class, "</div>"))
          }

          # ---- Global abundance ----
          if (want_abundance) {
            html <- .add_plot_section(html, "Non-A Abundance",
                                      ninetails::plot_nonA_abundance(residue_data = residue_data,
                                                                     grouping_factor = gf))
            html <- c(html, paste0("<div class='note'>", desc_abundance, "</div>"))
          }

          # ---- Global residue frequency ----
          if (want_residue) {
            html <- .add_plot_section(html, "Residue Frequency",
                                      ninetails::plot_residue_counts(residue_data = residue_data,
                                                                     grouping_factor = gf, frequency = TRUE, by_read = FALSE))
            html <- c(html, paste0("<div class='note'>", desc_residue, "</div>"))
          }

          # ---- Global rug density plots ----
          if (want_rug) {
            shiny::incProgress(0.10, detail = "Rug density plots...")

            html <- c(html, "<h2>Non-A Position Distribution</h2>",
                      paste0("<div class='note'>", desc_rug, "</div>"))

            sample_names_rug <- if ("sample_name" %in% names(residue_data)) {
              sort(unique(as.character(residue_data$sample_name)))
            } else { "all" }

            for (sn in sample_names_rug) {
              rd_s <- if (sn == "all") residue_data else {
                residue_data[residue_data$sample_name == sn, , drop = FALSE]
              }
              if (nrow(rd_s) == 0) next
              html <- c(html, paste0("<h3>",
                                     if (sn == "all") "All data" else sn, "</h3>"))

              for (base in c("C", "G", "U")) {
                n_total <- sum(rd_s$prediction == base, na.rm = TRUE)
                if (n_total == 0) {
                  html <- c(html, paste0(
                    "<div class='absent'>No ", base, " residues detected</div>"))
                  next
                }
                downsampled <- FALSE; rd_plot <- rd_s
                if (n_total > 1000) {
                  set.seed(42)
                  rd_base <- rd_s[rd_s$prediction == base, , drop = FALSE]
                  rd_other <- rd_s[rd_s$prediction != base, , drop = FALSE]
                  keep <- sample.int(n_total, 1000, replace = FALSE)
                  rd_plot <- rbind(rd_other, rd_base[keep, ])
                  downsampled <- TRUE
                }
                tryCatch({
                  p <- ninetails::plot_rug_density(residue_data = rd_plot,
                                                   base = base, max_length = rpt_max_len)
                  b64 <- .plot_to_base64(p, 5.5, 3.5)
                  if (!is.null(b64)) {
                    html <- c(html, paste0("<img src='", b64, "'>"))
                    if (downsampled) {
                      html <- c(html, paste0("<div class='note'>", base,
                                             " positions downsampled from ",
                                             format(n_total, big.mark = ","),
                                             " to 1,000 points.</div>"))
                    } else {
                      html <- c(html, paste0("<div class='note'>", base, ": ",
                                             format(n_total, big.mark = ","),
                                             " points (no downsampling).</div>"))
                    }
                  }
                }, error = function(e) {
                  html <<- c(html, paste0("<div class='absent'>", base,
                                          " rug error: ", e$message, "</div>"))
                })
              } # for base
            } # for sn
          } # want_rug

          # ---- Global poly(A) distribution ----
          if (want_polya) {
            shiny::incProgress(0.05, detail = "Poly(A) distribution...")

            base_data <- if (has_merged) merged_data else class_data
            tryCatch({
              p <- ninetails::plot_tail_distribution(input_data = base_data,
                                                     grouping_factor = gf, max_length = rpt_max_len,
                                                     value_to_show = rpt_center_val, ndensity = TRUE)
              p <- p + ggplot2::scale_color_manual(
                values = palettes[[rpt_palette]])
              b64 <- .plot_to_base64(p)
              if (!is.null(b64)) {
                html <- c(html, "<h2>Poly(A) Length Distribution</h2>")
                # Annotate settings
                center_label <- if (is.na(rpt_center_val)) "none" else rpt_center_val
                html <- c(html, paste0(
                  "<p class='setting'>Central tendency: ", center_label,
                  " | Max length: ", rpt_max_len,
                  " | Palette: ", rpt_palette, "</p>"))
                html <- c(html, paste0("<img src='", b64, "'>"),
                          paste0("<div class='note'>", desc_polya, "</div>"))
              }
            }, error = function(e) NULL)
          }

          # ---- Per-transcript sections ----
          if (!is.null(rpt_transcripts) && length(rpt_transcripts) > 0) {
            shiny::incProgress(0.10, detail = "Per-transcript sections...")

            for (tr in rpt_transcripts) {
              html <- c(html, paste0(
                "<h2 style='color: #2C7BB6;'>Transcript: ", tr, "</h2>"))

              # Filter data to this transcript
              cd_tr <- class_data[class_data[[transcript_col]] == tr, ,
                                  drop = FALSE]
              if (nrow(cd_tr) == 0) {
                html <- c(html, paste0(
                  "<div class='absent'>No reads for transcript ", tr, "</div>"))
                next
              }

              html <- c(html, paste0("<div class='note'>",
                                     format(nrow(cd_tr), big.mark = ","), " reads</div>"))

              # Classification
              if (want_class) {
                html <- .add_plot_section(html,
                                          paste0("Classification \u2014 ", tr),
                                          ninetails::plot_class_counts(class_data = cd_tr,
                                                                       grouping_factor = gf, frequency = FALSE, type = "N") +
                                            ggplot2::theme(plot.subtitle = ggplot2::element_blank()),
                                          w = 8, h = 5)
              }

              # Residue frequency
              if (want_residue && has_residue) {
                rd_tr <- residue_data[residue_data[[transcript_col]] == tr, ,
                                      drop = FALSE]
                if (nrow(rd_tr) > 0) {
                  html <- .add_plot_section(html,
                                            paste0("Residue Frequency \u2014 ", tr),
                                            ninetails::plot_residue_counts(residue_data = rd_tr,
                                                                           grouping_factor = gf, frequency = TRUE, by_read = FALSE),
                                            w = 8, h = 5)
                }
              }

              # Poly(A) distribution
              if (want_polya) {
                tr_data <- if (has_merged) {
                  merged_data[merged_data[[transcript_col]] == tr, ,
                              drop = FALSE]
                } else { cd_tr }
                if (nrow(tr_data) > 0) {
                  tryCatch({
                    p <- ninetails::plot_tail_distribution(
                      input_data = tr_data, grouping_factor = gf,
                      max_length = rpt_max_len, value_to_show = rpt_center_val,
                      ndensity = TRUE)
                    p <- p + ggplot2::scale_color_manual(
                      values = palettes[[rpt_palette]])
                    b64 <- .plot_to_base64(p, 8, 5)
                    if (!is.null(b64)) {
                      html <- c(html,
                                paste0("<h3>Poly(A) Length \u2014 ", tr, "</h3>"),
                                paste0("<img src='", b64, "'>"))
                    }
                  }, error = function(e) NULL)
                }
              }
            } # for tr
          } # per-transcript

          # ---- Signal plots (one per row, full width) ----
          if (want_signals) {
            shiny::incProgress(0.15, detail = "Signal plots...")

            html <- c(html,
                      "<h2>Example Tail Signals (5 random reads per category)</h2>")

            n_sig <- length(signal_config)
            sig_step <- 0.3 / max(n_sig, 1)

            for (sig_sn in names(signal_config)) {
              sc <- signal_config[[sig_sn]]
              if (!file.exists(sc$dorado_summary) || !dir.exists(sc$pod5_dir)) {
                html <- c(html, paste0("<h3>", sig_sn, "</h3>",
                                       "<div class='absent'>Signal files not accessible</div>"))
                next
              }

              html <- c(html, paste0("<h3>", sig_sn, "</h3>"))
              shiny::incProgress(sig_step,
                                 detail = paste0("Signals: ", sig_sn, "..."))

              dorado_sum <- tryCatch({
                ds <- vroom::vroom(sc$dorado_summary, show_col_types = FALSE)
                if (!"read_id" %in% names(ds) && "readname" %in% names(ds))
                  ds <- dplyr::rename(ds, read_id = readname)
                ds
              }, error = function(e) NULL)

              if (is.null(dorado_sum)) {
                html <- c(html,
                          "<div class='absent'>Could not load dorado summary</div>")
                next
              }

              # Get class/residue data for this sample
              cd_s <- if (sig_sn != "single" &&
                          "sample_name" %in% names(class_data)) {
                class_data[class_data$sample_name == sig_sn, , drop = FALSE]
              } else { class_data }

              rd_s <- if (has_residue && sig_sn != "single" &&
                          "sample_name" %in% names(residue_data)) {
                residue_data[residue_data$sample_name == sig_sn, , drop = FALSE]
              } else if (has_residue) { residue_data } else { NULL }

              # Define categories and their read IDs
              categories <- list()
              class_rn_col <- if ("readname" %in% names(cd_s)) {
                "readname"
              } else if ("read_id" %in% names(cd_s)) {
                "read_id"
              } else { NULL }

              if (!is.null(class_rn_col)) {
                blank_ids <- cd_s[[class_rn_col]][
                  cd_s$class == "blank" & !is.na(cd_s$class)]
                blank_ids <- blank_ids[!is.na(blank_ids)]
                if (length(blank_ids) > 0) categories$Blank <- blank_ids
              }

              if (!is.null(rd_s) && nrow(rd_s) > 0) {
                rn_col <- if ("readname" %in% names(rd_s)) {
                  "readname"
                } else { "read_id" }
                for (base in c("C", "G", "U")) {
                  bids <- unique(rd_s[[rn_col]][rd_s$prediction == base])
                  if (length(bids) > 0)
                    categories[[paste0("Decorated (", base, ")")]] <- bids
                }
              }

              if (length(categories) == 0) {
                html <- c(html,
                          "<div class='absent'>No classified reads available</div>")
                next
              }

              for (cat_name in names(categories)) {
                cat_ids <- categories[[cat_name]]
                n_avail <- length(cat_ids)
                set.seed(42)
                n_pick <- min(5, n_avail)
                picked <- if (n_avail > 5) {
                  sample(cat_ids, 5, replace = FALSE)
                } else { cat_ids }

                html <- c(html, paste0("<h4>", cat_name,
                                       " (", n_pick, " of ", format(n_avail, big.mark = ","),
                                       " reads)</h4>"))

                for (rid in picked) {
                  tryCatch({
                    p <- ninetails::plot_tail_range_pod5(
                      readname = rid,
                      dorado_summary = dorado_sum,
                      workspace = sc$pod5_dir,
                      flank = 150, rescale = FALSE,
                      residue_data = rd_s,
                      nonA_flank = 250)
                    b64 <- .plot_to_base64(p, 10, 4)
                    if (!is.null(b64))
                      html <- c(html, paste0("<img src='", b64, "'>"))
                  }, error = function(e) {
                    html <<- c(html, paste0(
                      "<div class='absent'>Error for ",
                      substr(rid, 1, 12), "...: ", e$message, "</div>"))
                  })
                } # for rid
              } # for cat_name
            } # for sig_sn
          } # want_signals

          shiny::incProgress(0.95, detail = "Finalizing...")

          html <- c(html, "<hr>",
                    paste0("<p style='color: #999; font-size: 11px; text-align: center;'>",
                           "Generated by ninetails v",
                           utils::packageVersion("ninetails"),
                           " | ", Sys.time(), "</p>"),
                    "</body></html>")

          writeLines(paste(html, collapse = "\n"), file)

        }) # withProgress
      },
      contentType = "text/html"
    ) # downloadHandler
  }


  ##############################################################################
  # SIGNAL VIEWER TAB
  ##############################################################################

  loaded_signal_data    <- shiny::reactiveVal(NULL)
  loaded_signal_residue <- shiny::reactiveVal(NULL)
  signal_error          <- shiny::reactiveVal(NULL)

  # ---- Data loading ----

  if (length(signal_config) > 0 && names(signal_config)[1] != "single") {
    # Multi mode with named samples
    shiny::observeEvent(input$signal_sample, {
      s <- signal_config[[input$signal_sample]]
      if (is.null(s)) { signal_error("Sample not found."); loaded_signal_data(NULL); return() }
      if (!file.exists(s$dorado_summary)) {
        signal_error(paste("Not found:", s$dorado_summary))
        loaded_signal_data(NULL); return()
      }
      tryCatch({
        data <- vroom::vroom(s$dorado_summary, show_col_types = FALSE)
        if (!"read_id" %in% names(data) && "readname" %in% names(data))
          data <- dplyr::rename(data, read_id = readname)
        signal_error(NULL); loaded_signal_data(data)
      }, error = function(e) {
        signal_error(paste("Error:", e$message)); loaded_signal_data(NULL)
      })

      # Residue overlay: use pre-loaded residue_data from config
      if (has_residue) {
        loaded_signal_residue(residue_data)
      } else { loaded_signal_residue(NULL) }
    }, ignoreNULL = FALSE)

  } else if (length(signal_config) > 0 && names(signal_config)[1] == "single") {
    # Single mode with pre-filled paths: auto-load
    shiny::observe({
      s <- signal_config[["single"]]
      if (is.null(s) || !file.exists(s$dorado_summary)) {
        signal_error("Summary file not found.")
        loaded_signal_data(NULL); return()
      }
      tryCatch({
        data <- vroom::vroom(s$dorado_summary, show_col_types = FALSE)
        if (!"read_id" %in% names(data) && "readname" %in% names(data))
          data <- dplyr::rename(data, read_id = readname)
        signal_error(NULL); loaded_signal_data(data)
      }, error = function(e) {
        signal_error(paste("Error:", e$message)); loaded_signal_data(NULL)
      })

      rf <- trimws(default_residue)
      if (nzchar(rf) && file.exists(rf)) {
        tryCatch({
          rd <- vroom::vroom(rf, show_col_types = FALSE)
          if (!"read_id" %in% names(rd) && "readname" %in% names(rd))
            rd <- dplyr::rename(rd, read_id = readname)
          loaded_signal_residue(rd)
        }, error = function(e) { loaded_signal_residue(NULL) })
      } else if (has_residue) {
        loaded_signal_residue(residue_data)
      } else { loaded_signal_residue(NULL) }
    })

  } else {
    # No signal config: manual load
    shiny::observeEvent(input$load_signal_data, {
      sf <- trimws(input$signal_summary_file)
      pd <- trimws(input$signal_pod5_dir)
      rf <- trimws(input$signal_residue_file)

      if (!nzchar(sf)) { signal_error("Enter summary file."); loaded_signal_data(NULL); return() }
      if (!file.exists(sf)) { signal_error(paste("Not found:", sf)); loaded_signal_data(NULL); return() }
      if (!nzchar(pd)) { signal_error("Enter POD5 directory."); loaded_signal_data(NULL); return() }
      if (!dir.exists(pd)) { signal_error(paste("Not found:", pd)); loaded_signal_data(NULL); return() }

      tryCatch({
        data <- vroom::vroom(sf, show_col_types = FALSE)
        if (!"read_id" %in% names(data) && "readname" %in% names(data))
          data <- dplyr::rename(data, read_id = readname)
        signal_error(NULL); loaded_signal_data(data)
      }, error = function(e) {
        signal_error(paste("Error:", e$message)); loaded_signal_data(NULL)
      })

      if (nzchar(rf) && file.exists(rf)) {
        tryCatch({
          rd <- vroom::vroom(rf, show_col_types = FALSE)
          if (!"read_id" %in% names(rd) && "readname" %in% names(rd))
            rd <- dplyr::rename(rd, read_id = readname)
          loaded_signal_residue(rd)
        }, error = function(e) { loaded_signal_residue(NULL) })
      } else { loaded_signal_residue(NULL) }
    })
  }

  # ---- Active pod5 dir ----
  active_pod5_dir <- shiny::reactive({
    if (length(signal_config) > 0 && names(signal_config)[1] != "single") {
      s <- signal_config[[input$signal_sample]]
      if (!is.null(s)) s$pod5_dir else ""
    } else if (length(signal_config) > 0 && names(signal_config)[1] == "single") {
      signal_config[["single"]]$pod5_dir %||% ""
    } else {
      trimws(input$signal_pod5_dir %||% "")
    }
  })

  # ---- Status display ----
  output$signal_data_status <- shiny::renderUI({
    err <- signal_error(); data <- loaded_signal_data()
    if (!is.null(err))
      shiny::tags$div(style = "background: #f8d7da; color: #721c24; padding: 10px; border-radius: 4px; margin-top: 10px;",
                      shiny::tags$strong("Error: "), err)
    else if (!is.null(data))
      shiny::tags$div(style = "background: #d4edda; color: #155724; padding: 10px; border-radius: 4px; margin-top: 10px;",
                      shiny::tags$strong("Loaded: "), format(nrow(data), big.mark = ","), " reads")
  })

  # ---- Dynamic filters ----
  output$genome_filter_ui <- shiny::renderUI({
    data <- loaded_signal_data(); shiny::req(data)
    if ("alignment_genome" %in% names(data)) {
      genomes <- sort(as.character(unique(stats::na.omit(data$alignment_genome))))
      shiny::selectInput("genome_filter", "Aligned transcript",
                         choices = c("All", genomes), selected = "All")
    }
  })
  output$mapq_filter_ui <- shiny::renderUI({
    data <- loaded_signal_data(); shiny::req(data)
    if ("alignment_mapq" %in% names(data)) {
      rng <- range(data$alignment_mapq, na.rm = TRUE)
      shiny::sliderInput("mapq_filter", "Mapping quality range",
                         min = rng[1], max = rng[2], value = rng, step = 1)
    }
  })

  # ---- Filtered reads ----
  sig_filtered_data <- shiny::reactive({
    data <- loaded_signal_data(); shiny::req(data)
    if (!is.null(input$polya_length_range) && "poly_tail_length" %in% names(data))
      data <- dplyr::filter(data,
                            poly_tail_length >= input$polya_length_range[1],
                            poly_tail_length <= input$polya_length_range[2])
    if (!is.null(input$comments_filter) && input$comments_filter != "All" && "comments" %in% names(data))
      data <- dplyr::filter(data, comments == input$comments_filter)
    if (!is.null(input$genome_filter) && input$genome_filter != "All" && "alignment_genome" %in% names(data))
      data <- dplyr::filter(data, alignment_genome == input$genome_filter)
    if (!is.null(input$mapq_filter) && "alignment_mapq" %in% names(data))
      data <- dplyr::filter(data, alignment_mapq >= input$mapq_filter[1], alignment_mapq <= input$mapq_filter[2])
    rt <- loaded_signal_residue()
    if (!is.null(input$sig_residue_filter) && input$sig_residue_filter != "All" && !is.null(rt)) {
      id_col <- if ("read_id" %in% names(rt)) "read_id" else if ("readname" %in% names(rt)) "readname" else NULL
      if (!is.null(id_col)) {
        rids <- unique(rt[[id_col]][rt$prediction == input$sig_residue_filter])
        data <- dplyr::filter(data, read_id %in% rids)
      }
    }
    dplyr::filter(data, !is.na(poly_tail_start), !is.na(poly_tail_end),
                  poly_tail_start > 0, poly_tail_end > poly_tail_start)
  })

  # ---- Read selection ----
  output$read_selection_ui <- shiny::renderUI({
    data <- sig_filtered_data()
    if (is.null(data) || nrow(data) == 0) return(shiny::helpText("No reads match filters."))
    shiny::tagList(
      shiny::selectizeInput("selected_read_id", "Select Read ID", choices = NULL,
                            options = list(placeholder = "Type to search...", maxOptions = 100)),
      shiny::helpText("Start typing to search"),
      shiny::uiOutput("read_nav_buttons"))
  })
  shiny::observe({
    data <- sig_filtered_data()
    if (!is.null(data) && nrow(data) > 0)
      shiny::updateSelectizeInput(session, "selected_read_id",
                                  choices = as.character(data$read_id),
                                  selected = as.character(data$read_id[1]), server = TRUE)
  })

  # ---- Navigation ----
  output$read_nav_buttons <- shiny::renderUI({
    data <- sig_filtered_data()
    if (is.null(data) || nrow(data) == 0) return(NULL)
    ids <- as.character(data$read_id)
    idx <- match(input$selected_read_id, ids); if (is.na(idx)) idx <- 1
    shiny::div(style = "display: flex; gap: 8px; margin-top: 10px;",
               shiny::tags$button(id = "prev_read", type = "button",
                                  class = paste("btn btn-default action-button", if (idx <= 1) "disabled"),
                                  disabled = if (idx <= 1) "disabled", shiny::HTML("&larr; Previous")),
               shiny::tags$button(id = "next_read", type = "button",
                                  class = paste("btn btn-default action-button", if (idx >= length(ids)) "disabled"),
                                  disabled = if (idx >= length(ids)) "disabled", shiny::HTML("Next &rarr;")),
               shiny::tags$span(style = "align-self: center; margin-left: 8px; color: #666; font-size: 12px;",
                                paste0("(", idx, " / ", length(ids), ")")))
  })
  shiny::observeEvent(input$prev_read, {
    data <- sig_filtered_data(); shiny::req(data, nrow(data) > 0)
    ids <- as.character(data$read_id); idx <- match(input$selected_read_id, ids)
    if (!is.na(idx) && idx > 1) shiny::updateSelectizeInput(session, "selected_read_id", selected = ids[idx - 1])
  }, ignoreInit = TRUE)
  shiny::observeEvent(input$next_read, {
    data <- sig_filtered_data(); shiny::req(data, nrow(data) > 0)
    ids <- as.character(data$read_id); idx <- match(input$selected_read_id, ids)
    if (!is.na(idx) && idx < length(ids)) shiny::updateSelectizeInput(session, "selected_read_id", selected = ids[idx + 1])
  }, ignoreInit = TRUE)

  output$filter_summary <- shiny::renderUI({
    data <- sig_filtered_data(); total <- loaded_signal_data()
    if (is.null(data) || is.null(total)) return(NULL)
    shiny::tags$p(style = "color: #666; font-size: 13px; margin-top: 10px;",
                  paste0("Showing ", nrow(data), " of ", nrow(total), " reads"))
  })

  # ---- Selected read + residue ----
  selected_read <- shiny::reactive({
    data <- sig_filtered_data(); shiny::req(data, input$selected_read_id)
    ri <- dplyr::filter(data, read_id == input$selected_read_id)
    if (nrow(ri) == 0) NULL else ri[1, ]
  })
  selected_residue <- shiny::reactive({
    rd <- loaded_signal_residue(); shiny::req(input$selected_read_id)
    if (is.null(rd)) return(NULL)
    if (!"read_id" %in% names(rd) && "readname" %in% names(rd))
      rd <- dplyr::rename(rd, read_id = readname)
    if (!"read_id" %in% names(rd)) return(NULL)
    rr <- rd[rd$read_id == input$selected_read_id, , drop = FALSE]
    if (nrow(rr) == 0) NULL else rr
  })

  # ---- Signal extraction ----
  signal_data <- shiny::reactive({
    ri <- selected_read(); shiny::req(ri)
    pd <- active_pod5_dir()
    if (!nzchar(pd) || !dir.exists(pd)) {
      shiny::showNotification(paste("POD5 dir not found:", pd), type = "error"); return(NULL)
    }
    fn <- if ("filename" %in% names(ri)) ri$filename else NULL
    pf <- ninetails:::.find_pod5_file(fn, pd)
    if (is.null(pf)) { shiny::showNotification("No POD5 file found", type = "error"); return(NULL) }

    shiny::withProgress(message = "Extracting signal...", value = 0.5, {
      result <- ninetails:::.extract_signal_pod5(read_id = ri$read_id, pod5_file = pf, winsorize = TRUE)
    })
    if (is.null(result)) { shiny::showNotification("Extraction failed", type = "error"); return(NULL) }

    sig <- result$signal; ps <- ri$poly_tail_start; pe <- ri$poly_tail_end
    if (length(sig) == 0) return(NULL)
    df <- data.frame(position = seq_along(sig), signal = sig)
    df$segment <- "Other"
    if (!is.na(ps) && !is.na(pe) && ps > 0 && pe > ps && pe <= length(sig)) {
      df$segment[df$position < ps] <- "Adapter"
      df$segment[df$position >= ps & df$position <= pe] <- "Poly(A)"
      df$segment[df$position > pe] <- "Transcript"
    }
    df$segment <- factor(df$segment, levels = c("Adapter", "Poly(A)", "Transcript", "Other"))
    attr(df, "read_id") <- ri$read_id; attr(df, "polya_start") <- ps; attr(df, "polya_end") <- pe
    if ("poly_tail_length" %in% names(ri)) attr(df, "polya_length") <- ri$poly_tail_length
    if ("comments" %in% names(ri)) attr(df, "comments") <- ri$comments
    if ("alignment_genome" %in% names(ri)) attr(df, "alignment_genome") <- ri$alignment_genome
    df
  })

  output$signal_loaded <- shiny::reactive({ !is.null(signal_data()) })
  shiny::outputOptions(output, "signal_loaded", suspendWhenHidden = FALSE)

  # ---- Read info helper ----
  .render_read_info <- function(df, read_res) {
    cc <- attr(df, "comments"); ct <- NULL
    if (!is.null(cc) && length(cc) == 1 && !is.na(cc)) {
      ct <- switch(cc,
                   "YAY" = "Decorated: Pseudomove present, non-A residue detected",
                   "MPU" = "Blank: Pseudomove present, but only A detected",
                   "MAU" = "Blank: No pseudomove, empty tail",
                   "IRL" = "Unclassified: Tail too short (<10 nt)",
                   "BAC" = "Unclassified: Bad coordinates",
                   "UNM" = "Unclassified: Unmapped read", NULL)
    }
    res_ui <- NULL
    if (!is.null(read_res) && nrow(read_res) > 0) {
      rs <- paste(vapply(seq_len(nrow(read_res)), function(i) {
        paste0(read_res$prediction[i], " (pos ", read_res$est_nonA_pos[i], " from 3')")
      }, character(1)), collapse = "; ")
      res_ui <- shiny::tags$div(class = "residue-info",
                                shiny::tags$strong("Non-A residues: "), rs)
    }
    shiny::tags$div(class = "read-info",
                    shiny::tags$p(shiny::tags$strong("Read ID: "), shiny::tags$code(attr(df, "read_id"))),
                    shiny::tags$p(shiny::tags$strong("Signal length: "), format(nrow(df), big.mark = ","), " positions"),
                    shiny::tags$p(shiny::tags$strong("Poly(A) start: "), format(attr(df, "polya_start"), big.mark = ",")),
                    shiny::tags$p(shiny::tags$strong("Poly(A) end: "), format(attr(df, "polya_end"), big.mark = ",")),
                    if (!is.null(attr(df, "polya_length")))
                      shiny::tags$p(shiny::tags$strong("Poly(A) length: "), attr(df, "polya_length"), " nt"),
                    if (!is.null(cc)) shiny::tags$p(shiny::tags$strong("Status: "), cc),
                    if (!is.null(attr(df, "alignment_genome")))
                      shiny::tags$p(shiny::tags$strong("Genome: "), attr(df, "alignment_genome")),
                    if (!is.null(ct)) shiny::tags$div(class = "classification-info",
                                                      shiny::tags$strong("Classification: "), ct),
                    res_ui)
  }

  output$read_info <- shiny::renderUI({
    df <- signal_data(); shiny::req(df); .render_read_info(df, selected_residue()) })
  output$read_info_explorer <- shiny::renderUI({
    df <- signal_data(); shiny::req(df); .render_read_info(df, selected_residue()) })

  # ---- Static plots ----

  # Helper: build non-A overlay layers without relying on internal function
  # (avoids alpha/RGB issues on some server graphics devices)
  .safe_nonA_overlay <- function(read_nonA_data, poly_tail_start,
                                 poly_tail_end, nonA_flank = 250) {
    layers <- list()
    if (is.null(read_nonA_data) || nrow(read_nonA_data) == 0) return(layers)

    rc <- c("C" = "#3a424f", "G" = "#50a675", "U" = "#b0bdd4")

    # Get poly(A) length from residue data
    pl <- if ("polya_length" %in% names(read_nonA_data)) {
      read_nonA_data$polya_length[1]
    } else { NULL }
    if (is.null(pl) || is.na(pl) || pl <= 0) return(layers)

    for (i in seq_len(nrow(read_nonA_data))) {
      pred <- as.character(read_nonA_data$prediction[i])
      pos  <- read_nonA_data$est_nonA_pos[i]
      if (is.na(pos) || is.na(pred)) next

      fc <- rc[pred]
      if (is.na(fc)) next

      # Convert nucleotide position to raw signal position
      raw_pos <- poly_tail_start +
        (pl - pos) * (poly_tail_end - poly_tail_start) / pl
      x0 <- raw_pos - nonA_flank
      x1 <- raw_pos + nonA_flank

      layers <- c(layers, list(
        ggplot2::annotate("rect", xmin = x0, xmax = x1,
                          ymin = -Inf, ymax = Inf, fill = fc, alpha = 0.12),
        ggplot2::annotate("text", x = (x0 + x1) / 2, y = Inf,
                          label = pred, vjust = 1.5, size = 3.5, color = fc,
                          fontface = "bold")
      ))
    }
    layers
  }

  output$full_signal_plot <- shiny::renderPlot({
    df <- signal_data(); shiny::req(df)
    ps <- attr(df, "polya_start"); pe <- attr(df, "polya_end")
    dp <- if (nrow(df) > 20000) df[round(seq(1, nrow(df), length.out = 20000)), ] else df
    res_data <- selected_residue()
    ov <- tryCatch(
      .safe_nonA_overlay(read_nonA_data = res_data,
                         poly_tail_start = ps, poly_tail_end = pe, nonA_flank = 250),
      error = function(e) list()
    )
    p <- ggplot2::ggplot(dp, ggplot2::aes(x = position, y = signal, color = segment))
    for (l in ov) p <- p + l
    p + ggplot2::geom_line(size = 0.3) +
      ggplot2::scale_color_manual(values = c("Adapter" = "#089bcc", "Poly(A)" = "#f56042",
                                             "Transcript" = "#3a414d", "Other" = "#95a5a6"), name = "Region") +
      ggplot2::geom_vline(xintercept = ps, color = "#700f25", linetype = "dashed", size = 0.8) +
      ggplot2::geom_vline(xintercept = pe, color = "#0f3473", linetype = "dashed", size = 0.8) +
      ggplot2::labs(x = "Position (samples)", y = "Signal (raw)", title = paste("Read:", attr(df, "read_id"))) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold", size = 14))
  })

  output$polya_zoom_plot <- shiny::renderPlot({
    df <- signal_data(); shiny::req(df)
    ps <- attr(df, "polya_start"); pe <- attr(df, "polya_end")
    zs <- max(1, ps - 250); ze <- min(nrow(df), pe + 250)
    dz <- dplyr::filter(df, position >= zs, position <= ze)
    res_data <- selected_residue()
    ov <- tryCatch(
      .safe_nonA_overlay(read_nonA_data = res_data,
                         poly_tail_start = ps, poly_tail_end = pe, nonA_flank = 250),
      error = function(e) list()
    )
    p <- ggplot2::ggplot(dz, ggplot2::aes(x = position, y = signal, color = segment))
    for (l in ov) p <- p + l
    p + ggplot2::geom_line(size = 0.5) +
      ggplot2::scale_color_manual(values = c("Adapter" = "#089bcc", "Poly(A)" = "#f56042",
                                             "Transcript" = "#3a414d", "Other" = "#95a5a6"), name = "Region") +
      ggplot2::geom_vline(xintercept = ps, color = "#700f25", linetype = "dashed", size = 1) +
      ggplot2::geom_vline(xintercept = pe, color = "#0f3473", linetype = "dashed", size = 1) +
      ggplot2::labs(x = "Position (samples)", y = "Signal (raw)",
                    title = paste("Poly(A) region:", ps, "-", pe), subtitle = paste("Positions", zs, "to", ze)) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(color = "#666666", size = 11))
  })

  # ---- Plotly explorer ----
  output$interactive_signal_plot <- plotly::renderPlotly({
    df <- signal_data(); shiny::req(df)
    ps <- attr(df, "polya_start"); pe <- attr(df, "polya_end")
    dp <- if (nrow(df) > 50000) df[round(seq(1, nrow(df), length.out = 50000)), ] else df
    rr <- selected_residue(); ns <- list()
    if (!is.null(rr) && nrow(rr) > 0) {
      rc <- c("C" = "#3a424f", "G" = "#50a675", "U" = "#b0bdd4")
      rp <- ninetails:::.estimate_nonA_signal_pos(est_nonA_pos = rr$est_nonA_pos,
                                                  poly_tail_length = rr$polya_length[1], poly_tail_start = ps, poly_tail_end = pe)
      for (i in seq_len(nrow(rr))) {
        fc <- rc[as.character(rr$prediction[i])]; if (is.na(fc)) fc <- "#999999"
        ns <- c(ns, list(list(type = "rect", x0 = rp[i] - 250, x1 = rp[i] + 250,
                              y0 = 0, y1 = 1, yref = "paper", fillcolor = fc, opacity = 0.12,
                              line = list(width = 0))))
      }
    }
    plotly::plot_ly(dp, x = ~position, y = ~signal, color = ~segment,
                    colors = c("Adapter" = "#089bcc", "Poly(A)" = "#f56042",
                               "Transcript" = "#3a414d", "Other" = "#95a5a6"),
                    type = "scatter", mode = "lines", line = list(width = 1)) %>%
      plotly::layout(
        title = list(text = paste("Read:", attr(df, "read_id")),
                     font = list(family = "Open Sans", size = 16)),
        font = list(family = "Open Sans"),
        xaxis = list(title = "Position (samples)"), yaxis = list(title = "Signal (raw)"),
        legend = list(orientation = "h", y = -0.15), dragmode = "zoom",
        shapes = c(list(
          list(type = "line", x0 = ps, x1 = ps, y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#700f25", dash = "dash", width = 2)),
          list(type = "line", x0 = pe, x1 = pe, y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#0f3473", dash = "dash", width = 2))),
          ns)) %>%
      plotly::config(displayModeBar = TRUE, displaylogo = FALSE)
  })

} # server

shiny::shinyApp(ui = ui, server = server)
