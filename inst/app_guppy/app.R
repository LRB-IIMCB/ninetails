################################################################################
# Ninetails Analysis Dashboard — Guppy Legacy Pipeline
#
# Interactive Shiny application for exploring poly(A) tail composition
# analysis results from the Guppy legacy pipeline (check_tails_guppy).
# Uses fast5 files and Nanopolish poly(A) coordinates for signal
# visualization. No Python dependency required.
#
# Location: inst/app_guppy/app.R
# Launch: ninetails::launch_signal_browser_guppy()
#
# Tabs:
#   1. Classification — value boxes, Nanopolish QC, class counts, non-A
#      abundance in responsive flex layout
#   2. Residues — residue counts, per-sample rug density plots, summary DT
#   3. Poly(A) length — tail distribution with condition filter, palette,
#      summary statistics table
#   4. Signal Viewer — sub-tabs: Static Viewer / Dynamic Explorer using
#      fast5 functions with moves toggle
#   5. Download — configurable report generator
#   6. About — package info, citation, credits
#
# Static assets (inst/app_guppy/www/):
#   - logo.png, favicon.ico, IIMCB_logo.png
#
################################################################################

library(shiny)
library(ggplot2)
library(dplyr)
library(vroom)

options(scipen = 999) # disable scientific notation


################################################################################
# DATA FROM LAUNCHER
################################################################################

class_data    <- shiny::getShinyOption("ninetails.class_data")
residue_data  <- shiny::getShinyOption("ninetails.residue_data")
merged_data   <- shiny::getShinyOption("ninetails.merged_data")
signal_config <- shiny::getShinyOption("ninetails.signal_config", list())
basecall_group <- shiny::getShinyOption("ninetails.basecall_group", "Basecall_1D_000")

has_class   <- !is.null(class_data) && nrow(class_data) > 0
has_residue <- !is.null(residue_data) && nrow(residue_data) > 0
has_merged  <- !is.null(merged_data) && nrow(merged_data) > 0
has_signal  <- length(signal_config) > 0

available_groups <- character(0)
if (has_class) {
  candidates <- c("sample_name", "group")
  available_groups <- candidates[candidates %in% names(class_data)]
  if (length(available_groups) == 0) available_groups <- "contig"
}

transcript_col <- "contig"
if (has_class && "symbol" %in% names(class_data)) transcript_col <- "symbol"

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
    shiny::tags$img(src = icon_url, class = "vbox-watermark",
                    style = "opacity: 0.06;"),
    shiny::div(class = "vbox-content",
               shiny::div(class = "vbox-value", style = paste0("color: ", color, ";"),
                          value),
               shiny::div(class = "vbox-label", style = paste0("color: ", color, ";"),
                          label))
  )
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
      .header {
        background-color: #dce2e8; color: #333333; padding: 16px 20px;
        margin: -15px -15px 20px -15px;
        display: flex; align-items: center; gap: 14px;
      }
      .header h2 { margin: 0; font-weight: bold; font-size: 22px; }
      .header p { margin: 4px 0 0 0; color: #555555; font-size: 14px; }
      .vbox-row { display: flex; gap: 12px; margin-bottom: 15px; flex-wrap: wrap; }
      .ninetails-vbox {
        background: #ffffff; border-radius: 6px; padding: 18px 20px;
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
      .card {
        background: #ffffff; border-radius: 6px; padding: 15px;
        margin-bottom: 15px; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
      }
      .card h4 {
        margin-top: 0; color: #333333; font-weight: 600;
        border-bottom: 2px solid #2C7BB6; padding-bottom: 8px;
      }
      .filter-section {
        background: #f5f7f9; padding: 15px; border-radius: 6px; margin-bottom: 15px;
      }
      .placeholder-msg { text-align: center; padding: 60px 20px; color: #888888; }
      .placeholder-msg h4 { color: #aaaaaa; }
      .read-info {
        background: #f0f5fa; padding: 12px 15px; border-radius: 6px;
        border-left: 4px solid #2C7BB6;
      }
      .read-info p { margin: 5px 0; }
      .classification-info {
        background: #e4eef7; padding: 8px 12px; border-radius: 4px;
        margin-top: 10px; font-size: 12px; color: #1a4971;
      }
      .plots-flex { display: flex; flex-wrap: wrap; gap: 15px; }
      .plots-flex > .card { flex: 1 1 400px; max-width: 100%; }
      @media (min-width: 1200px) {
        .plots-flex > .card { max-width: calc(50% - 8px); }
      }
      .plot-desc {
        background: #f5f7f9; padding: 10px 14px; border-radius: 4px;
        font-size: 13px; color: #555555; line-height: 1.5;
        border-left: 3px solid #2C7BB6; margin-top: 8px;
      }
      .plot-desc b { color: #333333; }
      .about-section { max-width: 800px; margin: 0 auto; padding: 20px 0; }
      .about-section h3 { color: #2C7BB6; border-bottom: 2px solid #dce2e8; padding-bottom: 6px; }
      .about-section a { color: #2C7BB6; }
      .about-logo { display: block; margin: 0 auto 20px auto; }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        color: #2C7BB6; border-bottom: 2px solid #2C7BB6; font-weight: 600;
      }
      .nav-tabs > li > a { color: #555555; }
      .selectize-input, .selectize-dropdown { font-size: 12px; }
      .btn-primary { background-color: #2C7BB6; border-color: #2569a0; }
      .btn-primary:hover { background-color: #2569a0; }
      .legacy-badge {
        display: inline-block; background: #ff6600; color: #ffffff;
        font-size: 10px; font-weight: 700; padding: 2px 6px;
        border-radius: 3px; margin-left: 8px; vertical-align: middle;
      }
    ")) # tags$style
  ), # tags$head

  #  Header
  shiny::div(class = "header",
             htmltools::img(src = "logo.png", height = 45, width = 45, alt = "ninetails"),
             shiny::div(
               shiny::h2(htmltools::strong(paste0(
                 "Ninetails analysis dashboard ",
                 as.character(utils::packageVersion("ninetails"))
               )), shiny::tags$span(class = "legacy-badge", "GUPPY LEGACY")),
               shiny::p("Poly(A) tail composition analysis — Guppy pipeline (fast5)")
             )
  ),

  #  Main tabs
  shiny::tabsetPanel(id = "main_tabs", type = "tabs",

                     ############################################################################
                     # TAB 1: Classification
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
                                                                                                       "Decorated (A)" = "A"), selected = "N"),
                                                                        shiny::hr(),
                                                                        shiny::checkboxInput("show_class_desc", "Show descriptions", FALSE)
                                                             )
                                         ),
                                         shiny::mainPanel(width = 10,
                                                          shiny::div(class = "vbox-row",
                                                                     .value_box_ui(format(n_samples_val, big.mark = ","),
                                                                                   "Samples analyzed", "#2C7BB6"),
                                                                     .value_box_ui(format(n_transcripts_val, big.mark = ","),
                                                                                   "Transcripts found", "#2C7BB6"),
                                                                     .value_box_ui(format(n_reads_val, big.mark = ","),
                                                                                   "Total reads", "#555555"),
                                                                     .value_box_ui(format(n_blank_val, big.mark = ","),
                                                                                   "Blank reads", "#089bcc"),
                                                                     .value_box_ui(format(n_decorated_val, big.mark = ","),
                                                                                   "Decorated reads", "#ff6600")
                                                          ),
                                                          # Nanopolish QC (Guppy-specific)
                                                          shiny::div(class = "card",
                                                                     shiny::h4("Nanopolish QC"),
                                                                     shiny::plotOutput("nanopolish_qc_plot", height = "350px"),
                                                                     shiny::conditionalPanel("input.show_class_desc",
                                                                                             shiny::div(class = "plot-desc",
                                                                                                        shiny::HTML(paste0(
                                                                                                          "Distribution of QC tags assigned by Nanopolish polya. ",
                                                                                                          "Reads marked as <b>PASS</b> meet quality criteria for ",
                                                                                                          "ninetails analysis. Other categories: ",
                                                                                                          "<b>ADAPTER</b>, <b>NOREGION</b>, <b>SUFFCLIP</b>, etc."
                                                                                                        ))))
                                                          ),
                                                          shiny::div(class = "plots-flex",
                                                                     shiny::div(class = "card",
                                                                                shiny::h4("Read Classification"),
                                                                                shiny::plotOutput("class_plot", height = "450px"),
                                                                                shiny::conditionalPanel("input.show_class_desc",
                                                                                                        shiny::div(class = "plot-desc",
                                                                                                                   shiny::HTML(paste0(
                                                                                                                     "Read classification summary. ",
                                                                                                                     "<b>YAY</b> \u2014 move present, non-A detected; ",
                                                                                                                     "<b>MPU</b> \u2014 move present, only A detected; ",
                                                                                                                     "<b>MAU</b> \u2014 no move, unmodified; ",
                                                                                                                     "<b>IRL</b> \u2014 tail too short; ",
                                                                                                                     "<b>QCF</b> \u2014 Nanopolish QC failed; ",
                                                                                                                     "<b>NIN</b> \u2014 not included (pass_only)."
                                                                                                                   ))))
                                                                     ),
                                                                     if (has_residue) {
                                                                       shiny::div(class = "card",
                                                                                  shiny::h4("Non-A Abundance (reads with 1, 2, 3+ non-As)"),
                                                                                  shiny::plotOutput("nonA_abundance_plot", height = "400px"),
                                                                                  shiny::conditionalPanel("input.show_class_desc",
                                                                                                          shiny::div(class = "plot-desc",
                                                                                                                     shiny::HTML(paste0(
                                                                                                                       "Frequency of reads with one, two, or three or more ",
                                                                                                                       "non-adenosine residues per read, relative to decorated reads."
                                                                                                                     )))))
                                                                     }
                                                          ) # plots-flex
                                         ) # mainPanel
                                       ) # sidebarLayout
                                     } else {
                                       shiny::div(class = "placeholder-msg",
                                                  shiny::h4("No classification data loaded"))
                                     }
                     ), # Classification

                     ############################################################################
                     # TAB 2: Residues
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
                                                             )
                                         ),
                                         shiny::mainPanel(width = 10,
                                                          shiny::div(class = "card",
                                                                     shiny::h4("Residue Counts"),
                                                                     shiny::plotOutput("residue_plot", height = "400px"),
                                                                     shiny::conditionalPanel("input.show_residue_desc",
                                                                                             shiny::div(class = "plot-desc",
                                                                                                        shiny::HTML("Distribution of C, G, and U residue types.")))
                                                          ),
                                                          shiny::div(class = "card",
                                                                     shiny::h4("Non-A Position Distribution (max 1000 points per residue per sample)"),
                                                                     shiny::uiOutput("rug_plots_ui"),
                                                                     shiny::conditionalPanel("input.show_residue_desc",
                                                                                             shiny::div(class = "plot-desc",
                                                                                                        shiny::HTML(paste0(
                                                                                                          "Positional distribution of non-A residues along poly(A) tails. ",
                                                                                                          "Points subsampled to max 1,000 per residue type per sample."))))
                                                          ),
                                                          if (has_merged) {
                                                            shiny::div(class = "card",
                                                                       shiny::h4("Summary Table"),
                                                                       DT::DTOutput("residue_summary_table"))
                                                          }
                                         )
                                       )
                                     } else {
                                       shiny::div(class = "placeholder-msg",
                                                  shiny::h4("No residue data loaded"))
                                     }
                     ), # Residues

                     ############################################################################
                     # TAB 3: Poly(A) length
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
                                                                                            class = "btn-default btn-sm", icon = shiny::icon("sync"),
                                                                                            style = "width: 100%;"),
                                                                        shiny::hr(),
                                                                        shiny::checkboxInput("show_polya_desc", "Show descriptions", FALSE)
                                                             )
                                         ),
                                         shiny::mainPanel(width = 10,
                                                          shiny::div(class = "card",
                                                                     shiny::h4("Poly(A) Tail Length Distribution"),
                                                                     shiny::plotOutput("polya_dist_plot", height = "500px"),
                                                                     shiny::conditionalPanel("input.show_polya_desc",
                                                                                             shiny::div(class = "plot-desc",
                                                                                                        shiny::HTML("Density distribution of poly(A) tail lengths.")))
                                                          ),
                                                          shiny::div(class = "card",
                                                                     shiny::h4("Poly(A) Length Summary"),
                                                                     DT::DTOutput("polya_summary_table"))
                                         )
                                       )
                                     } else {
                                       shiny::div(class = "placeholder-msg",
                                                  shiny::h4("No data loaded"))
                                     }
                     ), # Poly(A) length

                     ############################################################################
                     # TAB 4: Signal Viewer (fast5)
                     ############################################################################
                     shiny::tabPanel("Signal Viewer", shiny::br(),
                                     shiny::fluidRow(
                                       shiny::column(3,
                                                     shiny::div(class = "card",
                                                                shiny::h4("Data"),
                                                                if (length(signal_config) > 0 && names(signal_config)[1] != "single") {
                                                                  shiny::div(class = "filter-section",
                                                                             shiny::selectInput("signal_sample", "Select sample",
                                                                                                choices = names(signal_config),
                                                                                                selected = names(signal_config)[1]))
                                                                } else if (length(signal_config) == 0) {
                                                                  shiny::div(class = "filter-section",
                                                                             shiny::textInput("signal_nanopolish", "Nanopolish polya file",
                                                                                              placeholder = "/path/to/nanopolish_output.tsv"),
                                                                             shiny::textInput("signal_seqsum", "Sequencing summary file",
                                                                                              placeholder = "/path/to/sequencing_summary.txt"),
                                                                             shiny::textInput("signal_workspace", "Fast5 directory",
                                                                                              placeholder = "/path/to/fast5/"),
                                                                             shiny::actionButton("load_signal_data", "Load data",
                                                                                                 class = "btn-primary btn-sm", icon = shiny::icon("folder-open")))
                                                                },
                                                                shiny::uiOutput("signal_data_status")
                                                     ),
                                                     shiny::div(class = "card",
                                                                shiny::h4("Filters"),
                                                                shiny::div(class = "filter-section",
                                                                           shiny::sliderInput("sig_polya_range", "Poly(A) tail length range",
                                                                                              min = 0, max = 500, value = c(10, 500), step = 1),
                                                                           shiny::selectInput("comments_filter", "Decoration status",
                                                                                              choices = c("All", "YAY", "MPU", "MAU", "QCF", "NIN", "IRL"),
                                                                                              selected = "All"),
                                                                           shiny::checkboxInput("sig_show_moves", "Show moves", FALSE),
                                                                           shiny::checkboxInput("sig_rescale", "Rescale to pA", FALSE)
                                                                ),
                                                                shiny::hr(),
                                                                shiny::uiOutput("read_selection_ui"),
                                                                shiny::uiOutput("filter_summary")
                                                     )
                                       ), # column 3
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
                                                                                                                           shiny::h4("Poly(A) Region"),
                                                                                                                           shiny::plotOutput("polya_zoom_plot", height = "350px"))
                                                                                        ),
                                                                                        shiny::conditionalPanel("!output.signal_loaded",
                                                                                                                shiny::div(class = "card",
                                                                                                                           shiny::div(class = "placeholder-msg",
                                                                                                                                      shiny::h4("No signal loaded"),
                                                                                                                                      shiny::p("Load data and select a read."))))
                                                                        ), # Static Viewer
                                                                        shiny::tabPanel("Dynamic Explorer", shiny::br(),
                                                                                        shiny::conditionalPanel("output.signal_loaded",
                                                                                                                shiny::div(class = "card",
                                                                                                                           shiny::h4("Selected Read Information"),
                                                                                                                           shiny::uiOutput("read_info_explorer")),
                                                                                                                shiny::div(class = "card",
                                                                                                                           shiny::h4("Interactive Signal Explorer"),
                                                                                                                           shiny::p(style = "color: #666666; font-size: 13px;",
                                                                                                                                    "Click and drag to zoom. Double-click to reset."),
                                                                                                                           plotly::plotlyOutput("interactive_signal_plot", height = "500px"))
                                                                                        ),
                                                                                        shiny::conditionalPanel("!output.signal_loaded",
                                                                                                                shiny::div(class = "card",
                                                                                                                           shiny::div(class = "placeholder-msg",
                                                                                                                                      shiny::h4("No signal loaded"),
                                                                                                                                      shiny::p("Load data and select a read."))))
                                                                        ) # Dynamic Explorer
                                                     ) # tabsetPanel signal
                                       ) # column 9
                                     ) # fluidRow
                     ), # Signal Viewer

                     ############################################################################
                     # TAB 5: Download
                     ############################################################################
                     shiny::tabPanel("Download", shiny::br(),
                                     if (has_class) {
                                       shiny::div(style = "max-width: 800px; margin: 0 auto;",
                                                  shiny::div(class = "card",
                                                             shiny::h4("Report Configuration"),
                                                             shiny::checkboxInput("rpt_classification", "Classification plot", TRUE),
                                                             shiny::checkboxInput("rpt_qc", "Nanopolish QC plot", TRUE),
                                                             if (has_residue) shiny::checkboxInput("rpt_abundance",
                                                                                                   "Non-A abundance", TRUE),
                                                             if (has_residue) shiny::checkboxInput("rpt_residue",
                                                                                                   "Residue frequency", TRUE),
                                                             if (has_residue) shiny::checkboxInput("rpt_rug",
                                                                                                   "Rug density plots", TRUE),
                                                             shiny::checkboxInput("rpt_polya", "Poly(A) length distribution", TRUE),
                                                             shiny::hr(),
                                                             shiny::selectInput("rpt_center", "Central tendency",
                                                                                choices = c("mean", "median", "mode", "none"),
                                                                                selected = "none", width = "250px"),
                                                             shiny::selectInput("rpt_palette", "Color palette",
                                                                                choices = names(palettes), selected = "unova", width = "250px"),
                                                             shiny::hr(),
                                                             shiny::downloadButton("download_report", "Download report (.html)",
                                                                                   class = "btn-primary", style = "width: 100%; margin-top: 10px;")
                                                  )
                                       )
                                     } else {
                                       shiny::div(class = "placeholder-msg",
                                                  shiny::h4("No data loaded"))
                                     }
                     ), # Download

                     ############################################################################
                     # TAB 6: About
                     ############################################################################
                     shiny::tabPanel("About", shiny::br(),
                                     shiny::div(class = "about-section",
                                                htmltools::img(src = "logo.png", height = 120, width = 120,
                                                               class = "about-logo"),
                                                shiny::h3("Ninetails"),
                                                shiny::p(htmltools::strong(paste0(
                                                  "Version ", as.character(utils::packageVersion("ninetails"))
                                                )), shiny::tags$span(class = "legacy-badge", "GUPPY LEGACY")),
                                                shiny::p("An R package for finding non-adenosine residues in poly(A)",
                                                         "tails using convolutional neural networks on raw nanopore signal."),
                                                shiny::p(htmltools::HTML(paste0(
                                                  "<b>Note:</b> This dashboard is for the Guppy legacy pipeline ",
                                                  "(fast5 + Nanopolish). For Dorado DRS data, use ",
                                                  "<code>launch_signal_browser()</code> instead."
                                                ))),
                                                shiny::h3("Citation"),
                                                shiny::p(htmltools::HTML(paste0(
                                                  "Gumi\u0144ska, N., Matylla-Kuli\u0144ska, K., Krawczyk, P.S. ",
                                                  "<i>et al.</i> Direct profiling of non-adenosines in poly(A) tails. ",
                                                  "<i>Nat Commun</i> <b>16</b>, 2664 (2025). ",
                                                  "<a href='https://doi.org/10.1038/s41467-025-57787-6' target='_blank'>",
                                                  "https://doi.org/10.1038/s41467-025-57787-6</a>"
                                                ))),
                                                shiny::h3("Links"),
                                                shiny::tags$ul(
                                                  shiny::tags$li(shiny::tags$a(href = "https://github.com/LRB-IIMCB/ninetails",
                                                                               target = "_blank", "GitHub repository")),
                                                  shiny::tags$li(shiny::tags$a(href = "https://github.com/LRB-IIMCB/ninetails/wiki",
                                                                               target = "_blank", "Documentation (Wiki)")),
                                                  shiny::tags$li(shiny::tags$a(href = "https://LRB-IIMCB.github.io/ninetails/",
                                                                               target = "_blank", "Package website"))
                                                ),
                                                shiny::h3("Laboratory"),
                                                htmltools::img(src = "IIMCB_logo.png", height = 90, class = "about-logo"),
                                                shiny::p(htmltools::HTML(paste0(
                                                  "<a href='https://www.iimcb.gov.pl/en/research/",
                                                  "41-laboratory-of-rna-biology-era-chairs-group' target='_blank'>",
                                                  "Laboratory of RNA Biology</a> (ERA Chairs Group), IIMCB Warsaw"
                                                ))),
                                                shiny::h3("Developer"),
                                                shiny::p(htmltools::HTML(
                                                  "Natalia Gumi\u0144ska (<a href='mailto:nguminska@iimcb.gov.pl'>nguminska@iimcb.gov.pl</a>)"
                                                )),
                                                shiny::hr(),
                                                shiny::p(style = "color: #999999; font-size: 12px; text-align: center;",
                                                         paste0("Dashboard built with Shiny ",
                                                                utils::packageVersion("shiny"),
                                                                " | R ", R.version$major, ".", R.version$minor))
                                     )
                     ) # About

  ) # tabsetPanel
) # fluidPage


################################################################################
# SERVER
################################################################################

server <- function(input, output, session) {

  ##############################################################################
  # SHARED
  ##############################################################################

  if (has_class) {
    tc <- sort(unique(as.character(class_data[[transcript_col]])))
    shiny::updateSelectizeInput(session, "class_contig_filter",
                                choices = c("All" = "", tc), selected = "", server = TRUE)
    shiny::updateSelectizeInput(session, "polya_contig_filter",
                                choices = c("All" = "", tc), selected = "", server = TRUE)
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

    # Nanopolish QC plot (Guppy-specific)
    output$nanopolish_qc_plot <- shiny::renderPlot({
      cd <- class_data_filtered(); shiny::req(nrow(cd) > 0)
      shiny::req(input$class_grouping)
      gf <- if (input$class_grouping %in% names(cd)) input$class_grouping else NA
      tryCatch({
        qc_data <- ninetails::nanopolish_qc(class_data = cd,
                                            grouping_factor = gf)
        ninetails::plot_nanopolish_qc(qc_data,
                                      frequency = isTRUE(input$class_frequency))
      }, error = function(e) {
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = paste("Nanopolish QC error:", e$message),
                            size = 4, color = "#cc0000") +
          ggplot2::theme_void()
      })
    })

    output$class_plot <- shiny::renderPlot({
      cd <- class_data_filtered(); shiny::req(nrow(cd) > 0)
      shiny::req(input$class_grouping, input$class_plot_type,
                 !is.null(input$class_frequency))
      gf <- if (input$class_grouping %in% names(cd)) input$class_grouping else NA
      p <- ninetails::plot_class_counts(class_data = cd,
                                        grouping_factor = gf, frequency = input$class_frequency,
                                        type = input$class_plot_type)
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
                                label = paste("Error:", e$message),
                                size = 4, color = "#cc0000") +
              ggplot2::theme_void()
          }
        )
      })
    }
  }

  ##############################################################################
  # RESIDUES TAB
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

    # Per-sample rug plots (same pattern as Dorado app)
    .safe_id <- function(x) gsub("[^A-Za-z0-9]", "_", x)

    output$rug_plots_ui <- shiny::renderUI({
      rd <- residue_data_filtered(); shiny::req(nrow(rd) > 0)
      sample_names <- if ("sample_name" %in% names(rd)) {
        sort(unique(as.character(rd$sample_name)))
      } else { "all" }

      ui_list <- lapply(sample_names, function(sn) {
        sid <- .safe_id(sn)
        shiny::tagList(
          shiny::tags$h5(style = "color: #2C7BB6; font-weight: 600; margin-top: 15px;",
                         if (sn == "all") "All data" else sn),
          shiny::fluidRow(
            shiny::column(4, shiny::plotOutput(paste0("rug_", sid, "_C"), height = "350px")),
            shiny::column(4, shiny::plotOutput(paste0("rug_", sid, "_G"), height = "350px")),
            shiny::column(4, shiny::plotOutput(paste0("rug_", sid, "_U"), height = "350px")))
        )
      })

      for (sn in sample_names) {
        local({
          local_sn <- sn; sid <- .safe_id(local_sn)
          for (base_type in c("C", "G", "U")) {
            local({
              local_base <- base_type
              output[[paste0("rug_", sid, "_", local_base)]] <- shiny::renderPlot({
                rd <- residue_data_filtered(); shiny::req(nrow(rd) > 0, input$rug_max_length)
                if (local_sn != "all" && "sample_name" %in% names(rd))
                  rd <- rd[rd$sample_name == local_sn, , drop = FALSE]
                if (!local_base %in% rd$prediction || nrow(rd) == 0) {
                  return(ggplot2::ggplot() +
                           ggplot2::annotate("text", x = 0.5, y = 0.5,
                                             label = paste("No", local_base, "residues"),
                                             size = 4, color = "#999999") +
                           ggplot2::theme_void() + ggplot2::ggtitle(local_base))
                }
                rd_base <- rd[rd$prediction == local_base, , drop = FALSE]
                n_avail <- nrow(rd_base)
                if (n_avail > 1000) {
                  set.seed(42)
                  keep_idx <- sample.int(n_avail, 1000, replace = FALSE)
                  rd <- rbind(rd[rd$prediction != local_base, , drop = FALSE],
                              rd_base[keep_idx, ])
                }
                tryCatch(
                  ninetails::plot_rug_density(residue_data = rd,
                                              base = local_base, max_length = input$rug_max_length),
                  error = function(e) {
                    ggplot2::ggplot() +
                      ggplot2::annotate("text", x = 0.5, y = 0.5,
                                        label = paste(local_base, "error:", e$message),
                                        size = 3.5, color = "#cc0000") +
                      ggplot2::theme_void()
                  })
              })
            })
          }
        })
      }
      do.call(shiny::tagList, ui_list)
    })

    if (has_merged) {
      output$residue_summary_table <- DT::renderDT({
        shiny::req(input$residue_grouping)
        md <- .filter_by_contig(merged_data, input$residue_contig_filter)
        shiny::req(nrow(md) > 0)
        gf <- if (input$residue_grouping %in% names(md)) input$residue_grouping else "contig"
        tryCatch({
          tbl <- ninetails::summarize_nonA(merged_nonA_tables = md,
                                           summary_factors = gf, transcript_id_column = NULL)
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
  # POLY(A) LENGTH TAB
  ##############################################################################

  if (has_class) {

    output$polya_condition_filter_ui <- shiny::renderUI({
      shiny::req(input$polya_grouping)
      base <- if (has_merged) merged_data else class_data
      if (!input$polya_grouping %in% names(base)) return(NULL)
      conditions <- sort(unique(as.character(base[[input$polya_grouping]])))
      shiny::selectizeInput("polya_condition_filter", "Select condition(s)",
                            choices = conditions, selected = conditions, multiple = TRUE,
                            options = list(placeholder = "All conditions",
                                           plugins = list("remove_button")))
    })

    shiny::observeEvent(input$polya_reset, {
      shiny::updateSelectInput(session, "polya_grouping", selected = available_groups[1])
      shiny::updateSelectizeInput(session, "polya_contig_filter", selected = "")
      shiny::updateSelectInput(session, "polya_center", selected = "none")
      shiny::updateSliderInput(session, "polya_max_length", value = 200)
      shiny::updateCheckboxInput(session, "polya_ndensity", value = TRUE)
      shiny::updateSelectInput(session, "polya_palette", selected = "unova")
    })

    polya_data_filtered <- shiny::reactive({
      base <- if (has_merged) merged_data else class_data
      base <- .filter_by_contig(base, input$polya_contig_filter)
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
      p + ggplot2::scale_color_manual(values = palettes[[input$polya_palette]])
    })

    output$polya_summary_table <- DT::renderDT({
      pd <- polya_data_filtered(); shiny::req(nrow(pd) > 0)
      shiny::req(input$polya_grouping)
      gf <- if (input$polya_grouping %in% names(pd)) input$polya_grouping else NULL
      len_col <- if ("polya_length" %in% names(pd)) "polya_length" else NULL
      shiny::req(len_col)
      if (!is.null(gf)) {
        tbl <- pd %>%
          dplyr::group_by(!!rlang::sym(gf)) %>%
          dplyr::summarise(
            n_reads = dplyr::n(),
            mean_length = mean(!!rlang::sym(len_col), na.rm = TRUE),
            median_length = stats::median(!!rlang::sym(len_col), na.rm = TRUE),
            sd = stats::sd(!!rlang::sym(len_col), na.rm = TRUE),
            sem = stats::sd(!!rlang::sym(len_col), na.rm = TRUE) / sqrt(dplyr::n()),
            .groups = "drop")
      } else {
        tbl <- pd %>%
          dplyr::summarise(
            n_reads = dplyr::n(),
            mean_length = mean(!!rlang::sym(len_col), na.rm = TRUE),
            median_length = stats::median(!!rlang::sym(len_col), na.rm = TRUE),
            sd = stats::sd(!!rlang::sym(len_col), na.rm = TRUE),
            sem = stats::sd(!!rlang::sym(len_col), na.rm = TRUE) / sqrt(dplyr::n()))
      }
      tbl <- dplyr::mutate(tbl, dplyr::across(
        tidyselect::vars_select_helpers$where(is.numeric), ~ round(., 2)))
      DT::datatable(tbl,
                    options = list(dom = "t", paging = FALSE, searching = FALSE, ordering = FALSE),
                    rownames = FALSE,
                    colnames = c(if (!is.null(gf)) gf else NULL,
                                 "Reads", "Mean length", "Median length", "SD", "SEM"))
    })
  }

  ##############################################################################
  # SIGNAL VIEWER TAB (fast5-based)
  ##############################################################################

  loaded_nanopolish <- shiny::reactiveVal(NULL)
  loaded_seqsum     <- shiny::reactiveVal(NULL)
  active_workspace  <- shiny::reactiveVal(NULL)
  signal_error      <- shiny::reactiveVal(NULL)

  # Helper: ensure nanopolish data has both readname and read_id columns
  # (plot_squiggle_fast5 uses read_id, plot_tail_range_fast5 uses readname)
  .normalize_nanopolish <- function(np) {
    if ("readname" %in% names(np) && !"read_id" %in% names(np))
      np$read_id <- np$readname
    if ("read_id" %in% names(np) && !"readname" %in% names(np))
      np$readname <- np$read_id
    np
  }

  #  Data loading
  if (length(signal_config) > 0 && names(signal_config)[1] != "single") {
    shiny::observeEvent(input$signal_sample, {
      s <- signal_config[[input$signal_sample]]
      if (is.null(s)) { signal_error("Sample not found."); return() }
      tryCatch({
        np <- vroom::vroom(s$nanopolish_path, show_col_types = FALSE)
        np <- .normalize_nanopolish(np)
        ss <- vroom::vroom(s$sequencing_summary_path, show_col_types = FALSE)
        loaded_nanopolish(np); loaded_seqsum(ss)
        active_workspace(s$workspace); signal_error(NULL)
      }, error = function(e) {
        signal_error(paste("Error:", e$message))
        loaded_nanopolish(NULL); loaded_seqsum(NULL)
      })
    }, ignoreNULL = FALSE)
  } else if (length(signal_config) > 0 && names(signal_config)[1] == "single") {
    shiny::observe({
      s <- signal_config[["single"]]
      tryCatch({
        np <- vroom::vroom(s$nanopolish_path, show_col_types = FALSE)
        np <- .normalize_nanopolish(np)
        ss <- vroom::vroom(s$sequencing_summary_path, show_col_types = FALSE)
        loaded_nanopolish(np); loaded_seqsum(ss)
        active_workspace(s$workspace); signal_error(NULL)
      }, error = function(e) {
        signal_error(paste("Error:", e$message))
        loaded_nanopolish(NULL); loaded_seqsum(NULL)
      })
    })
  } else {
    shiny::observeEvent(input$load_signal_data, {
      np_f <- trimws(input$signal_nanopolish)
      ss_f <- trimws(input$signal_seqsum)
      ws_f <- trimws(input$signal_workspace)
      if (!nzchar(np_f) || !file.exists(np_f)) {
        signal_error("Nanopolish file not found."); return() }
      if (!nzchar(ss_f) || !file.exists(ss_f)) {
        signal_error("Sequencing summary not found."); return() }
      if (!nzchar(ws_f) || !dir.exists(ws_f)) {
        signal_error("Fast5 directory not found."); return() }
      tryCatch({
        np <- vroom::vroom(np_f, show_col_types = FALSE)
        np <- .normalize_nanopolish(np)
        ss <- vroom::vroom(ss_f, show_col_types = FALSE)
        loaded_nanopolish(np); loaded_seqsum(ss)
        active_workspace(ws_f); signal_error(NULL)
      }, error = function(e) {
        signal_error(paste("Error:", e$message))
        loaded_nanopolish(NULL); loaded_seqsum(NULL)
      })
    })
  }

  output$signal_data_status <- shiny::renderUI({
    err <- signal_error(); np <- loaded_nanopolish()
    if (!is.null(err))
      shiny::tags$div(style = "background: #f8d7da; color: #721c24; padding: 10px; border-radius: 4px; margin-top: 10px;",
                      shiny::tags$strong("Error: "), err)
    else if (!is.null(np))
      shiny::tags$div(style = "background: #d4edda; color: #155724; padding: 10px; border-radius: 4px; margin-top: 10px;",
                      shiny::tags$strong("Loaded: "), format(nrow(np), big.mark = ","), " reads")
  })

  #  Filtered reads
  sig_filtered_data <- shiny::reactive({
    np <- loaded_nanopolish(); shiny::req(np)
    # Normalize read ID column
    if (!"readname" %in% names(np) && "read_id" %in% names(np))
      np <- dplyr::rename(np, readname = read_id)
    # Enrich with class/comments from class_data if available
    if (has_class && !"comments" %in% names(np)) {
      cd_join <- class_data[, c("readname", "class", "comments"), drop = FALSE]
      cd_join <- dplyr::distinct(cd_join, readname, .keep_all = TRUE)
      np <- dplyr::left_join(np, cd_join, by = "readname")
    }
    # Apply filters
    if (!is.null(input$sig_polya_range) && "polya_length" %in% names(np))
      np <- dplyr::filter(np, polya_length >= input$sig_polya_range[1],
                          polya_length <= input$sig_polya_range[2])
    if (!is.null(input$comments_filter) && input$comments_filter != "All" &&
        "comments" %in% names(np))
      np <- dplyr::filter(np, comments == input$comments_filter)
    # Keep only PASS reads with valid coords
    if ("qc_tag" %in% names(np))
      np <- dplyr::filter(np, qc_tag == "PASS")
    np
  })

  #  Read selection
  output$read_selection_ui <- shiny::renderUI({
    data <- sig_filtered_data()
    if (is.null(data) || nrow(data) == 0) return(shiny::helpText("No reads match filters."))
    shiny::tagList(
      shiny::selectizeInput("selected_read_id", "Select Read", choices = NULL,
                            options = list(placeholder = "Type to search...", maxOptions = 100)),
      shiny::uiOutput("read_nav_buttons"))
  })
  shiny::observe({
    data <- sig_filtered_data()
    if (!is.null(data) && nrow(data) > 0)
      shiny::updateSelectizeInput(session, "selected_read_id",
                                  choices = as.character(data$readname),
                                  selected = as.character(data$readname[1]), server = TRUE)
  })

  output$read_nav_buttons <- shiny::renderUI({
    data <- sig_filtered_data()
    if (is.null(data) || nrow(data) == 0) return(NULL)
    ids <- as.character(data$readname)
    idx <- match(input$selected_read_id, ids); if (is.na(idx)) idx <- 1
    shiny::div(style = "display: flex; gap: 8px; margin-top: 10px;",
               shiny::tags$button(id = "prev_read", type = "button",
                                  class = paste("btn btn-default action-button", if (idx <= 1) "disabled"),
                                  disabled = if (idx <= 1) "disabled", shiny::HTML("&larr; Previous")),
               shiny::tags$button(id = "next_read", type = "button",
                                  class = paste("btn btn-default action-button", if (idx >= length(ids)) "disabled"),
                                  disabled = if (idx >= length(ids)) "disabled", shiny::HTML("Next &rarr;")),
               shiny::tags$span(style = "align-self: center; margin-left: 8px; color: #666666; font-size: 12px;",
                                paste0("(", idx, " / ", length(ids), ")")))
  })
  shiny::observeEvent(input$prev_read, {
    data <- sig_filtered_data(); shiny::req(data, nrow(data) > 0)
    ids <- as.character(data$readname); idx <- match(input$selected_read_id, ids)
    if (!is.na(idx) && idx > 1) shiny::updateSelectizeInput(session, "selected_read_id", selected = ids[idx - 1])
  }, ignoreInit = TRUE)
  shiny::observeEvent(input$next_read, {
    data <- sig_filtered_data(); shiny::req(data, nrow(data) > 0)
    ids <- as.character(data$readname); idx <- match(input$selected_read_id, ids)
    if (!is.na(idx) && idx < length(ids)) shiny::updateSelectizeInput(session, "selected_read_id", selected = ids[idx + 1])
  }, ignoreInit = TRUE)

  output$filter_summary <- shiny::renderUI({
    data <- sig_filtered_data(); total <- loaded_nanopolish()
    if (is.null(data) || is.null(total)) return(NULL)
    shiny::tags$p(style = "color: #666666; font-size: 13px; margin-top: 10px;",
                  paste0("Showing ", nrow(data), " of ", nrow(total), " reads"))
  })

  #  Signal plots
  output$signal_loaded <- shiny::reactive({
    !is.null(loaded_nanopolish()) && !is.null(input$selected_read_id)
  })
  shiny::outputOptions(output, "signal_loaded", suspendWhenHidden = FALSE)

  .render_read_info <- function() {
    np <- loaded_nanopolish(); shiny::req(np, input$selected_read_id)
    ri <- np[np$readname == input$selected_read_id, , drop = FALSE]
    if (nrow(ri) == 0) return(NULL); ri <- ri[1, ]
    cc <- if ("comments" %in% names(ri)) ri$comments else NULL
    ct <- if (!is.null(cc) && !is.na(cc)) {
      switch(as.character(cc),
             "YAY" = "Decorated: Move present, non-A detected",
             "MPU" = "Blank: Move present, only A detected",
             "MAU" = "Blank: No move, unmodified",
             "IRL" = "Unclassified: Tail too short",
             "QCF" = "Unclassified: Nanopolish QC failed",
             "NIN" = "Unclassified: Not included (pass_only)", NULL)
    } else { NULL }
    # Non-A residue info
    res_data <- selected_residue()
    res_ui <- NULL
    if (!is.null(res_data) && nrow(res_data) > 0) {
      rs <- paste(vapply(seq_len(nrow(res_data)), function(i) {
        paste0(res_data$prediction[i], " (pos ", res_data$est_nonA_pos[i], " from 3')")
      }, character(1)), collapse = "; ")
      res_ui <- shiny::tags$div(
        style = "background: #e8f5e9; padding: 8px 12px; border-radius: 4px; margin-top: 8px; font-size: 12px; color: #2e7d32; border-left: 3px solid #50a675;",
        shiny::tags$strong("Non-A residues: "), rs)
    }
    shiny::tags$div(class = "read-info",
                    shiny::tags$p(shiny::tags$strong("Read ID: "), shiny::tags$code(input$selected_read_id)),
                    if ("polya_length" %in% names(ri))
                      shiny::tags$p(shiny::tags$strong("Poly(A) length: "), ri$polya_length, " nt"),
                    if ("contig" %in% names(ri))
                      shiny::tags$p(shiny::tags$strong("Transcript: "), ri$contig),
                    if (!is.null(cc)) shiny::tags$p(shiny::tags$strong("Status: "), cc),
                    if (!is.null(ct)) shiny::tags$div(class = "classification-info",
                                                      shiny::tags$strong("Classification: "), ct),
                    res_ui)
  }
  output$read_info <- shiny::renderUI({ .render_read_info() })
  output$read_info_explorer <- shiny::renderUI({ .render_read_info() })

  #  Selected residue data for overlay
  selected_residue <- shiny::reactive({
    shiny::req(input$selected_read_id)
    if (!has_residue) return(NULL)
    rd <- residue_data
    if (!"readname" %in% names(rd) && "read_id" %in% names(rd))
      rd <- dplyr::rename(rd, readname = read_id)
    if (!"readname" %in% names(rd)) return(NULL)
    rr <- rd[rd$readname == input$selected_read_id, , drop = FALSE]
    if (nrow(rr) == 0) NULL else rr
  })

  #  Helper: get poly(A) boundaries from nanopolish for selected read
  .get_polya_bounds <- function() {
    np <- loaded_nanopolish(); shiny::req(np, input$selected_read_id)
    id_col <- if ("read_id" %in% names(np)) "read_id" else "readname"
    ri <- np[np[[id_col]] == input$selected_read_id, , drop = FALSE]
    if (nrow(ri) == 0) return(NULL)
    ri <- ri[1, ]
    ps <- if ("polya_start" %in% names(ri)) ri$polya_start else NULL
    pe <- if ("transcript_start" %in% names(ri)) ri$transcript_start else NULL
    pl <- if ("polya_length" %in% names(ri)) ri$polya_length else NULL
    if (is.null(ps) || is.null(pe) || is.na(ps) || is.na(pe)) return(NULL)
    list(polya_start = ps, polya_end = pe, polya_length = pl)
  }

  #  Non-A overlay builder (same approach as Dorado app)
  .safe_nonA_overlay <- function(read_nonA_data, poly_tail_start,
                                 poly_tail_end, polya_length,
                                 nonA_flank = 250) {
    layers <- list()
    if (is.null(read_nonA_data) || nrow(read_nonA_data) == 0) return(layers)
    if (is.null(polya_length) || is.na(polya_length) || polya_length <= 0) return(layers)

    rc <- c("C" = "#3a424f", "G" = "#50a675", "U" = "#b0bdd4")

    for (i in seq_len(nrow(read_nonA_data))) {
      pred <- as.character(read_nonA_data$prediction[i])
      pos  <- read_nonA_data$est_nonA_pos[i]
      if (is.na(pos) || is.na(pred)) next

      fc <- rc[pred]
      if (is.na(fc)) next

      raw_pos <- poly_tail_start +
        (polya_length - pos) * (poly_tail_end - poly_tail_start) / polya_length
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
    shiny::req(input$selected_read_id, loaded_nanopolish(), loaded_seqsum(),
               active_workspace())
    shiny::withProgress(message = "Extracting signal...", value = 0.5, {
      tryCatch({
        p <- ninetails::plot_squiggle_fast5(
          readname = input$selected_read_id,
          nanopolish = loaded_nanopolish(),
          sequencing_summary = loaded_seqsum(),
          workspace = active_workspace(),
          basecall_group = basecall_group,
          moves = isTRUE(input$sig_show_moves),
          rescale = isTRUE(input$sig_rescale))

        # Add non-A overlay (only for raw signal, rescale=FALSE)
        if (!isTRUE(input$sig_rescale)) {
          bounds <- .get_polya_bounds()
          res_data <- selected_residue()
          if (!is.null(bounds) && !is.null(res_data)) {
            ov <- tryCatch(
              .safe_nonA_overlay(
                read_nonA_data = res_data,
                poly_tail_start = bounds$polya_start,
                poly_tail_end = bounds$polya_end,
                polya_length = bounds$polya_length,
                nonA_flank = 250),
              error = function(e) list())
            for (l in ov) p <- p + l
          }
        }
        p
      },
      error = function(e) {
        err_msg <- conditionMessage(e)
        if (!nzchar(err_msg)) err_msg <- as.character(e)
        cat("[Signal Viewer] Full signal ERROR:", err_msg, "\n")
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = paste("Signal error:", err_msg),
                            size = 3.5, color = "#cc0000", hjust = 0.5) +
          ggplot2::theme_void()
      })
    })
  })

  output$polya_zoom_plot <- shiny::renderPlot({
    shiny::req(input$selected_read_id, loaded_nanopolish(), loaded_seqsum(),
               active_workspace())
    shiny::withProgress(message = "Extracting tail region...", value = 0.5, {
      tryCatch({
        p <- ninetails::plot_tail_range_fast5(
          readname = input$selected_read_id,
          nanopolish = loaded_nanopolish(),
          sequencing_summary = loaded_seqsum(),
          workspace = active_workspace(),
          basecall_group = basecall_group,
          moves = isTRUE(input$sig_show_moves),
          rescale = isTRUE(input$sig_rescale))

        # Add non-A overlay (only for raw signal, rescale=FALSE)
        if (!isTRUE(input$sig_rescale)) {
          bounds <- .get_polya_bounds()
          res_data <- selected_residue()
          if (!is.null(bounds) && !is.null(res_data)) {
            ov <- tryCatch(
              .safe_nonA_overlay(
                read_nonA_data = res_data,
                poly_tail_start = bounds$polya_start,
                poly_tail_end = bounds$polya_end,
                polya_length = bounds$polya_length,
                nonA_flank = 250),
              error = function(e) list())
            for (l in ov) p <- p + l
          }
        }
        p
      },
      error = function(e) {
        err_msg <- conditionMessage(e)
        if (!nzchar(err_msg)) err_msg <- as.character(e)
        cat("[Signal Viewer] Tail range ERROR:", err_msg, "\n")
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = paste("Tail range error:", err_msg),
                            size = 3.5, color = "#cc0000", hjust = 0.5) +
          ggplot2::theme_void()
      })
    })
  })

  # Plotly explorer
  output$interactive_signal_plot <- plotly::renderPlotly({
    shiny::req(input$selected_read_id, loaded_nanopolish(), loaded_seqsum(),
               active_workspace())
    tryCatch({
      p <- ninetails::plot_squiggle_fast5(
        readname = input$selected_read_id,
        nanopolish = loaded_nanopolish(),
        sequencing_summary = loaded_seqsum(),
        workspace = active_workspace(),
        basecall_group = basecall_group,
        moves = FALSE,
        rescale = isTRUE(input$sig_rescale))

      # Build plotly shapes for non-A overlay
      ns <- list()
      if (!isTRUE(input$sig_rescale)) {
        bounds <- .get_polya_bounds()
        res_data <- selected_residue()
        if (!is.null(bounds) && !is.null(res_data) && nrow(res_data) > 0 &&
            !is.null(bounds$polya_length) && bounds$polya_length > 0) {
          rc <- c("C" = "#3a424f", "G" = "#50a675", "U" = "#b0bdd4")
          for (i in seq_len(nrow(res_data))) {
            pos <- res_data$est_nonA_pos[i]
            pred <- as.character(res_data$prediction[i])
            if (is.na(pos) || is.na(pred)) next
            fc <- rc[pred]; if (is.na(fc)) next
            raw_pos <- bounds$polya_start +
              (bounds$polya_length - pos) *
              (bounds$polya_end - bounds$polya_start) / bounds$polya_length
            ns <- c(ns, list(list(
              type = "rect", x0 = raw_pos - 250, x1 = raw_pos + 250,
              y0 = 0, y1 = 1, yref = "paper",
              fillcolor = fc, opacity = 0.12,
              line = list(width = 0))))
          }
        }
      }

      plotly::ggplotly(p) %>%
        plotly::layout(dragmode = "zoom", shapes = ns) %>%
        plotly::config(displayModeBar = TRUE, displaylogo = FALSE)
    }, error = function(e) {
      err_msg <- conditionMessage(e)
      if (!nzchar(err_msg)) err_msg <- as.character(e)
      cat("[Signal Viewer] Plotly ERROR:", err_msg, "\n")
      plotly::plot_ly(type = "scatter", mode = "markers") %>%
        plotly::layout(title = list(
          text = paste("Error:", err_msg),
          font = list(size = 12, color = "#cc0000")))
    })
  })

  ##############################################################################
  # DOWNLOAD TAB
  ##############################################################################

  if (has_class) {
    output$download_report <- shiny::downloadHandler(
      filename = function() {
        paste0("ninetails_guppy_report_", format(Sys.Date(), "%Y%m%d"), ".html")
      },
      content = function(file) {
        shiny::withProgress(message = "Generating report...", value = 0, {

          tmp_dir <- tempdir()
          gf <- if ("sample_name" %in% names(class_data)) "sample_name" else NA

          .plot_to_base64 <- function(p, w = 10, h = 6) {
            png_file <- tempfile(tmpdir = tmp_dir, fileext = ".png")
            tryCatch({
              ggplot2::ggsave(png_file, p, width = w, height = h, dpi = 150)
              b64 <- base64enc::dataURI(file = png_file, mime = "image/png")
              unlink(png_file); b64
            }, error = function(e) { unlink(png_file); NULL })
          }

          rpt_center_val <- if (!is.null(input$rpt_center) &&
                                input$rpt_center != "none") input$rpt_center else NA

          html <- c(
            "<!DOCTYPE html><html><head><meta charset='UTF-8'>",
            paste0("<title>Ninetails Guppy Report - ", Sys.Date(), "</title>"),
            "<style>",
            "body { font-family: 'Helvetica Neue', sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; }",
            "h1 { color: #2C7BB6; } h2 { color: #333333; border-bottom: 2px solid #dce2e8; padding-bottom: 6px; }",
            "h3 { color: #2C7BB6; margin-top: 30px; }",
            "img { max-width: 100%; margin: 10px 0; display: block; }",
            ".stat { display: inline-block; padding: 10px 20px; margin: 5px; background: #f0f5fa; border-radius: 6px; border-left: 4px solid #2C7BB6; }",
            ".stat b { font-size: 20px; color: #2C7BB6; }",
            ".note { background: #f5f7f9; padding: 8px 12px; border-radius: 4px; font-size: 13px; color: #666666; margin: 5px 0 15px 0; }",
            ".absent { background: #f5f5f5; padding: 20px; text-align: center; color: #999999; border-radius: 6px; margin: 5px 0; }",
            ".legacy { background: #fff3e0; padding: 8px 12px; border-radius: 4px; border-left: 4px solid #ff6600; margin-bottom: 15px; }",
            "</style></head><body>",
            "<h1>Ninetails Analysis Report (Guppy Legacy)</h1>",
            "<div class='legacy'>This report was generated from the Guppy legacy pipeline (fast5 + Nanopolish).</div>",
            paste0("<p>Generated: ", Sys.time(), " | v", utils::packageVersion("ninetails"), "</p>"),
            "<h2>Summary</h2>",
            paste0("<div class='stat'><b>", n_samples_val, "</b><br>Samples</div>"),
            paste0("<div class='stat'><b>", format(n_reads_val, big.mark = ","), "</b><br>Total reads</div>"),
            paste0("<div class='stat'><b>", format(n_blank_val, big.mark = ","), "</b><br>Blank</div>"),
            paste0("<div class='stat'><b>", format(n_decorated_val, big.mark = ","), "</b><br>Decorated</div>")
          )

          shiny::incProgress(0.1, detail = "Plots...")

          # Nanopolish QC
          if (isTRUE(input$rpt_qc)) {
            tryCatch({
              qc <- ninetails::nanopolish_qc(class_data, grouping_factor = gf)
              p <- ninetails::plot_nanopolish_qc(qc, frequency = TRUE)
              b64 <- .plot_to_base64(p)
              if (!is.null(b64)) html <- c(html, "<h2>Nanopolish QC</h2>",
                                           paste0("<img src='", b64, "'>"),
                                           "<div class='note'>Nanopolish QC tag distribution.</div>")
            }, error = function(e) NULL)
          }

          # Classification
          if (isTRUE(input$rpt_classification)) {
            tryCatch({
              p <- ninetails::plot_class_counts(class_data = class_data,
                                                grouping_factor = gf, frequency = FALSE, type = "N") +
                ggplot2::theme(plot.subtitle = ggplot2::element_blank())
              b64 <- .plot_to_base64(p)
              if (!is.null(b64)) html <- c(html, "<h2>Read Classification</h2>",
                                           paste0("<img src='", b64, "'>"))
            }, error = function(e) NULL)
          }

          # Abundance
          if (isTRUE(input$rpt_abundance) && has_residue) {
            tryCatch({
              p <- ninetails::plot_nonA_abundance(residue_data, grouping_factor = gf)
              b64 <- .plot_to_base64(p)
              if (!is.null(b64)) html <- c(html, "<h2>Non-A Abundance</h2>",
                                           paste0("<img src='", b64, "'>"))
            }, error = function(e) NULL)
          }

          # Residue frequency
          if (isTRUE(input$rpt_residue) && has_residue) {
            tryCatch({
              p <- ninetails::plot_residue_counts(residue_data, grouping_factor = gf,
                                                  frequency = TRUE, by_read = FALSE)
              b64 <- .plot_to_base64(p)
              if (!is.null(b64)) html <- c(html, "<h2>Residue Frequency</h2>",
                                           paste0("<img src='", b64, "'>"))
            }, error = function(e) NULL)
          }

          # Rug plots
          if (isTRUE(input$rpt_rug) && has_residue) {
            shiny::incProgress(0.2, detail = "Rug plots...")
            html <- c(html, "<h2>Non-A Position Distribution</h2>")
            sample_names_rug <- if ("sample_name" %in% names(residue_data)) {
              sort(unique(as.character(residue_data$sample_name)))
            } else { "all" }
            for (sn in sample_names_rug) {
              rd_s <- if (sn == "all") residue_data else {
                residue_data[residue_data$sample_name == sn, , drop = FALSE]
              }
              if (nrow(rd_s) == 0) next
              html <- c(html, paste0("<h3>", if (sn == "all") "All data" else sn, "</h3>"))
              for (base in c("C", "G", "U")) {
                n_total <- sum(rd_s$prediction == base, na.rm = TRUE)
                if (n_total == 0) {
                  html <- c(html, paste0("<div class='absent'>No ", base, " residues</div>"))
                  next
                }
                rd_plot <- rd_s; downsampled <- FALSE
                if (n_total > 1000) {
                  set.seed(42)
                  rd_base <- rd_s[rd_s$prediction == base, , drop = FALSE]
                  rd_other <- rd_s[rd_s$prediction != base, , drop = FALSE]
                  rd_plot <- rbind(rd_other, rd_base[sample.int(n_total, 1000, replace = FALSE), ])
                  downsampled <- TRUE
                }
                tryCatch({
                  p <- ninetails::plot_rug_density(rd_plot, base = base, max_length = 200)
                  b64 <- .plot_to_base64(p, 5.5, 3.5)
                  if (!is.null(b64)) {
                    html <- c(html, paste0("<img src='", b64, "'>"))
                    note <- if (downsampled) {
                      paste0(base, ": downsampled from ", format(n_total, big.mark = ","), " to 1,000 points.")
                    } else {
                      paste0(base, ": ", format(n_total, big.mark = ","), " points (no downsampling).")
                    }
                    html <- c(html, paste0("<div class='note'>", note, "</div>"))
                  }
                }, error = function(e) {
                  html <<- c(html, paste0("<div class='absent'>", base, " error: ", e$message, "</div>"))
                })
              }
            }
          }

          # Poly(A)
          if (isTRUE(input$rpt_polya)) {
            shiny::incProgress(0.1, detail = "Poly(A)...")
            tryCatch({
              base_data <- if (has_merged) merged_data else class_data
              p <- ninetails::plot_tail_distribution(input_data = base_data,
                                                     grouping_factor = gf, max_length = 200,
                                                     value_to_show = rpt_center_val, ndensity = TRUE) +
                ggplot2::scale_color_manual(values = palettes[[input$rpt_palette]])
              b64 <- .plot_to_base64(p)
              if (!is.null(b64)) html <- c(html, "<h2>Poly(A) Length Distribution</h2>",
                                           paste0("<img src='", b64, "'>"))
            }, error = function(e) NULL)
          }

          html <- c(html, "<hr>",
                    paste0("<p style='color: #999999; font-size: 11px; text-align: center;'>",
                           "Generated by ninetails v", utils::packageVersion("ninetails"),
                           " (Guppy legacy) | ", Sys.time(), "</p>"),
                    "</body></html>")
          writeLines(paste(html, collapse = "\n"), file)
        }) # withProgress
      },
      contentType = "text/html"
    )
  }

} # server

shiny::shinyApp(ui = ui, server = server)
