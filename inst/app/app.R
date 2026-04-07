################################################################################
# Ninetails Analysis Dashboard
#
# Location: inst/app/app.R
# Launch:   ninetails::launch_signal_browser()
#
# Layout: fluidPage with header bar + tabsetPanel
#   Tab 1: Summary        — value boxes (samples, transcripts, reads)
#   Tab 2: Classification — plot_class_counts, plot_nonA_abundance
#   Tab 3: Residues       — plot_residue_counts, summary DT table
#   Tab 4: Poly(A)        — plot_tail_distribution with palette picker
#   Tab 5: Signal Viewer  — squiggle plots with sub-tabs (Viewer / Explorer)
#
# Static assets (inst/app/www/):
#   - logo.png    (copy from man/figures/logo.png)
#   - favicon.ico (copy from pkgdown/favicon/)
#
################################################################################

library(shiny)
library(ggplot2)
library(dplyr)
library(vroom)
library(reticulate)


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

# Startup diagnostics (visible in R console)
cat("--- Ninetails Dashboard Startup ---\n")
cat("  class_data:    ", if (has_class) paste(nrow(class_data), "rows,",
                                              ncol(class_data), "cols") else "NULL", "\n")
cat("  residue_data:  ", if (has_residue) paste(nrow(residue_data), "rows")
    else "NULL", "\n")
cat("  merged_data:   ", if (has_merged) paste(nrow(merged_data), "rows")
    else "NULL", "\n")
cat("  signal_config: ", length(signal_config), "sample(s)\n")
if (has_class) cat("  columns:       ",
                   paste(names(class_data), collapse = ", "), "\n")

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
        background-color: #dce2e8; color: #333; padding: 16px 20px;
        margin: -15px -15px 20px -15px;
        display: flex; align-items: center; gap: 14px;
      }
      .header h2 { margin: 0; font-weight: bold; font-size: 22px; }
      .header p { margin: 4px 0 0 0; color: #555; font-size: 14px; }
      .card {
        background: #fff; border-radius: 6px; padding: 15px;
        margin-bottom: 15px; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
      }
      .card h4 {
        margin-top: 0; color: #333; font-weight: 600;
        border-bottom: 2px solid #2C7BB6; padding-bottom: 8px;
      }
      .value-box {
        background: #fff; border-radius: 6px; padding: 20px;
        text-align: center; box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        border-left: 5px solid #2C7BB6; margin-bottom: 15px;
      }
      .value-box .vb-value { font-size: 28px; font-weight: 700; color: #2C7BB6; }
      .value-box .vb-label { font-size: 13px; color: #666; margin-top: 4px; }
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
    # TAB 1: Summary
    ############################################################################
    shiny::tabPanel("Summary", shiny::br(),
                    if (has_class) {
                      shiny::fluidRow(
                        shiny::column(4,
                                      shiny::div(class = "value-box",
                                                 shiny::div(class = "vb-value",
                                                            shiny::textOutput("n_samples", inline = TRUE)),
                                                 shiny::div(class = "vb-label", "Samples analyzed"))),
                        shiny::column(4,
                                      shiny::div(class = "value-box",
                                                 shiny::div(class = "vb-value",
                                                            shiny::textOutput("n_transcripts", inline = TRUE)),
                                                 shiny::div(class = "vb-label", "Transcripts found"))),
                        shiny::column(4,
                                      shiny::div(class = "value-box",
                                                 shiny::div(class = "vb-value",
                                                            shiny::textOutput("n_reads", inline = TRUE)),
                                                 shiny::div(class = "vb-label", "Total reads")))
                      ) # fluidRow
                    } else {
                      shiny::div(class = "placeholder-msg",
                                 shiny::h4("No analysis data loaded"),
                                 shiny::p("Provide class_file or a YAML config to enable this tab."))
                    }
    ), # tabPanel Summary

    ############################################################################
    # TAB 2: Classification
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
                                                                          selected = "N")
                                            ) # div filter-section
                        ), # sidebarPanel
                        shiny::mainPanel(width = 10,
                                         shiny::div(class = "card",
                                                    shiny::h4("Read Classification"),
                                                    shiny::plotOutput("class_plot", height = "450px")),
                                         if (has_residue) {
                                           shiny::div(class = "card",
                                                      shiny::h4("Non-A Abundance (reads with 1, 2, 3+ non-As)"),
                                                      shiny::plotOutput("nonA_abundance_plot", height = "400px"))
                                         }
                        ) # mainPanel
                      ) # sidebarLayout
                    } else {
                      shiny::div(class = "placeholder-msg",
                                 shiny::h4("No classification data loaded"),
                                 shiny::p("Provide class_file or a YAML config."))
                    }
    ), # tabPanel Classification

    ############################################################################
    # TAB 3: Residues
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
                                                       shiny::checkboxInput("residue_by_read", "Count by read", FALSE)
                                            ) # div filter-section
                        ), # sidebarPanel
                        shiny::mainPanel(width = 10,
                                         shiny::div(class = "card",
                                                    shiny::h4("Residue Counts"),
                                                    shiny::plotOutput("residue_plot", height = "450px")),
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
    # TAB 4: Poly(A) Distribution
    ############################################################################
    shiny::tabPanel("Poly(A)", shiny::br(),
                    if (has_class) {
                      shiny::sidebarLayout(
                        shiny::sidebarPanel(width = 2,
                                            shiny::div(class = "filter-section",
                                                       shiny::selectInput("polya_grouping", "Grouping variable",
                                                                          choices = available_groups, selected = available_groups[1]),
                                                       shiny::selectizeInput("polya_contig_filter", "Filter by transcript",
                                                                             choices = NULL,
                                                                             options = list(placeholder = "All transcripts", maxOptions = 200)),
                                                       shiny::selectInput("polya_center", "Central tendency",
                                                                          choices = c("mean", "median", "mode", "none"), selected = "none"),
                                                       shiny::sliderInput("polya_max_length", "Max tail length",
                                                                          min = 0, max = 500, value = 200),
                                                       shiny::checkboxInput("polya_ndensity", "Normalized density", TRUE),
                                                       shiny::selectInput("polya_palette", "Color palette",
                                                                          choices = names(palettes), selected = "unova")
                                            ) # div filter-section
                        ), # sidebarPanel
                        shiny::mainPanel(width = 10,
                                         shiny::div(class = "card",
                                                    shiny::h4("Poly(A) Tail Length Distribution"),
                                                    shiny::plotOutput("polya_dist_plot", height = "500px"))
                        ) # mainPanel
                      ) # sidebarLayout
                    } else {
                      shiny::div(class = "placeholder-msg",
                                 shiny::h4("No data loaded"),
                                 shiny::p("Provide class_file or a YAML config."))
                    }
    ), # tabPanel Poly(A)

    ############################################################################
    # TAB 5: Signal Viewer
    ############################################################################
    shiny::tabPanel("Signal Viewer", shiny::br(),
                    shiny::fluidRow(

                      # ---- Left panel ----
                      shiny::column(3,

                                    shiny::div(class = "card",
                                               shiny::h4("Data"),
                                               if (length(signal_config) > 0 && names(signal_config)[1] != "single") {
                                                 # Multi mode: sample dropdown only (paths pre-loaded)
                                                 shiny::div(class = "paths-section",
                                                            shiny::selectInput("signal_sample", "Select sample",
                                                                               choices = names(signal_config),
                                                                               selected = names(signal_config)[1]))
                                               } else if (length(signal_config) == 0) {
                                                 # No pre-loaded data: manual path entry
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
                                               # Single pre-loaded: no paths shown, auto-loads in server
                                               shiny::uiOutput("signal_data_status")
                                    ), # div card data

                                    shiny::div(class = "card",
                                               shiny::h4("Filters"),
                                               shiny::div(class = "filter-section",
                                                          shiny::numericInput("min_polya_length",
                                                                              "Minimum poly(A) length", value = 10, min = 0, step = 1),
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

                                                       shiny::tabPanel("Signal Viewer", shiny::br(),
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
                                                       ), # tabPanel Signal Viewer

                                                       shiny::tabPanel("Signal Explorer", shiny::br(),
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
                                                       ) # tabPanel Signal Explorer
                                    ) # tabsetPanel signal_sub_tabs
                      ) # column 9
                    ) # fluidRow
    ) # tabPanel Signal Viewer

  ) # tabsetPanel main_tabs
) # fluidPage


################################################################################
# SERVER
################################################################################

server <- function(input, output, session) {

  ##############################################################################
  # SUMMARY TAB
  ##############################################################################

  if (has_class) {
    output$n_samples <- shiny::renderText({
      if ("sample_name" %in% names(class_data)) {
        length(unique(class_data$sample_name))
      } else { 1 }
    })
    output$n_transcripts <- shiny::renderText({
      length(unique(class_data[[transcript_col]]))
    })
    output$n_reads <- shiny::renderText({
      format(nrow(class_data), big.mark = ",")
    })
  }

  ##############################################################################
  # SHARED: contig filters (server-side selectize)
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

    output$class_plot <- shiny::renderPlot({
      cd <- class_data_filtered(); shiny::req(nrow(cd) > 0)
      shiny::req(input$class_grouping, input$class_plot_type, !is.null(input$class_frequency))
      gf <- if (input$class_grouping %in% names(cd)) input$class_grouping else NA
      ninetails::plot_class_counts(class_data = cd,
                                   grouping_factor = gf, frequency = input$class_frequency,
                                   type = input$class_plot_type)
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
  # RESIDUES TAB
  ##############################################################################

  if (has_residue) {

    residue_data_filtered <- shiny::reactive({
      .filter_by_contig(residue_data, input$residue_contig_filter)
    })

    output$residue_plot <- shiny::renderPlot({
      rd <- residue_data_filtered(); shiny::req(nrow(rd) > 0)
      shiny::req(input$residue_grouping, !is.null(input$residue_frequency), !is.null(input$residue_by_read))
      gf <- if (input$residue_grouping %in% names(rd)) input$residue_grouping else NA
      ninetails::plot_residue_counts(residue_data = rd,
                                     grouping_factor = gf, frequency = input$residue_frequency,
                                     by_read = input$residue_by_read)
    })

    if (has_merged) {
      output$residue_summary_table <- DT::renderDT({
        shiny::req(input$residue_grouping)
        md <- .filter_by_contig(merged_data, input$residue_contig_filter)
        shiny::req(nrow(md) > 0)
        tid_col <- if ("ensembl_transcript_id_short" %in% names(md)) {
          "ensembl_transcript_id_short"
        } else { NULL }
        gf <- if (input$residue_grouping %in% names(md)) input$residue_grouping else "contig"
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

    polya_data_filtered <- shiny::reactive({
      base <- if (has_merged) merged_data else class_data
      .filter_by_contig(base, input$polya_contig_filter)
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
  }

  ##############################################################################
  # SIGNAL VIEWER TAB
  ##############################################################################

  loaded_signal_data    <- shiny::reactiveVal(NULL)
  loaded_signal_residue <- shiny::reactiveVal(NULL)
  signal_error          <- shiny::reactiveVal(NULL)

  # ---- Data loading ----

  # Multi mode with named samples
  if (length(signal_config) > 0 && names(signal_config)[1] != "single") {

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
      shiny::selectInput("genome_filter", "Alignment genome",
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
    if (!is.null(input$min_polya_length) && "poly_tail_length" %in% names(data))
      data <- dplyr::filter(data, poly_tail_length >= input$min_polya_length)
    if (!is.null(input$comments_filter) && input$comments_filter != "All" && "comments" %in% names(data))
      data <- dplyr::filter(data, comments == input$comments_filter)
    if (!is.null(input$genome_filter) && input$genome_filter != "All" && "alignment_genome" %in% names(data))
      data <- dplyr::filter(data, alignment_genome == input$genome_filter)
    if (!is.null(input$mapq_filter) && "alignment_mapq" %in% names(data))
      data <- dplyr::filter(data, alignment_mapq >= input$mapq_filter[1], alignment_mapq <= input$mapq_filter[2])
    rt <- loaded_signal_residue()
    if (!is.null(input$sig_residue_filter) && input$sig_residue_filter != "All" && !is.null(rt)) {
      # Normalize column name
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
    # Normalize column name (ninetails output uses readname, dorado uses read_id)
    if (!"read_id" %in% names(rd) && "readname" %in% names(rd)) {
      rd <- dplyr::rename(rd, read_id = readname)
    }
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
  output$full_signal_plot <- shiny::renderPlot({
    df <- signal_data(); shiny::req(df)
    ps <- attr(df, "polya_start"); pe <- attr(df, "polya_end")
    dp <- if (nrow(df) > 20000) df[round(seq(1, nrow(df), length.out = 20000)), ] else df
    ov <- ninetails:::.build_nonA_overlay(read_nonA_data = selected_residue(),
                                          poly_tail_start = ps, poly_tail_end = pe, nonA_flank = 250)
    p <- ggplot2::ggplot(dp, ggplot2::aes(x = position, y = signal, color = segment))
    for (l in ov) p <- p + l
    p + ggplot2::geom_line(linewidth = 0.3) +
      ggplot2::scale_color_manual(values = c("Adapter" = "#089bcc", "Poly(A)" = "#f56042",
                                             "Transcript" = "#3a414d", "Other" = "#95a5a6"), name = "Region") +
      ggplot2::geom_vline(xintercept = ps, color = "#700f25", linetype = "dashed", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = pe, color = "#0f3473", linetype = "dashed", linewidth = 0.8) +
      ggplot2::labs(x = "Position (samples)", y = "Signal (raw)", title = paste("Read:", attr(df, "read_id"))) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(text = ggplot2::element_text(family = "Open Sans"),
                     legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold", size = 14))
  })

  output$polya_zoom_plot <- shiny::renderPlot({
    df <- signal_data(); shiny::req(df)
    ps <- attr(df, "polya_start"); pe <- attr(df, "polya_end")
    zs <- max(1, ps - 250); ze <- min(nrow(df), pe + 250)
    dz <- dplyr::filter(df, position >= zs, position <= ze)
    ov <- ninetails:::.build_nonA_overlay(read_nonA_data = selected_residue(),
                                          poly_tail_start = ps, poly_tail_end = pe, nonA_flank = 250)
    p <- ggplot2::ggplot(dz, ggplot2::aes(x = position, y = signal, color = segment))
    for (l in ov) p <- p + l
    p + ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::scale_color_manual(values = c("Adapter" = "#089bcc", "Poly(A)" = "#f56042",
                                             "Transcript" = "#3a414d", "Other" = "#95a5a6"), name = "Region") +
      ggplot2::geom_vline(xintercept = ps, color = "#700f25", linetype = "dashed", linewidth = 1) +
      ggplot2::geom_vline(xintercept = pe, color = "#0f3473", linetype = "dashed", linewidth = 1) +
      ggplot2::labs(x = "Position (samples)", y = "Signal (raw)",
                    title = paste("Poly(A) region:", ps, "-", pe), subtitle = paste("Positions", zs, "to", ze)) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(text = ggplot2::element_text(family = "Open Sans"),
                     legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold", size = 14),
                     plot.subtitle = ggplot2::element_text(color = "#666", size = 11))
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
        fc <- rc[as.character(rr$prediction[i])]; if (is.na(fc)) fc <- "#999"
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
