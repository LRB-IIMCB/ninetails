################################################################################
# Ninetails - cDNA orientation validation browser
#
# Minimum-useful validation tool for verifying cDNA orientation
# classifications against raw signal layout. Companion to the main
# launch_signal_browser() app; not a replacement.
#
# Launch via: ninetails::launch_cdna_signal_browser(
#   dorado_summary, pod5_dir, orientation_data
# )
#
# Inputs are populated by the launcher via shiny::shinyOptions().
# Direct app.R invocation is not supported.
################################################################################

library(shiny)

# #### Inputs from launcher ################################################

cdna_data <- shiny::getShinyOption("ninetails.cdna_data", default = NULL)
pod5_dir <- shiny::getShinyOption("ninetails.cdna_pod5_dir", default = NULL)

if (is.null(cdna_data) || is.null(pod5_dir)) {
  stop("App must be launched via ninetails::launch_cdna_signal_browser().",
       call. = FALSE)
}

# #### Colour palette (ninetails signal-segment scheme) ###################
# Matches the colours used across ninetails plotting functions
# (plot_squiggle_*, plot_tail_range_*) so the cDNA browser is visually
# consistent with the rest of the package.

COLOUR_ADAPTER <- "#089bcc" # adapter
COLOUR_TAIL <- "#f56042" # poly(A)/poly(T) tail
COLOUR_BODY <- "#3a414d" # transcript body
COLOUR_AMBIG <- "#95a5a6" # ambiguous / other

region_colours <- function(inferred_layout) {
  switch(inferred_layout,
    "polyA_layout" = list(pre = COLOUR_BODY, tail = COLOUR_TAIL, post = COLOUR_ADAPTER),
    "polyT_layout" = list(pre = COLOUR_ADAPTER, tail = COLOUR_TAIL, post = COLOUR_BODY),
    list(pre = COLOUR_AMBIG, tail = COLOUR_TAIL, post = COLOUR_AMBIG))
}

region_labels <- function(inferred_layout) {
  switch(inferred_layout,
    "polyA_layout" = list(pre = "transcript body", tail = "polyA tail", post = "adapter"),
    "polyT_layout" = list(pre = "adapter", tail = "polyT tail", post = "transcript body"),
    list(pre = "?", tail = "tail", post = "?"))
}

# #### UI ##################################################################

ui <- fluidPage(
  titlePanel("Ninetails - cDNA orientation validation browser"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("tail_type_filter",
                  "Filter by algorithm tail_type",
                  choices = c("All", "polyA", "polyT", "unidentified"),
                  selected = "All"),

      hr(),
      uiOutput("read_selector_ui"),

      fluidRow(
        column(6, actionButton("prev_btn", "- Previous", width = "100%")),
        column(6, actionButton("next_btn", "Next -", width = "100%"))
      ),

      hr(),
      htmlOutput("read_meta")
    ),

    mainPanel(
      width = 9,
      htmlOutput("agreement_header"),
      plotOutput("signal_plot", height = "420px"),
      htmlOutput("layout_legend"),
      hr(),
      htmlOutput("interpretation_hint")
    )
  )
)

# #### Server ##############################################################

server <- function(input, output, session) {

  filtered_reads <- reactive({
    if (input$tail_type_filter == "All") return(cdna_data)
    cdna_data[cdna_data$tail_type == input$tail_type_filter, , drop = FALSE]
  })

  output$read_selector_ui <- renderUI({
    fr <- filtered_reads()
    if (nrow(fr) == 0) {
      return(div(style = "color: #c33; padding: 8px;",
                 "No reads match this filter."))
    }
    choices <- setNames(
      fr$read_id,
      sprintf("%s... [%s]", substr(fr$read_id, 1, 12), fr$tail_type)
    )
    selectInput("read_select",
                "Read",
                choices = choices,
                selectize = TRUE)
  })

  observeEvent(input$prev_btn, {
    fr <- filtered_reads()
    if (is.null(input$read_select) || nrow(fr) == 0) return()
    idx <- match(input$read_select, fr$read_id)
    if (!is.na(idx) && idx > 1L) {
      updateSelectInput(session, "read_select", selected = fr$read_id[idx - 1L])
    }
  })

  observeEvent(input$next_btn, {
    fr <- filtered_reads()
    if (is.null(input$read_select) || nrow(fr) == 0) return()
    idx <- match(input$read_select, fr$read_id)
    if (!is.na(idx) && idx < nrow(fr)) {
      updateSelectInput(session, "read_select", selected = fr$read_id[idx + 1L])
    }
  })

  selected_row <- reactive({
    req(input$read_select)
    fr <- filtered_reads()
    hits <- fr[fr$read_id == input$read_select, , drop = FALSE]
    if (nrow(hits) == 0) return(NULL)
    hits[1, , drop = FALSE]
  })

  signal_data <- reactive({
    row <- selected_row()
    req(row)

    pod5_path <- ninetails:::.find_pod5_file(row$filename, pod5_dir)
    if (is.null(pod5_path)) {
      return(list(error = sprintf("POD5 file '%s' not found under %s",
                                  row$filename, pod5_dir)))
    }

    res <- tryCatch(
      ninetails:::.extract_signal_pod5(read_id = row$read_id,
                                       pod5_file = pod5_path,
                                       winsorize = TRUE),
      error = function(e) list(error = e$message)
    )
    if (is.null(res)) return(list(error = "extract_signal_pod5 returned NULL"))
    if (!is.null(res$error)) return(list(error = res$error))

    pA_signal <- (as.numeric(res$signal) + res$calibration_offset) *
      res$calibration_scale

    list(signal = pA_signal,
         signal_length = res$signal_length,
         sample_rate = res$sample_rate)
  })

  layout_and_marker <- reactive({
    row <- selected_row()
    sig <- signal_data()
    req(row)
    if (!is.null(sig$error)) return(NULL)

    inferred <- ninetails::infer_cdna_layout(
      poly_tail_start = row$poly_tail_start,
      poly_tail_end = row$poly_tail_end,
      signal_length = sig$signal_length
    )
    marker <- ninetails::cdna_layout_agreement_marker(row$tail_type, inferred)
    list(inferred = inferred, marker = marker)
  })

  output$agreement_header <- renderUI({
    row <- selected_row()
    lm <- layout_and_marker()
    req(row)
    if (is.null(lm)) return(HTML("<div style='padding:10px 0;'>(no layout - signal extraction failed)</div>"))

    marker_colour <- switch(lm$marker,
      "=" = "#228833", # check
      "x" = "#c33", # not-equal
      "?" = "#cc8800",
      "-" = "#888", # em-dash
      "#888")

    HTML(sprintf(
      paste0(
        "<div style='padding:10px 0; font-size:14px;'>",
        "<b>read</b> <code>%s</code> &nbsp;|&nbsp; ",
        "<b>reference</b> %s &nbsp;|&nbsp; ",
        "<b>algorithm</b> <code>%s</code> &nbsp;|&nbsp; ",
        "<b>layout</b> <code>%s</code> &nbsp;|&nbsp; ",
        "<b>agreement</b> ",
        "<span style='color:%s; font-size:22px; font-weight:bold;'>%s</span>",
        "</div>"),
      row$read_id,
      if ("reference" %in% colnames(row) && !is.na(row$reference)) row$reference else "(none)",
      row$tail_type,
      lm$inferred,
      marker_colour,
      lm$marker
    ))
  })

  output$read_meta <- renderUI({
    row <- selected_row()
    sig <- signal_data()
    req(row)

    if (!is.null(sig$error)) {
      return(HTML(sprintf(
        "<div style='color:#c33; font-size:12px;'>Signal extraction failed:<br>%s</div>",
        sig$error)))
    }

    pre_tail <- row$poly_tail_start - 1L
    post_tail <- sig$signal_length - row$poly_tail_end
    ratio_text <- if (min(pre_tail, post_tail) > 0) {
      sprintf("%.2f", max(pre_tail, post_tail) / min(pre_tail, post_tail))
    } else {
      "n/a"
    }

    HTML(sprintf(
      paste0("<table style='font-size:12px; width:100%%;'>",
             "<tr><td><b>signal length</b></td><td>%s</td></tr>",
             "<tr><td><b>poly_tail_start</b></td><td>%s</td></tr>",
             "<tr><td><b>poly_tail_end</b></td><td>%s</td></tr>",
             "<tr><td><b>pre-tail samples</b></td><td>%s</td></tr>",
             "<tr><td><b>post-tail samples</b></td><td>%s</td></tr>",
             "<tr><td><b>larger / smaller ratio</b></td><td>%s</td></tr>",
             "</table>"),
      format(sig$signal_length, big.mark = ","),
      format(row$poly_tail_start, big.mark = ","),
      format(row$poly_tail_end, big.mark = ","),
      format(pre_tail, big.mark = ","),
      format(post_tail, big.mark = ","),
      ratio_text
    ))
  })

  output$signal_plot <- renderPlot({
    row <- selected_row()
    sig <- signal_data()
    lm <- layout_and_marker()
    req(row)
    if (!is.null(sig$error) || is.null(lm)) {
      plot.new()
      title(main = "Signal extraction failed - see metadata panel.")
      return(invisible(NULL))
    }

    colours <- region_colours(lm$inferred)
    df <- data.frame(x = seq_along(sig$signal), y = sig$signal)
    y_range <- range(df$y, na.rm = TRUE)

    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::annotate("rect",
        xmin = 1L, xmax = row$poly_tail_start,
        ymin = y_range[1], ymax = y_range[2],
        fill = colours$pre, alpha = 0.15) +
      ggplot2::annotate("rect",
        xmin = row$poly_tail_start, xmax = row$poly_tail_end,
        ymin = y_range[1], ymax = y_range[2],
        fill = colours$tail, alpha = 0.30) +
      ggplot2::annotate("rect",
        xmin = row$poly_tail_end, xmax = sig$signal_length,
        ymin = y_range[1], ymax = y_range[2],
        fill = colours$post, alpha = 0.15) +
      ggplot2::geom_line(size = 0.2, colour = "black") +
      ggplot2::labs(x = "signal sample index",
                    y = "current (pA)",
                    title = sprintf("Inferred layout: %s", lm$inferred)) +
      ggplot2::theme_minimal(base_size = 12)
  })

  output$layout_legend <- renderUI({
    lm <- layout_and_marker()
    if (is.null(lm)) return(NULL)
    labs <- region_labels(lm$inferred)
    cols <- region_colours(lm$inferred)

    HTML(sprintf(
      paste0("<div style='padding:6px 0; font-size:12px;'>",
             "<span style='background:%s; padding:3px 8px; color:white; border-radius:3px;'>%s</span> ",
             "<span style='background:%s; padding:3px 8px; color:white; border-radius:3px;'>%s</span> ",
             "<span style='background:%s; padding:3px 8px; color:white; border-radius:3px;'>%s</span>",
             "</div>"),
      cols$pre, labs$pre,
      cols$tail, labs$tail,
      cols$post, labs$post
    ))
  })

  output$interpretation_hint <- renderUI({
    HTML(paste0(
      "<div style='font-size:11px; color:#666; padding:8px 0;'>",
      "Agreement markers: ",
      "<b>=</b> algorithm and layout agree &nbsp;|&nbsp; ",
      "<b>x</b> they disagree &nbsp;|&nbsp; ",
      "<b>?</b> layout ambiguous (pre/post regions too small or too similar) &nbsp;|&nbsp; ",
      "<b>-</b> algorithm returned unidentified, nothing to compare against.",
      "<br>Layout inference defaults: minimum region 100 samples, strictness ratio 1.5.",
      "</div>"
    ))
  })
}

shinyApp(ui, server)
