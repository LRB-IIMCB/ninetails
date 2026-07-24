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

# #### Primer schematic helpers (sequence-space tail/primer view) #########
# Locates SSP/VNP primers in the read sequence and draws a proportional
# schematic. Orientation-aware: polyA reads expect SSP at the 5' end and the
# reverse complement of VNP at the 3' end; polyT reads use the mirror layout.
# Primers are located with utils::adist (Levenshtein, tolerant of the indels
# real primers carry); reverse complements reuse ninetails::reverse_complement()
# so the layout matches detect_orientation_single().

# Dorado primers, copied verbatim from detect_orientation_single()
PRIMER_FRONT <- "TTTCTGTTGGTGCTGATATTGCTTT"  # SSP
PRIMER_REAR <- "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAGTTTT"  # VNP
PRIMER_WINDOW <- 150L

# Best single match of `primer` within [search_start, search_end] of `sequence`.
# Returns position (full-read coords), edit distance, and the matched substring.
find_primer_best_match <- function(sequence, primer,
                                   search_start = 1L,
                                   search_end = nchar(sequence)) {
  seq_len <- nchar(sequence)
  primer_len <- nchar(primer)
  search_start <- max(1L, as.integer(search_start))
  search_end <- min(seq_len, as.integer(search_end))

  # Guard: window shorter than primer -> no valid placement
  if (is.na(sequence) || search_end - search_start + 1L < primer_len) {
    return(list(position = NA_integer_, distance = NA_integer_, matched_seq = NA_character_))
  }

  starts <- seq.int(search_start, search_end - primer_len + 1L)
  windows <- substring(sequence, starts, starts + primer_len - 1L)
  dists <- as.integer(utils::adist(primer, windows)[1, ])  # single vectorized call
  best <- which.min(dists)

  list(position = starts[best], distance = dists[best], matched_seq = windows[best])
}

# Orientation-aware extraction. `tail_type` is the algorithm call
# ("polyA"/"polyT"/anything else -> best hypothesis, flagged unclassified).
extract_oriented_primers <- function(sequence, tail_type,
                                     front_primer = PRIMER_FRONT,
                                     rear_primer = PRIMER_REAR,
                                     search_window = PRIMER_WINDOW) {
  seq_len <- nchar(sequence)
  rear_trimmed <- sub("T+$", "", rear_primer)  # trim trailing Ts, as in source
  front_rc <- ninetails::reverse_complement(front_primer)
  rear_rc <- ninetails::reverse_complement(rear_trimmed)

  top_start <- 1L
  top_end <- min(search_window, seq_len)
  bot_start <- max(1L, seq_len - search_window + 1L)
  bot_end <- seq_len

  score_hypothesis <- function(five_primer, three_primer) {
    five <- find_primer_best_match(sequence, five_primer, top_start, top_end)
    three <- find_primer_best_match(sequence, three_primer, bot_start, bot_end)
    list(five = five, three = three,
         total_dist = sum(five$distance, three$distance, na.rm = TRUE))
  }
  hypo_a <- score_hypothesis(front_primer, rear_rc)   # SSP@5', RC(VNP)@3'
  hypo_t <- score_hypothesis(rear_trimmed, front_rc)  # VNP@5', RC(SSP)@3'

  if (identical(tail_type, "polyA")) {
    chosen <- hypo_a; five_name <- "SSP"; three_name <- "RC(VNP)"; kind <- "polyA"
  } else if (identical(tail_type, "polyT")) {
    chosen <- hypo_t; five_name <- "VNP"; three_name <- "RC(SSP)"; kind <- "polyT"
  } else if (hypo_a$total_dist <= hypo_t$total_dist) {
    chosen <- hypo_a; five_name <- "SSP"; three_name <- "RC(VNP)"; kind <- "unclassified"
  } else {
    chosen <- hypo_t; five_name <- "VNP"; three_name <- "RC(SSP)"; kind <- "unclassified"
  }

  list(
    classification = kind,
    seq_length = seq_len,
    five = list(name = five_name, position = chosen$five$position,
                distance = chosen$five$distance, matched = chosen$five$matched_seq),
    three = list(name = three_name, position = chosen$three$position,
                 distance = chosen$three$distance, matched = chosen$three$matched_seq)
  )
}

# Two-panel proportional schematic. Top: read drawn strictly to scale. Bottom:
# both ends to scale with the long body compressed behind an axis break.
plot_primer_schematic <- function(sequence, primers,
                                  search_window = PRIMER_WINDOW,
                                  gap_render = 40) {
  seq_len <- nchar(sequence)
  p5 <- primers$five$position
  p5e <- if (!is.na(p5)) p5 + nchar(primers$five$matched) - 1L else NA_integer_
  p3 <- primers$three$position
  p3e <- if (!is.na(p3)) p3 + nchar(primers$three$matched) - 1L else NA_integer_

  # Real segments [a, b], inclusive; drop empty/NA ones
  segs <- list()
  push <- function(role, a, b, name = NA_character_, ed = NA_integer_) {
    if (is.na(a) || is.na(b) || b < a) return(invisible())
    segs[[length(segs) + 1L]] <<- data.frame(
      role = role, a = a, b = b, name = name, ed = ed,
      stringsAsFactors = FALSE)
  }
  push("flank", 1L, p5 - 1L)
  push("primer", p5, p5e, primers$five$name, primers$five$distance)
  push("body", p5e + 1L, p3 - 1L)
  push("primer", p3, p3e, primers$three$name, primers$three$distance)
  push("flank", p3e + 1L, seq_len)
  seg_df <- do.call(rbind, segs)
  if (is.null(seg_df)) return(NULL)

  # Bottom-panel piecewise remap: left zone at scale, body compressed, right at scale
  body_lo <- min(search_window, seq_len)
  body_hi <- max(1L, seq_len - search_window + 1L) - 1L
  has_break <- body_hi > body_lo
  remap <- function(x) {
    if (!has_break) return(x)
    if (x <= body_lo) return(x)
    if (x >= body_hi) return(body_lo + gap_render + (x - body_hi))
    body_lo + gap_render * (x - body_lo) / (body_hi - body_lo)
  }

  panel_levels <- c("True scale (proportional)", "Zoomed ends (body compressed)")

  prop <- data.frame(
    panel = panel_levels[1], role = seg_df$role,
    xmin = seg_df$a - 1L, xmax = seg_df$b,
    stringsAsFactors = FALSE)
  brk <- data.frame(
    panel = panel_levels[2], role = seg_df$role,
    xmin = vapply(seg_df$a - 1L, remap, numeric(1)),
    xmax = vapply(seg_df$b, remap, numeric(1)),
    stringsAsFactors = FALSE)
  rect_df <- rbind(prop, brk)
  rect_df$panel <- factor(rect_df$panel, levels = panel_levels)
  rect_df <- rect_df[order(rect_df$role == "primer"), ]  # primers drawn on top

  prim <- seg_df[seg_df$role == "primer", , drop = FALSE]
  prim$mid <- (prim$a + prim$b) / 2
  prim$label <- paste0(prim$name, "\ned=", prim$ed)
  lab_df <- rbind(
    data.frame(panel = panel_levels[1], x = prim$mid, label = prim$label,
               stringsAsFactors = FALSE),
    data.frame(panel = panel_levels[2],
               x = vapply(prim$mid, remap, numeric(1)), label = prim$label,
               stringsAsFactors = FALSE))
  lab_df$panel <- factor(lab_df$panel, levels = panel_levels)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = rect_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = role),
      colour = "white") +
    ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(x = x, y = 1.18, label = label),
      size = 3, lineheight = 0.9, vjust = 0) +
    ggplot2::scale_fill_manual(
      values = c(primer = COLOUR_ADAPTER, body = COLOUR_BODY, flank = COLOUR_AMBIG),
      breaks = c("primer", "body", "flank"),
      labels = c("primer", "transcript body", "flank outside primers"),
      name = NULL) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x") +
    ggplot2::scale_y_continuous(limits = c(-0.4, 1.6), expand = c(0, 0)) +
    ggplot2::labs(x = "sequence position (bp) / rendered position", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom")

  # Axis-break glyph + body-length annotation (bottom panel only)
  if (has_break) {
    body_bp <- (p3 - 1L) - (p5e + 1L) + 1L
    bx <- body_lo + gap_render / 2
    seg_glyph <- data.frame(
      panel = factor(panel_levels[2], levels = panel_levels),
      x = c(bx - 3, bx - 1), xend = c(bx + 1, bx + 3),
      y = c(-0.1, -0.1), yend = c(1.1, 1.1),
      stringsAsFactors = FALSE)
    txt_glyph <- data.frame(
      panel = factor(panel_levels[2], levels = panel_levels),
      x = bx, y = -0.3,
      label = sprintf("~ %s bp body", format(body_bp, big.mark = ",")),
      stringsAsFactors = FALSE)
    p <- p +
      ggplot2::geom_segment(
        data = seg_glyph,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        colour = "#000000") +
      ggplot2::geom_text(
        data = txt_glyph,
        ggplot2::aes(x = x, y = y, label = label), size = 2.8)
  }

  p
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
      htmlOutput("interpretation_hint"),
      hr(),
      h4("Primer schematic (sequence space)"),
      plotOutput("primer_schematic", height = "300px"),
      htmlOutput("primer_table")
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
      return(div(style = "color: #cc3333; padding: 8px;",
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
                            "x" = "#cc3333", # not-equal
                            "?" = "#cc8800",
                            "-" = "#888888", # em-dash
                            "#888888")

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
        "<div style='color:#cc3333; font-size:12px;'>Signal extraction failed:<br>%s</div>",
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
      "<div style='font-size:11px; color:#666666; padding:8px 0;'>",
      "Agreement markers: ",
      "<b>=</b> algorithm and layout agree &nbsp;|&nbsp; ",
      "<b>x</b> they disagree &nbsp;|&nbsp; ",
      "<b>?</b> layout ambiguous (pre/post regions too small or too similar) &nbsp;|&nbsp; ",
      "<b>-</b> algorithm returned unidentified, nothing to compare against.",
      "<br>Layout inference defaults: minimum region 100 samples, strictness ratio 1.5.",
      "</div>"
    ))
  })

  # #### Primer schematic (sequence-space view) ###########################

  has_sequence <- reactive("sequence" %in% colnames(cdna_data))

  primer_extract <- reactive({
    row <- selected_row()
    req(row)
    if (!has_sequence()) return(NULL)
    seqv <- row$sequence
    if (is.na(seqv) || nchar(seqv) < 50L) return(list(too_short = TRUE))
    extract_oriented_primers(seqv, row$tail_type)
  })

  output$primer_schematic <- renderPlot({
    validate(need(has_sequence(),
                  "Primer schematic unavailable: orientation_data has no 'sequence' column."))
    pe <- primer_extract()
    validate(need(!is.null(pe) && !isTRUE(pe$too_short),
                  "Sequence missing or too short (< 50 bp) for primer search."))
    row <- selected_row()
    plot_primer_schematic(row$sequence, pe)
  })

  output$primer_table <- renderUI({
    if (!has_sequence()) return(NULL)
    pe <- primer_extract()
    if (is.null(pe) || isTRUE(pe$too_short)) return(NULL)
    HTML(sprintf(
      paste0("<table style='font-size:12px; width:100%%;'>",
             "<tr><td><b>primer-based call</b></td><td><code>%s</code></td></tr>",
             "<tr><td><b>%s (5')</b></td><td>pos %s, edit dist %s</td></tr>",
             "<tr><td style='padding-left:14px;'>matched</td><td><code>%s</code></td></tr>",
             "<tr><td><b>%s (3')</b></td><td>pos %s, edit dist %s</td></tr>",
             "<tr><td style='padding-left:14px;'>matched</td><td><code>%s</code></td></tr>",
             "</table>",
             "<div style='font-size:11px; color:#666666; padding:4px 0;'>",
             "Schematic is in sequence (nucleotide) space; the poly(A/T) tail ",
             "is not drawn here (its coordinates are signal-space - see the ",
             "signal plot above).</div>"),
      pe$classification,
      pe$five$name, pe$five$position, pe$five$distance, pe$five$matched,
      pe$three$name, pe$three$position, pe$three$distance, pe$three$matched
    ))
  })
}

shinyApp(ui, server)
