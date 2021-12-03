#' Given a numeric matrix, assign to each cell a color (character) value based on linearly
#' interpolating a given vector of colors.
#'
#' This is an improved version of paintmap::colour_matrix() function. Original one threw
#' "‘breaks’ are not unique" error in cut() if the quantile values are the same in multiple
#' cols.
#'
#' @param gasf numeric matrix.
#' @param colors character vector of colors.
#'
#' @return character matrix - gasf with hex colors.
#' @export
#'
#' @examples
#'\dontrun{
#'
#' color_matrix_hex(gasf=gasf, colors=rainbow(100))
#'
#'}
#'
#'
color_matrix_hex <- function(gasf, colors=rainbow(100)) {
  centres <- seq(from=min(gasf), to=max(gasf), length.out=length(colors))
  interval_width <- 1/length(colors)
  matrix(
    as.character(factor(levels=colors, cut(as.numeric(gasf), labels=colors, breaks=unique(c(centres[1] - interval_width/2, centres + interval_width/2))))),
    nrow(gasf),
    ncol(gasf),
    dimnames=dimnames(gasf)
  )
}

