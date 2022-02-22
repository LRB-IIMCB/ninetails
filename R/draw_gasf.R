#' Plotting gramian angular summation field.
#'
#' @param chunkname character string. Name of the given signal chunk within the
#' analyzed dataset.
#'
#' @param tail_chunk_list list object produced by create_tail_chunk_list function.
#'
#' @return an array (100,100,1) with values (greyscale) representing ONT signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_gasf(chunkname = '1234-jhgfdr54-io643gju-ouy4378989', tail_chunk_list = tail_chunk_list)
#'
#'}
#'
#'
#'
draw_gasf <- function(chunkname, gasf_list){

  gasf <- gasf_list[[chunkname]]

  gasf <- reshape2::melt(gasf)

  plt <- gasf %>% ggplot2::ggplot(aes(Var2, Var1, fill=value)) + ggplot2::geom_raster() +ggplot2::scale_fill_gradientn(colours = rainbow(100), guide=FALSE) +
    ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   panel.background = element_blank(),
                   plot.margin=grid::unit(c(0,0,0,0), "mm"))

  ggsave(paste0(chunkname, ".png"), device = "png", bg = "transparent")

  #return(plt)
}
