#' Plot Identity/Divergence
#'
#' As in SONAR's id/div plots, comparing each sequence to divergence from a
#' germline sequence (x-axis) and identity to a mature reference (y-axis).
#'
#' @param al named character vector of already-aligned sequences
#' @param mature name of mature reference within al
#' @param germline name of germline reference within al
#' @param bins number of bins for density plot in either dimension
#'
#' @return ggplot2 plot object
ggplot_id_div <- function(al, mature, germline, bins=60) {
  dists_mature <- levenR::leven(al, al[mature])
  dists_germline <- levenR::leven(al, al[germline])
  pident_germline <- (nchar(al[1]) - dists_germline) / nchar(al[1])
  pident_mature <- (nchar(al[1]) - dists_mature) / nchar(al[1])

  df <- data.frame(
    Divergence = 100 * (1 - pident_germline),
    Identity = 100 * pident_mature
  )

  max_log <- ceiling(log10(max(table(paste(
    cut(df$Divergence, breaks = bins),
    cut(df$Identity, breaks = bins))))))
  breaks <- 10 ^ (0:max_log)

  ggplot2::ggplot(df, aes(x = Divergence, y = Identity)) +
    ggplot2::geom_bin2d(ggplot2::aes(fill = stat(count)), bins = bins) +
    ggplot2::scale_fill_gradientn(
      colors = c("blue", "green", "red"),
      breaks = breaks,
      labels = breaks,
      limits = range(breaks),
      trans = "log") +
    ggplot2::labs(
      x = paste("% divergence from", names(al[germline])),
      y = paste("% ID to", names(al)[mature]),
      fill = "number of\nsequences") +
    ggplot2::theme_bw()
}
