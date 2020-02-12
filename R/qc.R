load_qualgrid <- function(fp_csv, trimblanks=TRUE) {
  grid <- read.csv(fp_csv, row.names = 1, header = TRUE, check.names = FALSE)
  grid <- grid[seq(nrow(grid), 1, -1), ]
  if (trimblanks) {
    # trim off an all-zero portion on the right, if there is one.
    trim <- cumsum(rev(colSums(grid))) == 0
    if (trim[1]) {
      idx <- match(FALSE, trim)
      grid <- grid[, 1:(ncol(grid) - idx + 1)]
    } 
  }
  grid
}

plot_qualgrid <- function(grid, ...) {
  grid[grid == 0] <- NA
  pheatmap::pheatmap(
    log10(grid),
    cluster_rows = FALSE, cluster_cols = FALSE,
    fontsize = 6,
    show_colnames = FALSE,
    ...)
}