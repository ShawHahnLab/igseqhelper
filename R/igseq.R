#' Parse gene names from allele names
#' @param txt character vector of immunoglobulin allele names
#' @return character vector of gene names
#' @export
parse_gene <- function(txt) {
  gsub("^I?G?([^\\*]+)\\*.*", "\\1", txt)
}

#' Fetch metadata from Google Sheets
#'
#' This takes a nested list structure with identifiers for Google Sheets and
#' saves each sheet to a CSV file.  Each entry in the list should have a short
#' workbook name (used as a directory name for the contained sheets) and should
#' itself be a list containing two things: one entry named "slug" set to the
#' long alphanumeric string for the workbook, and a second entry named "sheets"
#' containing the integer GID for each sheet to download with names set to the
#' names of the sheets.
#'
#' @param google_sheets nested list structure like in
#'   \code{\link{GOOGLE_SHEETS}}
#' @param prefix directory to put CSV files in.  If \code{NULL}, defaults to
#' "metadata/" in the current working directory. 
#' @export
update_metadata <- function(google_sheets, prefix=NULL) {
  if (is.null(prefix)) {
    prefix <- "metadata"
  }
  for (workbook in names(google_sheets)) {
    slug <- google_sheets[[workbook]][["slug"]]
    path_dir <- file.path(prefix, workbook)
    if (! dir.exists(path_dir)) {
      dir.create(path_dir, recursive = TRUE)
    }
    for (sheet in names(google_sheets[[workbook]][["sheets"]])) {
      gid <- google_sheets[[workbook]][["sheets"]][[sheet]]
      url <- paste0(
        "https://docs.google.com/spreadsheets/d/e/", slug,
        "/pub?gid=", gid,
        "&single=true&output=csv")
      path_file <- file.path(path_dir, paste0(sheet, ".csv"))
      data <- read.csv(file = url)
      write.csv(data, path_file, quote = FALSE, row.names = FALSE)
    }
  }
}
