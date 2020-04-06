# General data and metadata handling


# Named Metadata ----------------------------------------------------------


#' Load sequence metadata CSV
#' @export
load_sequences <- function(fp) {
  sequences <- load_csv(fp, "Name")
  sequences$Seq <- gsub(" ", "", sequences$Seq)
  if (nrow(sequences) != length(unique(sequences$Seq))) {
    stop("Non-unique sequences in input CSV")
  }
  sequences
}

#' Load run metadata CSV
#' @export
load_runs <- function(fp) {
  load_csv(fp, "Run")
}

#' Load specimen metadata CSV
#' @export
load_specimens <- function(fp) {
  data <- load_csv(fp, "Specimen")
  data
}

#' Load sample metadata CSV
#' @export
load_samples <- function(fp, specimens=NULL, runs=NULL, sequences=NULL) {
  samples <- load_csv(fp, "Sample")
  if (! is.null(specimens)) {
    samples <- merge_df(samples, specimens, "Specimen")
  }
  if (! is.null(runs)) {
    samples <- merge_df(samples, runs, "Run")
  }
  if (! is.null(sequences)) {
    samples <- merge_df(samples, sequences, "BarcodeFwd", "Name")
    samples <- merge_df(samples, sequences, "BarcodeRev", "Name")
  }
  samples
}

load_antibody_lineages <- function(fp) {
  antibody_lineages <- load_csv(fp, "AntibodyLineage")
  antibody_lineages
}

load_antibody_isolates <- function(fp, antibody_lineages=NULL) {
  antibody_isolates <- load_csv(fp, "AntibodyIsolate")
  if (! is.null(antibody_lineages)) {
    antibody_isolates <- merge_df(
      antibody_isolates, antibody_lineages, "AntibodyLineage")
  }
  antibody_isolates
}

merge_df <- function(df1, df2, key, name=key) {
  idx <- match(df1[[key]], df2[[name]])
  df2 <- df2[idx, ]
  df2 <- df2[, -match(name, colnames(df2))]
  colnames(df2) <- paste0(key, colnames(df2))
  if (any(colnames(df2) %in% colnames(df1))) {
    stop("Column name collision between data frames")
  }
  df1 <- cbind(df1, df2)
  df1
}


# Google Sheets -----------------------------------------------------------


#' Fetch metadata from Google Sheets
#'
#' This takes a nested list structure with identifiers for Google Sheets and
#' saves each sheet to a CSV file.  Each entry in the list should have a short
#' workbook name (used as a directory name for the contained sheets) and should
#' itself be a list containing two things: one entry named "slug" set to the
#' long alphanumeric string for the workbook, and a second entry named "sheets"
#' containing the integer GID for each sheet to download with names set to the
#' names of the sheets.  (The slug and the sheet GID are to the left and right
#' of "/pub?gid=" in the published URL.
#'
#' @param google_sheets nested list structure
#' @param prefix directory to put CSV files in.  If \code{NULL}, defaults to
#' "metadata/" in the current working directory.
#' @export
#' @describeIn update_metadata Update metadata based on list
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
      data <- read.csv(
        file = url,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE)
      write.csv(data, path_file, row.names = FALSE)
    }
  }
}

#' @describeIn update_metadata Update metadata based on list
#' @param fp file path to YAML providing Google Sheet identifiers
#' @param prefix directory to put CSV files in.  If \code{NULL}, defaults to
#' "metadata/" in the current working directory.
#' @export
update_metadata_via_yaml <- function(fp, prefix=NULL) {
  text <- readChar(fp, file.info(fp)$size)
  google_sheets <- yaml::yaml.load(text)
  update_metadata(google_sheets, prefix)
}


# CSV ---------------------------------------------------------------------


#' Generic CSV loader
#' @export
load_csv <- function(fp, key=NULL) {
  data <- read.csv(
    fp, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(key)) {
    key <- colnames(data)[1]
  }
  if (! is.na(key)) {
    if (nrow(data) != length(unique(data[[key]]))) {
      stop(paste("Non-unique keys in input CSV from column:", key))
    }
    rownames(data) <- data[[key]]
  }
  data
}

#' Save to CSV
#'
#' @param data data frame to save
#' @param fp file path to save to
#'
#' Write a data frame to a CSV file with some common settings.
save_csv <- function(data, fp) {
  write.csv(x = data, file = fp, na = "", row.names = FALSE)
}


# FASTA -------------------------------------------------------------------


#' Load FASTA file
#'
#' Load a FASTA file into a named character vector
#'
#' @param fp FASTA file path
#'
#' @return named character vector
load_fasta <- function(fp) {
  al <- dnar::read.fa(fp)
  txt <- al$seq
  names(txt) <- al$name
  txt
}

#' Save FASTA file
#'
#' Save a FASTA file from a named character vector
#'
#' @param txt named character vector.  names will be used for sequence
#'   description lines.
#' @param fp FASTA file path
#'
#' @return named character vector
save_fasta <- function(txt, fp) {
  dnar::write.fa(names = names(txt), dna = txt, fileName = fp)
}


# Misc --------------------------------------------------------------------


load_blast_table <- function(fp, columns=NULL) {
  if (is.null(columns)) {
    columns <- c(
      "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  }
  data <- read.table(
    fp,
    sep = "\t",
    stringsAsFactors = FALSE,
    header = FALSE,
    col.names = columns)
  data
}

#' Parse gene names from allele names
#' @param txt character vector of immunoglobulin allele names
#' @return character vector of gene names
#' @export
parse_gene <- function(txt) {
  gsub("^I?G?([^\\*]+)\\*.*", "\\1", txt)
}
