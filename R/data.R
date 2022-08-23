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
