#' Load primer sequence metadata
#'
#' @param fp file path to sample attributes CSV (defaults to the file built into
#'   the package).
#'
#' @return data frame of primer names and sequences
load_primers <- function(fp=NULL) {
  if (is.null(fp)) {
    fp <- system.file("metadata", "primers.csv", package = "igseq")
  }
  primers <- read.csv(fp, stringsAsFactors = FALSE)
  primers
}

#' Load sample attributes
#'
#' @param primers data frame of primer metadata to include
#' @param fp file path to sample attributes CSV (defaults to the file built into
#'   the package)
#'
#' @return data frame of sample attributes include primers used
load_s_attrs <- function(primers, fp=NULL) {
  if (is.null(fp)) {
    fp <- system.file("metadata", "samples.csv", package = "igseq")
  }
  s_attrs <- read.csv(fp, stringsAsFactors = FALSE)
  s_attrs$PrimerFwdSeq <- primers$PrimerSeq[match(
    s_attrs$PrimerFwd,
    primers$PrimerName)]
  s_attrs$PrimerRevSeq <- primers$PrimerSeq[match(
    s_attrs$PrimerRev,
    primers$PrimerName)]
  s_attrs$BarcodeFwd <- sub("^(N+[ACTG]{8}).*$", "\\1", s_attrs$PrimerFwdSeq)
  s_attrs$BarcodeRev <- substr(s_attrs$PrimerRevSeq, 25, 25 + 7)
  s_attrs$BarcodePair <- with(s_attrs, factor(paste(BarcodeFwd, BarcodeRev)))
  s_attrs
}

load_specimens <- function(fp=NULL) {
  if (is.null(fp)) {
    fp <- system.file("metadata", "specimens.csv", package = "igseq")
  }
  specimens <- read.csv(fp, stringsAsFactors = FALSE)
  specimens
}

load_counts_table <- function(fp="counts.csv", ...) {
  read.csv(fp, stringsAsFactors = FALSE, ...)
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
