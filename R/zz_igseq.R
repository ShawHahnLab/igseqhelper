# TODO this file starts with zz_ so we can use functions from within the
# package while the package is loading.  This is almostly certainly a dumb
# idea. Fix this.

#' Parse gene names from allele names
#' @param txt character vector of immunoglobulin allele names
#' @return character vector of gene names
#' @export
parse_gene <- function(txt) {
  gsub("^I?G?([^\\*]+)\\*.*", "\\1", txt)
}

#' MiSeq Runs
#'
#' MiSeq run IDs for the raw data.  At first they're split by heavy/light, but
#' in the latest set they're mixed.  I'll have to see if this breaks any of my
#' existing assumptions.
RUNS <- list(
  "dataset1" = c(
  light = "190816_M00281_0522_000000000-CL9YN",
  heavy = "190816_M05588_0227_000000000-CL8PC"),
  "dataset2" = c(
    light = "190827_M00281_0525_000000000-CL479",
    heavy = "190903_M05588_0232_000000000-CLL8M"),
  "dataset3" = c(
    run1 = "191220_M00281_0569_000000000-CN8J9",
    run2 = "191220_M05588_0264_000000000-CN8JR"))

#' Samples
#'
#' Sample (timepoint) identifiers for heavy and light chains
SAMPLES <- list(
  light = c(
    WK12 = "WK12L",
    WK24 = "WK24L",
    WK48 = "WK48L"),
  heavy = c(
    WK12 = "WK12H",
    WK24 = "WK24H",
    WK48 = "WK48H")
)

#' Germline Allele Names
#'
#' Germline allele names that produced the antibodies of interest.
ALLELES <- list(
  heavy = c(
    V = "IGHV4-ABB-S*01_S8200",
    D = "IGHD3-1*01",
    J = "IGHJ2*01"
  ),
  light = c(
    V = "IGLV1-7*01",
    J = "IGLJ2*01"
  )
)

#' Germline Gene Names
#'
#' Gene names for each of the expected germline allele names.
GENES <- list()
GENES$heavy <- sapply(ALLELES$heavy, parse_gene)
GENES$light <- sapply(ALLELES$light, parse_gene)

#' File paths for Reference Sequences
#'
#' See REFERENCES for more information.
PATHS_REFERENCES <- list(
  heavy = list(
    combined = "heavy.fasta",
    mature = "heavy.mature.fasta",
    checkpoints = "heavy.checkpoints.fasta"
  ),
  light = list(
    combined = "light.fasta",
    mature = "light.mature.fasta",
    germline = "light.germline.fasta"
  )
)

#' Reference Sequences
#'
#' Small sets of sequences for each chain for easy reference.  Each is a named
#' character vector, stored within a nested list and grouped by chain.
REFERENCES <- lapply(PATHS_REFERENCES, function(set) {
  lapply(set, function(fp) {
    load_fasta(system.file("reference", fp, package = "igseq"))
  })
})

#' CDR[HL]3 amino acid sequences
#'
#' Amino acid sequences for the CDR3 region of each chain's mature sequences.
#' Here we've picked a representative sequence in cases where the mature
#' antibody sequences aren't identical.
CDR3 <- c(
  # The four for light chain:
  #      ..
  # CQCYDSSVLF (1)
  # CQCYDSNVLF (2)
  # CQCYDTTVLF (3)
  # CQCYDTTVLF (3)
  # This (2) is the only one I see in WK48L.
  light = "CQCYDSNVLF",
  # Approximate (just took one of the four) CDRH3 AA sequence.  In this region
  # they differ by at most one AA: "YFDLW" versus "FFDLW" at end.
  heavy = "CARKGEDFYEDDYGQYFTAGWYFDLW")

#' Antibody Annotations
#'
#' Sets of named position ranges within one antibody sequence, currently just
#' for the heavy chain.
ANNOTATIONS <- list(
  light = list(),
  heavy = list(
    `5695_P12_E6H` = list(
      IMGT = list(
        FR1 = c(1, 75),
        CDR1 = c(76, 105),
        FR2 = c(106, 156),
        CDR2 = c(157, 180),
        FR3 = c(181, 294),
        CDR3 = c(295, 366)
      ),
      # TODO check these
      Germline = list(
        V = c(1, 302),
        # NOTE D will win in this case
        D = c(310, 354),
        J = c(348, 366)
      )
    )
  )
)

#' Google Sheets IDs
#'
#' A list of workbook names and associated sheet data to track.  Each contained
#' list should have one entry named "slug" set to the long alphanumeric string
#' for the workbook, and a second entry named "sheets" containing the integer
#' GID for each sheet to download with names set to the names of the sheets.
GOOGLE_SHEETS <- list(
  metadata = list(
    slug = "2PACX-1vRsKo1rnHyhUcJUdP-fO6ACBPtkNCXTJZHoGmX81hHCFrvlfZ2pMNZBIp7vZIHD9N32GN8Ufn4qQAA4",
    sheets = c(
      samples = "0",
      specimens = "2071037194",
      primers = "1722819559",
      runs = "899418000")
  )
)

#' Fetch metadata from Google Sheets
#'
#' @param google_sheets nested list structure like in
#'   \code{\link{GOOGLE_SHEETS}}
#' @param prefix directory to put CSV files in.  If \code{NULL}, defaults to the
#'   package's inst directory.
#' @export
update_metadata <- function(google_sheets=GOOGLE_SHEETS, prefix=NULL) {
  if (is.null(prefix)) {
    prefix <- system.file(package = "igseq")
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
