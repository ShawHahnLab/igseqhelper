# IMGT's reference sequences for Rhesus Macaque IgG

URL_IMGT <- "http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Macaca_mulatta/IG/"

fetch_imgt_ref <- function(region) {
  fa <- dnar::read.fa(file.path(URL_IMGT, paste0(region, ".fasta")))
}

download_all_imgt_refs <- function(dir_out=NULL) {
  if (is.null(dir_out)) {
    dir_out <- file.path(system.file("reference", package = "igseq"), "imgt")
  }
  if (! dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }
  for (region in c("IGHV", "IGHD", "IGHJ", "IGLV", "IGLJ")) {
    fp_out <- file.path(dir_out, paste0(region, ".fasta"))
    fa <- fetch_imgt_ref(region)
    dnar::write.fa(fa$name, toupper(fa$seq), fp_out)
  }
}
