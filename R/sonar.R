# The allele sequences used in the SONAR analysis.
load_germline <- function(dp_root=".") {
  dp <- file.path(dp_root, "sonar-analysis")
  list(
    light = list(
      V = load_fasta(file.path(dp, "light/germline.V.fasta")),
      J = load_fasta(file.path(dp, "light/germline.J.fasta"))),
    heavy = list(
      V = load_fasta(file.path(dp, "heavy/germline.V.fasta")),
      D = load_fasta(file.path(dp, "heavy/germline.D.fasta")),
      J = load_fasta(file.path(dp, "heavy/germline.J.fasta")))
  )
}
