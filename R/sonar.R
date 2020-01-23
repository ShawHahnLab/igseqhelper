load_module_1 <- function(dp_root=".") {
  module1 <- list()
  module1$counts <- list()
  module1$rearrs <- list()
  for (chain in names(SAMPLES)) {
    for (samp in names(SAMPLES[[chain]])) {
      fp <- file.path(
        dp_root,
        "sonar-analysis",
        chain,
        SAMPLES[[chain]][[samp]],
        "output/tables",
        paste0(SAMPLES[[chain]][[samp]], "_rearrangements.tsv"))
      cat("Loading ", fp, "\n", sep = "", file = 2)
      rearr <- airr::read_rearrangement(fp)
      idxl_good <- rearr$status == "good"
      idxl_v_match <- grepl(ALLELES[[chain]][["V"]], rearr$v_call, fixed = TRUE)
      idxl_j_match <- grepl(ALLELES[[chain]][["J"]], rearr$j_call, fixed = TRUE)
      x <- rearr[idxl_good & idxl_v_match & idxl_j_match, ]
      x$Chain = chain
      x$Sample = samp
      module1$rearrs <- c(module1$rearrs, list(x))
      module1$counts <- c(module1$counts, list(data.frame(
        Chain = chain,
        Sample = samp,
        Total = sum(rearr$duplicate_count),
        Dedup = nrow(rearr),
        Good = sum(idxl_good),
        VJMatch = sum(idxl_good & idxl_v_match & idxl_j_match),
        stringsAsFactors = FALSE)))
    }
  }
  module1$rearrs <- do.call(rbind, module1$rearrs)
  module1$counts <- do.call(rbind, module1$counts)
  module1
}

load_module_2 <- function(dp_root=".") {
  module2 <- list()
  module2$islands <- list()
  for (chain in names(SAMPLES)) {
    module2$islands[[chain]] <- lapply(SAMPLES[[chain]], function(samp) {
      load_fasta(file.path(
        dp_root,
        "sonar-analysis",
        chain,
        samp,
        "output/sequences/nucleotide",
        paste0(samp, "_islandSeqs.fa")))
    })
  }
  
  module2$lineages <- lapply(names(SAMPLES), function(chain) {
    lapply(SAMPLES[[chain]], function(samp) {
      fp <- file.path(
        dp_root,
        "sonar-analysis",
        chain,
        samp,
        "output/tables",
        paste0(samp, "_lineages.txt"))
      if (file.exists(fp)) {
        read.table(fp, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      }
    })
  })
  names(module2$lineages) <- names(SAMPLES)
  
  module2$lineages_long <- do.call(rbind, lapply(
    names(module2$lineages), function(chain) {
    do.call(rbind, lapply(names(module2$lineages[[chain]]), function(samp) {
      chunk <- module2$lineages[[chain]][[samp]]
      if (! is.null(chunk)) {
        chunk$Chain <- factor(chain, levels = names(SAMPLES))
        chunk$Sample <- factor(samp, levels = names(SAMPLES[[1]]))
        chunk$v_match <- chunk$v_call %in% GENESETS[[chain]][["V"]]
        chunk$j_match <- chunk$j_call %in% GENESETS[[chain]][["J"]]
        chunk$junction_dist <- levenR::leven(chunk$junction_aa, CDR3[[chain]])
        chunk
      }
    }))
  }))
  module2
}

load_module_3 <- function(dp_root=".") {
  module3 <- list()
  for (chain in names(SAMPLES)) {
    module3$collected[[chain]] <- load_fasta(file.path(
      dp_root,
      "sonar-analysis",
      chain,
      "longitudinal/output/sequences/nucleotide/longitudinal-collected.fa"))
    # There's one sequence here for every node (including tips) present in the
    # tree, but the order isn't consistent.  see phylo.R in igseq for helpers.
    module3$inferred[[chain]] <- load_fasta(file.path(
      dp_root,
      "sonar-analysis",
      chain,
      "longitudinal/output/sequences/nucleotide",
      "longitudinal_inferredAncestors.fa"))
    module3$trees[[chain]] <- read.tree(file.path(
      dp_root,
      "sonar-analysis",
      chain,
      "longitudinal/output/longitudinal_igphyml.tree"))
  }
  module3
}

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

# Some of the reference allele sequences are very close to our alleles even 
# though they're supposedly from different genes.  The lineage tables only
# report the first match per gene so we'll make a note of the similar alleles to
# track those too.
make_allele_sets <- function(genes, alleles, germline, maxdist=4) {
  allelesets <- list()
  for (chain in c("light", "heavy")) {
    allelesets[[chain]] <- list()
    for (part in names(genes[[chain]])) {
      allele_name <- alleles[[chain]][[part]]
      seqs <- gsub("\\.", "", germline[[chain]][[part]])
      dists <- levenR::leven(seqs, seqs[allele_name])
      allelesets[[chain]][[part]] <- names(seqs[dists <= maxdist])
    }
  }
  allelesets
}

# given a tree and a named vector of inferred sequences from SONAR/IgPhyML,
# match up the inferred sequences with the tree's nodes.  There is one that's
# inherently missing for the root of the tree itself; this will always be NA
# (but is left in place so the vector indexing matches up).
match_inferred_ancestors <- function(tree, seqs) {
  # Character vector version of the clades for each node in the tree
  tree_clades <- sapply(tree_get_clades(tree), paste, collapse = ",")
  # These look like they might be sorted already but I'm not going to chance it
  inf_clades <- sapply(strsplit(names(seqs), ";"), function(vec) {
    paste(sort(strsplit(vec[2], ",")[[1]]), collapse=",")
    })
  # These are the inferred sequences for every node (following the tree's
  # numbering)
  module3$inferred$heavy[match(tree_clades, inf_clades)]
}
