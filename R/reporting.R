prep_counts_table <- function(fp="counts.csv") {
  chains <- c(light="light", heavy="heavy")
  runs_by_chain <- lapply(chains, function(chain) sapply(RUNS, '[', chain))
  data <- load_counts_table(fp)
  data$Timepoint <- label_timepoint(data$Filename)
  data$Run <- label_run(data$Filename)
  data$Chain <- label_chain(data, runs_by_chain)
  data$Stage <- label_stage(data$Filename)
  data
}

label_timepoint <- function(txt) {
  re_factor(txt, "WK..[HL]")
}

label_run <- function(txt) {
  re_factor(txt, "[0-9]{6}_[A-Z0-9]+_[0-9]{4}_0{9}-[A-Z0-9]{5}")
}

label_chain <- function(data, runs_by_chain) {
  with(data,
  factor(
    ((Run %in% runs_by_chain[["light"]]) %in% TRUE ) |
    (is.na(Run) & Timepoint %in% SAMPLES$light),
    labels = c("Light", "Heavy"),
    levels = c(TRUE, FALSE)))
}

re_factor <- function(txt, pat) {
  factor(sapply(regmatches(txt, gregexpr(pat, txt)), `[`, 1))
}

label_stage <- function(txt) {
  prefix <- sub("/.*", "", txt)
  stage <- factor(
    prefix,
    levels = c(
      "0-data", "1-primerfilt", "2-primerfilt-peared", "4-demux", "6-derep", "7-blast"),
    labels = c(
      "Original", "Filtered for Primer", "Reads Merged", "Samples Assigned", "Deduplicated", "Matched to References"))
  stage
}

summarize_counts <- function(counts_raw) {
  counts <- reshape2::dcast(
    data = counts_raw,
    formula = Stage ~ Chain + Timepoint,
    fun.aggregate = sum,
    value.var = "NumSequences")
  colnames(counts)[match("Light_NA", colnames(counts))] <- "Light_Total"
  colnames(counts)[match("Heavy_NA", colnames(counts))] <- "Heavy_Total"
  stages_post_demux <- c("Samples Assigned", "Deduplicated", "Matched to References")
  stages_post_demux <- stages_post_demux[stages_post_demux %in% counts$Stage]
  for (stage in stages_post_demux) {
    idxrow <- match(stage, counts$Stage)
    counts[idxrow, "Light_Total"] <- sum(subset(counts_raw, Stage == stage & Chain == "Light")$NumSequences)
    counts[idxrow, "Heavy_Total"] <- sum(subset(counts_raw, Stage == stage & Chain == "Heavy")$NumSequences)
  }
  rownames(counts) <- counts$Stage
  # drop any completely zero columns
  idxl_remove <- sapply(counts, function(vec) if (is.numeric(vec)) {sum(vec, na.rm = TRUE) == 0} else {FALSE})
  counts <- counts[, ! idxl_remove]
  counts
}