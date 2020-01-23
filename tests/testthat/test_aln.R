context("Test alignment functions")

make_example_aln <- function() {
  c(
    A = paste0(c("GGCT-ACGA", rep("ACTT", 10), "CGGCACGTAC"), collapse = ""),
    B = paste0(c("GGCTTACGA--------", rep("ACTT", 8), "CG-CACGTAC"),
               collapse = ""),
    C = paste0(c("GCCTTACGA----", rep("ACTT", 9), "CGGCACTTAC"), collapse = "")
  )
}


# make_aln_mat ------------------------------------------------------------


test_that("make_aln_mat makes matrix version", {
  aln_vec <- make_example_aln()
  aln_mat <- make_aln_mat(aln_vec)
  expect_equal(rownames(aln_mat), names(aln_vec))
  expect_equal(ncol(aln_mat), max(nchar(aln_vec)))
  expect_equal(aln_mat[, 14], c(A = "A", B = "-", C = "A"))
})

test_that("make_aln_mat handles length mismatch", {
  aln_vec <- make_example_aln()
  aln_vec[1] <- paste0(aln_vec[1], "A")
  expect_warning(aln_mat <- make_aln_mat(aln_vec), "length mismatch")
  expect_equal(rownames(aln_mat), names(aln_vec))
  expect_equal(ncol(aln_mat), max(nchar(aln_vec)))
  expect_equal(aln_mat[, 14], c(A = "A", B = "-", C = "A"))
  expect_equal(aln_mat[, max(nchar(aln_vec))], c(A = "A", B = NA, C = NA))
})

test_that("make_aln_df handles duplicated names", {
  aln_vec <- make_example_aln()
  names(aln_vec) <- c("A", "B", "B")
  expect_warning(aln_mat <- make_aln_mat(aln_vec), "duplicate names")
  expect_equal(nrow(aln_mat), length(aln_vec))
  expect_equal(ncol(aln_mat), max(nchar(aln_vec)))
})


# make_aln_df -------------------------------------------------------------


check_aln_df <- function(aln_df, shape = c(177, 6)) {
  expect_equal(dim(aln_df), shape)
  expect_equal(
    sapply(aln_df, class),
    c(SeqNum = "integer",
      Position = "integer",
      Base = "factor",
      MatchesReference = "logical",
      Name = "factor",
      BaseMasked = "character"))
}

test_that("make_aln_df makes data frame version", {
  aln_vec <- make_example_aln()
  aln_df <- make_aln_df(aln_vec)
  check_aln_df(aln_df)
})

test_that("make_aln_df handles nameless version", {
  aln_vec <- make_example_aln()
  names(aln_vec) <- NULL
  aln_df <- make_aln_df(aln_vec)
  check_aln_df(aln_df)
})

test_that("make_aln_df handles length mismatch", {
  aln_vec <- make_example_aln()
  aln_vec[1] <- paste0(aln_vec[1], "A")
  expect_warning(aln_df <- make_aln_df(aln_vec), "length mismatch")
  check_aln_df(aln_df, shape = c(180, 6))
  chunk <- subset(aln_df, Position == 60)
  expect_equal(chunk$Base,
               factor(c("A", NA, NA), levels = c("-", "A", "C", "G", "T")))
})

test_that("make_aln_df handles duplicated names", {
  aln_vec <- make_example_aln()
  names(aln_vec) <- c("A", "B", "B")
  expect_warning(aln_df <- make_aln_df(aln_vec), "duplicate names")
  check_aln_df(aln_df)
})

test_that("make_aln_df respects name order", {
  aln_vec <- rev(make_example_aln())
  aln_df <- make_aln_df(aln_vec)
  check_aln_df(aln_df)
  expect_equal(levels(aln_df$Name), c("C", "B", "A"))
})


# plot_aln ----------------------------------------------------------------


test_that("plot_aln plots an alignment", {
  aln_vec <- make_example_aln()
  plt <- plot_aln(aln_vec)
  vdiffr::expect_doppelganger("Alignment Plot", plt)
})

test_that("plot_aln handles length mismatch", {
  aln_vec <- make_example_aln()
  aln_vec[1] <- paste0(aln_vec[1], "A")
  expect_warning(plt <- plot_aln(aln_vec), "length mismatch")
  expect_warning(
    vdiffr::expect_doppelganger("Alignment Plot - Length Mismatch", plt),
    "missing values")
})


test_that("plot_aln handles duplicate names", {
  aln_vec <- make_example_aln()
  names(aln_vec) <- c("A", "B", "B")
  expect_warning(plt <- plot_aln(aln_vec), "duplicate names")
  vdiffr::expect_doppelganger("Alignment Plot - Duplicate Names", plt)
})

test_that("plot_aln handles nameless version", {
  aln_vec <- make_example_aln()
  names(aln_vec) <- NULL
  plt <- plot_aln(aln_vec)
  vdiffr::expect_doppelganger("Alignment Plot - Nameless", plt)
})

test_that("plot_aln respects name order", {
  aln_vec <- rev(make_example_aln())
  plt <- plot_aln(aln_vec)
  vdiffr::expect_doppelganger("Alignment Plot - Reversed Names", plt)
})

test_that("plot_aln can use data frame", {
  aln_vec <- make_example_aln()
  plt1 <- plot_aln(aln_vec)
  aln_df <- make_aln_df(aln_vec)
  plt2 <- plot_aln(aln_df = aln_df)
  expect_equal(
    plt1[-match("plot_env", names(plt1))],
    plt2[-match("plot_env", names(plt2))])
})

test_that("plot_aln can facet sequences", {
  aln_vec <- make_example_aln()
  aln_df <- make_aln_df(aln_vec)
  aln_df$Facet <- factor(
    aln_df$Name == "A",
    labels = c("Alternate", "Reference"),
    levels = c(FALSE, TRUE))
  plt <- plot_aln(aln_df = aln_df)
  vdiffr::expect_doppelganger("Alignment Plot - Faceted", plt)
})


# plot_aln styles ---------------------------------------------------------

test_that("plot_aln can style as vsref", {
  aln_vec <- make_example_aln()
  plt <- plot_aln(aln_vec, style = "vsref")
})

test_that("plot_aln can style as dnaplotr", {
  aln_vec <- make_example_aln()
  plt <- plot_aln(aln_vec, style = "dnaplotr")
})


test_that("plot_aln stops with unknown style", {
  aln_vec <- make_example_aln()
  expect_error(plt <- plot_aln(aln_vec, style = "wrong"), "not recognized")
})
