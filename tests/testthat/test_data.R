context("Test data and metadata functions")

testdir <- function() {
  # TODO should probably reorganize now that I'm using this in both packages.
  system.file("python", "test_igseq", package="igseq")
}

expected_metadata <- function() {
  expected <- list(
    sequences = list(
      cols = c("Name", "Seq", "Use", "Annotation", "Direction", "Notes"),
      rows = c(
        c("P5_Graft", "P5_Seq", "5PIIA", "P7",
          "IgG", "IgD", "IgK", "IgL", "RhIgM"),
          paste0("BC_", 1:20), paste0("i7_", 1:20))
    ),
    runs = list(
      cols = c(
        "Run", "ReverseComplement", "Comments", "URL", "URLR1", "URLR2",
        "URLI1", "MD5R1", "MD5R2", "MD5I1"),
      rows = c(
        "200000_M05588_0227_000000000-XYZ12",
        "200000_M00281_0522_000000000-ABC34")
    ),
    specimens = list(
      cols = c(
        "Specimen", "Subject", "Timepoint",
        "CellCount", "CellType", "Comments"),
      rows = c("specimen1", "specimen2", "specimen3")
    ),
    samples = list(
      cols = c(
        "Sample", "Run", "Specimen", "BarcodeFwd", "BarcodeRev", "Replicate",
        "Chain", "Type", "Comments"),
      rows = paste0("sample", 1:6)
    )
  )
  expected
}


# Named Metadata ----------------------------------------------------------


test_that("load_sequences can load a sequence metadata CSV", {
  sequences <- load_sequences(file.path(testdir(), "metadata", "sequences.csv"))
  expected <- expected_metadata()$sequences
  expect_equal(rownames(sequences), expected$rows)
  expect_equal(colnames(sequences), expected$cols)
})

test_that("load_runs can load a run metadata CSV", {
  runs <- load_runs(file.path(testdir(), "metadata", "runs.csv"))
  expected <- expected_metadata()$runs
  expect_equal(rownames(runs), expected$rows)
  expect_equal(colnames(runs), expected$cols)
})

test_that("load_specimens can load a specimen metadata CSV", {
  specimens <- load_specimens(file.path(testdir(), "metadata", "specimens.csv"))
  expected <- expected_metadata()$specimens
  expect_equal(rownames(specimens), expected$rows)
  expect_equal(colnames(specimens), expected$cols)
})

test_that("load_samples can load a basic sample metadata CSV", {
  samples <- load_samples(file.path(testdir(), "metadata", "samples.csv"))
  expected <- expected_metadata()$samples
  expect_equal(rownames(samples), expected$rows)
  expect_equal(colnames(samples), expected$cols)
})

test_that("load_samples can load a joined sample metadata CSV", {
  specimens <- load_specimens(file.path(testdir(), "metadata", "specimens.csv"))
  runs <- load_runs(file.path(testdir(), "metadata", "runs.csv"))
  sequences <- load_sequences(file.path(testdir(), "metadata", "sequences.csv"))
  samples <- load_samples(
    file.path(testdir(), "metadata", "samples.csv"), specimens, runs, sequences)
  expected <- expected_metadata()$samples
  expected$cols <- c(
    expected$cols,
    paste0("Specimen", colnames(specimens)),
    paste0("Run", colnames(runs)),
    paste0("BarcodeFwd", colnames(sequences)),
    paste0("BarcodeRev", colnames(sequences)))
  expected$cols <- expected$cols[
    -match(
      c("SpecimenSpecimen", "RunRun", "BarcodeFwdName", "BarcodeRevName"),
      expected$cols)]
  expect_equal(rownames(samples), expected$rows)
  expect_equal(colnames(samples), expected$cols)
})
