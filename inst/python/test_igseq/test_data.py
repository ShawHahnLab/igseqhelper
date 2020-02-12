"""
Test data and metadata functions.
"""

import unittest
import logging
import igseq.data
from test_igseq.util import ROOT


class TestMetadata(unittest.TestCase):
    """Basic tests for metadata loading from spreadsheets."""

    def setUp(self):
        self.path = ROOT / "metadata"
        self.expected = {
            "sequences": {
                "cols": ["Name", "Seq", "Use", "Annotation", "Direction", "Notes"],
                "rows":
                    ["P5_Graft", "P5_Seq", "5PIIA", "P7", "IgG", "IgD", "IgK", "IgL", "RhIgM"] +
                    ["BC_%d" % val for val in range(1, 21)] +
                    ["i7_%d" % val for val in range(1, 21)]
                },
            "runs": {
                "cols": ["Run", "ReverseComplement", "Comments",
                         "URL", "URLR1", "URLR2", "URLI1",
                         "MD5R1", "MD5R2", "MD5I1"],
                "rows": ["200000_M05588_0227_000000000-XYZ12", "200000_M00281_0522_000000000-ABC34"]
                },
            "specimens": {
                "cols": ["Specimen", "Subject", "Timepoint", "CellCount", "CellType", "Comments"],
                "rows": ["specimen1", "specimen2", "specimen3"]
                },
            "samples": {
                "cols": ["Sample", "Run", "Specimen",
                         "BarcodeFwd", "BarcodeRev", "Replicate", "Chain", "Type", "Comments"],
                "rows": ["sample%d" % num for num in range(1, 7)]
                }
            }

    def test_load_sequences(self):
        """Test loading sequences CSV."""
        sequences = igseq.data.load_sequences(self.path / "sequences.csv")
        self.check_metadata(sequences, self.expected["sequences"], "Name")

    def test_load_runs(self):
        """Test loading runs CSV."""
        runs = igseq.data.load_runs(self.path / "runs.csv")
        self.check_metadata(runs, self.expected["runs"], "Run")

    def test_load_specimens(self):
        """Test loading specimens CSV."""
        specimens = igseq.data.load_specimens(self.path / "specimens.csv")
        self.check_metadata(specimens, self.expected["specimens"], "Specimen")

    def test_load_samples_basic(self):
        """Test loading samples CSV, by itself."""
        samples = igseq.data.load_samples(self.path / "samples.csv")
        self.check_metadata(samples, self.expected["samples"], "Sample")

    def check_metadata(self, obs, exp, key):
        """Check that the rows and columns are as expected.
        obs: observed metadata dict structure
        exp: expected metadata row and columns names
        key: expected key for rows (should also be column name)
        """
        # every entry is a row and keys are the row names.
        self.assertEqual(list(obs.keys()), exp["rows"])
        # Every row should have an entry for every column, where keys are
        # column names.  The row name should be present as the expected column
        # name.
        for name, attrs in obs.items():
            self.assertEqual(attrs[key], name)
            self.assertEqual(list(attrs.keys()), exp["cols"])

    def test_load_samples_joined(self):
        """Test loading samples CSV, joined to other metadata.

        This is a more complicated version, where we have all the same sample
        metadata plus additional nested dictionaries filled in from the other
        CSV files.
        """
        specimens = igseq.data.load_specimens(self.path / "specimens.csv")
        runs = igseq.data.load_runs(self.path / "runs.csv")
        sequences = igseq.data.load_sequences(self.path / "sequences.csv")
        samples = igseq.data.load_samples(self.path / "samples.csv", specimens, runs, sequences)
        # Check the same sort of rows/columns things as above, but also the nested metadata
        nested_items = {
            "SpecimenAttrs": "specimens",
            "BarcodeFwdAttrs": "sequences",
            "BarcodeRevAttrs": "sequences",
            "RunAttrs": "runs"}
        self.assertEqual(list(samples.keys()), self.expected["samples"]["rows"])
        for name, attrs in samples.items():
            self.assertEqual(name, attrs["Sample"])
            keys = attrs.keys()
            # Each of the expected nested metadata entries should be present,
            # with the same keys on those as expected for their standalone
            # spreadsheets.
            for key, sheet in nested_items.items():
                self.assertIn(key, keys)
                self.assertEqual(list(attrs[key].keys()), self.expected[sheet]["cols"])
            # Excluding those nested keys, the sample attributes should still match up as before.
            keys = [key for key in keys if key not in nested_items.keys()]
            self.assertEqual(keys, self.expected["samples"]["cols"])


class TestMetadataDuplicates(TestMetadata):
    """Test metadata spreadsheets that have duplicated entries.

    This shouldn't be allowed, since each type of thing should have exactly one row per thing.
    There are a few additional types of duplicates to test for: for sequences,
    neither the name nor the sequence itself should occur more than once.  For
    samples, the combination of barcodes and run ID should be unique.
    """

    def setUp(self):
        self.path = ROOT / "metadata_duplicates"

    def test_load_sequences(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(ValueError):
                igseq.data.load_sequences(self.path / "sequences.csv")
        # But also, we shouldn't have any duplicated sequences.
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(ValueError):
                igseq.data.load_sequences(self.path / "sequences_dup_seq.csv")

    def test_load_runs(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(ValueError):
                igseq.data.load_runs(self.path / "runs.csv")

    def test_load_specimens(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(ValueError):
                igseq.data.load_specimens(self.path / "specimens.csv")

    def test_load_samples_basic(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(ValueError):
                igseq.data.load_samples(self.path / "samples.csv")
        # But also, we shouldn't have any duplicated combinations of
        # BarcodeFwd, BarcodeRev, and Run.  In this case the supposed sample7
        # is indistinguishable from sample1 because  it's in the same run with
        # the same barcodes.
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(ValueError):
                igseq.data.load_samples(self.path / "samples_dup_barcodes.csv")

    def test_load_samples_joined(self):
        self.skipTest("not yet implemented")
