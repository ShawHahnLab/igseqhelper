"""
Test data and metadata functions.
"""

import re
import unittest
import builtins
import gzip
import logging
from tempfile import NamedTemporaryFile
import igseq.data
from igseq.data import MetadataError
from test_igseq.util import ROOT


def cat(file_path, opener=builtins.open, **kwargs):
    """Read in contents of file."""
    with opener(file_path, **kwargs) as f_in:
        return f_in.read()

def zcat(file_path, **kwargs):
    """Read in contents of gzipped text file."""
    kwargs["mode"] = kwargs.get("mode", "rt")
    return cat(file_path, gzip.open, **kwargs)

class TestMetadataBase(unittest.TestCase):
    """Base test class for metadata loading from spreadsheets.

    These tests will run for any metadata scenario.
    """

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
                },
            "antibody_isolates": {
                "cols": ["AntibodyIsolate", "AntibodyLineage", "HeavySeq", "LightSeq", "Comments"],
                "rows": ["ISO_S1_1_A", "ISO_S1_1_B", "ISO_S3_1_A", "ISO_S4_1_A"]
                },
            "antibody_lineages": {
                "cols": [
                    "AntibodyLineage", "Subject", "HeavyConsensus", "LightConsensus",
                    "VH", "DH", "JH", "VL", "JL", "Comments"],
                "rows": ["ISO_S1_1", "ISO_S3_1", "ISO_S4_1"]
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

    def test_load_antibody_lineages(self):
        """Test loading antibody lineages CSV."""
        lineages = igseq.data.load_antibody_lineages(self.path / "antibody_lineages.csv")
        self.check_metadata(lineages, self.expected["antibody_lineages"], "AntibodyLineage")

    def test_load_antibody_isolates_basic(self):
        """Test loading antibody isolates CSV, by itself."""
        isolates = igseq.data.load_antibody_isolates(self.path / "antibody_isolates.csv")
        self.check_metadata(isolates, self.expected["antibody_isolates"], "AntibodyIsolate")

    def test_load_antibody_isolates_joined(self):
        """Test loading antibody isolates CSV, joined to lineage metadata.

        See a similar test (just more intricate) in test_load_samples_joined.
        """
        lineages = igseq.data.load_antibody_lineages(self.path / "antibody_lineages.csv")
        isolates = igseq.data.load_antibody_isolates(self.path / "antibody_isolates.csv", lineages)
        nested_items = {"AntibodyLineageAttrs": "antibody_lineages"}
        self.assertEqual(list(isolates.keys()), self.expected["antibody_isolates"]["rows"])
        for name, attrs in isolates.items():
            self.assertEqual(name, attrs["AntibodyIsolate"])
            keys = attrs.keys()
            for key, sheet in nested_items.items():
                self.assertIn(key, keys)
                self.assertEqual(list(attrs[key].keys()), self.expected[sheet]["cols"])
            keys = [key for key in keys if key not in nested_items.keys()]
            self.assertEqual(keys, self.expected["antibody_isolates"]["cols"])

    def test_load_samples_basic(self):
        """Test loading samples CSV, by itself."""
        samples = igseq.data.load_samples(self.path / "samples.csv")
        self.check_metadata(samples, self.expected["samples"], "Sample")

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

    def test_load_csv(self):
        """Test generic CSV loading."""
        # This should behave similarly to the load_sequences wrapper, except
        # for the extra sequence content processing.
        sequences = igseq.data.load_csv(self.path / "sequences.csv")
        for _, seq_data in sequences.items():
            seq_data["Seq"] = re.sub(" ", "", seq_data["Seq"])
        self.check_metadata(sequences, self.expected["sequences"], "Name")
        # Or, we can use a different key, so long as it's unique.
        sequences = igseq.data.load_csv(self.path / "sequences.csv", "Seq")
        # Non-unique should give an error.
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_csv(self.path / "sequences.csv", "Use")

    def test_get_data(self):
        """Test getting data from local disk or URLs."""
        self.skipTest("not yet implemented.")

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


class TestMetadata(TestMetadataBase):
    """Basic tests for metadata loading from spreadsheets.

    This does some additioal checks for correctly-formatted metadata sheets.
    """

    def test_get_samples_per_run(self):
        """Test making a dictionary of run IDs to sample names."""
        samples = igseq.data.load_samples(self.path / "samples.csv")
        samples_per_run = igseq.data.get_samples_per_run(samples)
        # If we tally up pairs of run Id and sample name we should get the same
        # list either way
        pairs_exp = [(entry["Run"], entry["Sample"]) for entry in samples.values()]
        pairs_exp = sorted(pairs_exp)
        pairs_obs = []
        for run_id, samp_list in samples_per_run.items():
            for samp in samp_list:
                pairs_obs.append((run_id, samp))
        pairs_obs = sorted(pairs_obs)
        self.assertEqual(pairs_obs, pairs_exp)


class TestMetadataDuplicates(TestMetadataBase):
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
            with self.assertRaises(MetadataError):
                igseq.data.load_sequences(self.path / "sequences.csv")
        # But also, we shouldn't have any duplicated sequences.
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_sequences(self.path / "sequences_dup_seq.csv")

    def test_load_runs(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_runs(self.path / "runs.csv")

    def test_load_specimens(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_specimens(self.path / "specimens.csv")

    def test_load_antibody_lineages(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_antibody_lineages(self.path / "antibody_lineages.csv")

    def test_load_antibody_isolates_basic(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_antibody_isolates(self.path / "antibody_isolates.csv")

    def test_load_antibody_isolates_joined(self):
        self.skipTest("not yet implemented")

    def test_load_samples_basic(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_samples(self.path / "samples.csv")
        # But also, we shouldn't have any duplicated combinations of
        # BarcodeFwd, BarcodeRev, and Run.  In this case the supposed sample7
        # is indistinguishable from sample1 because  it's in the same run with
        # the same barcodes.
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_samples(self.path / "samples_dup_barcodes.csv")

    def test_load_samples_joined(self):
        self.skipTest("not yet implemented")

    def test_load_csv(self):
        with self.assertLogs(level=logging.CRITICAL):
            with self.assertRaises(MetadataError):
                igseq.data.load_csv(self.path / "sequences.csv")


class TestData(unittest.TestCase):
    """Basic tests for data helper functions."""

    def setUp(self):
        self.path = ROOT / "data"

    def test_get_data(self):
        """Test getting run data from local disk or URL."""
        self.skipTest("not yet implemented")

    def test_md5(self):
        """Test MD5 checksum on a file."""
        # /dev/null is an empty file when reading.
        self.assertEqual(
            igseq.data.md5("/dev/null"),
            "d41d8cd98f00b204e9800998ecf8427e")

    def test_chunk_fqgz(self):
        """Test divvying up a fastq.gz into a fixed number of files."""

        path_in = self.path / "stub.fastq.gz"
        # With 1:1 you get the same thing out as you put in
        with NamedTemporaryFile() as path_out:
            igseq.data.chunk_fqgz([path_in], [path_out.name])
            self.assertEqual(zcat(path_in), zcat(path_out.name))
        # With 1:N it's split up and will leave empty files for any extras
        paths_out = [NamedTemporaryFile() for _ in range(6)]
        igseq.data.chunk_fqgz([path_in], paths_out)
        for path_out in paths_out:
            path_out.flush()
        self.assertEqual(
            zcat(path_in),
            "".join([zcat(path_out.name) for path_out in paths_out]))
        self.assertEqual("", zcat(paths_out[-1].name))
        for path_out in paths_out:
            path_out.close()

    def test_amplicon_files(self):
        """Test filename helper for files per specimen per target chain type."""
        # I think we can probably remove this outright eventually, and just use
        # snakemake's optional function argument to expand() for the same
        # functionality.
        self.skipTest("not yet implemented")
