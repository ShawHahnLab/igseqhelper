"""
Test SONAR helper functions.

We have some wrappers to help format the input as SONAR expects, and we'll test
to ensure that it logs problems as expected.
"""

import unittest
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from Bio import SeqIO
from igseq import R_PKG_INST as INST
from igseq.sonar import (gather_germline, munge_seqid_for_sonar, LOGGER)

class TestSonar(unittest.TestCase):
    """Test SONAR helper functions."""

    def test_gather_germline(self):
        """Test gathering germline sequences for SONAR.

        This should take in a FASTA file and output a FASTA with sequence IDs
        formatted the way SONAR wants.
        """
        input_fps = [INST/"reference/imgt/IGHV.fasta"]
        with TemporaryDirectory() as outdir:
            output_fp = Path(outdir) / "output.fasta"
            gather_germline(input_fps, output_fp)
            parser = SeqIO.parse(output_fp, "fasta")
            # first and last records
            record1 = next(parser)
            for record in parser:
                record2 = record
        # We should have short simple IDs with just the gene+allele
        self.assertEqual(record1.id, "IGHV1-1*01")
        self.assertEqual(record2.id, "IGHV7-1*01")
        # We should NOT have periods (like IMGT likes to use to ensure
        # standardized positions)
        self.assertNotIn(".", str(record1.seq))

    def test_munge_seqid_for_sonar(self):
        """Test modifying sequence IDs for SONAR."""
        seqids = set()
        seqid = munge_seqid_for_sonar("NW_001121240|IGHV1-1*01|Macaca mulatta", seqids)
        self.assertEqual(seqid, "IGHV1-1*01")
        self.assertEqual(seqids, set(("IGHV1-1*01",)))
        # unrecognized format
        with self.assertLogs(LOGGER, logging.ERROR):
            seqid = munge_seqid_for_sonar("XYZ", seqids)
        # already seen this ID
        with self.assertLogs(LOGGER, logging.ERROR):
            seqid = munge_seqid_for_sonar("IGHV1-1*01", seqids)
