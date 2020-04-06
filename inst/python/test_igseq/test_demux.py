"""
Test demultiplexing functions.
"""
import unittest
import gzip
from pathlib import Path
from io import StringIO
from tempfile import TemporaryDirectory
from Bio import SeqIO
from igseq.demux import (assign_barcode_fwd, assign_barcode_rev, demux, DemuxError)
from test_igseq.util import DATA


class TestDemux(unittest.TestCase):
    """Test demultiplexing.

    Here we have a single read pair which is a perfect match to sample2.
    """

    def setUp(self):
        self.record = "perfect"
        self.samples = {
            "sample1": {
                "Sample": "sample1",
                "BarcodeFwdAttrs": {"Seq": "NNNNNAACTCTAA"},
                "BarcodeRevAttrs": {"Seq": "CGGGAAGA"}},
            "sample2": {
                "Sample": "sample2",
                "BarcodeFwdAttrs": {"Seq": "NNNNNNAAGGCCCT"},
                "BarcodeRevAttrs": {"Seq": "TGCAACAA"}},
            "sample3": {
                "Sample": "sample3",
                "BarcodeFwdAttrs": {"Seq": "NNNNNNNAATATGTC"},
                "BarcodeRevAttrs": {"Seq": "CATCTGTG"}}}
        prefix = ["F", "M05588:232:000000000-CLL8M:1:1101:19922:1593", "CGGGCTAAGGCCCT"]
        stream_fwd = "\n".join([
            "\t".join(prefix + ["NNNNNNAAGGCCCT", "8.0", "*"]),
            "\t".join(prefix + ["NNNNNNNAATATGTC", "4.0", ""]),
            "\t".join(prefix + ["NNNNNAACTCTAA", "4.0", ""])]) + "\n"
        prefix = ["R", "M05588:232:000000000-CLL8M:1:1101:19922:1593", "TGCAACAA"]
        stream_rev = "\n".join([
            "\t".join(prefix + ["TGCAACAA", "8.0", "*"]),
            "\t".join(prefix + ["CGGGAAGA", "4.0", ""]),
            "\t".join(prefix + ["CATCTGTG", "3.0", ""])]) + "\n"
        names = sorted(self.samples.keys()) + ["unassigned"]
        self.expected = {
            "sample": "sample2",
            "demuxed_data":  {
                "sample1": [],
                "sample2": [[
                    "M05588:232:000000000-CLL8M:1:1101:19922:1593",
                    "AAGCAGTGGTATCAACGCAGAGTACATGGGTTTCCCCCCCCTCTCC" # barcode removed
                    "TCGTGTTCCTCATTTTCCTCCCCTTCAGTCTTCGTGGTTTTATCCTCCCACTTTTCTTCT"
                    "TCGTCACCCAGCGCCCATCCCTCTTCCCCCTCGCTTGTTTCTATCTCGTCTTCCGTCTTC"
                    "TCTTTTCCTCTTCCTCTTCCTTCCCCCTCTTCCCTTCCTTCTTCCCCCTTTCTTTTTTCT"
                    "TTTTCGTCTTCCGTCTTCTTCTTTTCCTCTCCCCCCCTCTCTCTTCTTCTCTCTCTCTCT"
                    "TTCTTCCTC",
                    "M05588:232:000000000-CLL8M:1:1101:19922:1593",
                    "CCCTTCTGCTTCCTTGGCTGGTATTTCCTCCATCACTTCTCTGGGGGACAATGACCATCA"
                    "CGAGGGGGGGGGCGGGCCTCCCCTGTACTCTTCCTTTATACCACTCCCTCGGGCCTTTGC"
                    "CCCGGATCCGTAGTGCCTCCCTTTTGGCAAGTGTTTTGATCTTGCTTGTCGCCGTTTCCT"
                    "TTAAAAAAAAAATTTTCCCCCTTTCTTCCTTCCTTCGTCTTCTTTCTGTCTCTTCTCTCT"
                    "CTTTGTCGCCGTCTACTTTACATCCTACTACGTCAATTGCTCTCTTTTCTCTTTCTTCTC"
                    "TTCTTTGTC",
                    "M05588:232:000000000-CLL8M:1:1101:19922:1593",
                    "TTGTTGCA"]],
                "sample3": []},
            "files": sorted(
                ["%s.R1.fastq.gz" % sample for sample in names] +
                ["%s.R2.fastq.gz" % sample for sample in names] +
                ["%s.I1.fastq.gz" % sample for sample in names]),
            "stream_fwd": stream_fwd,
            "stream_rev": stream_rev
            }

    def test_assign_barcode_fwd(self):
        """Try assigning a barcode to a single R1 read.

        We should get the second barcode in the list returned as the assigned
        value, and should see TSV logged to the supplied stream for each
        barcode considered.
        """
        keys = sorted(self.samples.keys())
        barcodes_fwd = [self.samples[s]["BarcodeFwdAttrs"]["Seq"] for s in keys]
        stream = StringIO()
        record = self.load_input()
        assigned = assign_barcode_fwd(record["R1"][0], barcodes_fwd, send_stats=stream)
        for samp_name, sample in self.samples.items():
            if samp_name == self.expected["sample"]:
                barcode_expected = sample["BarcodeFwdAttrs"]["Seq"]
        self.assertEqual(assigned, barcode_expected)
        self.assertEqual(stream.getvalue(), self.expected["stream_fwd"])

    def test_assign_barcode_rev(self):
        """Try assigning a barcode to a single I1 read.

        We should get the second barcode in the list returned as the assigned
        value, and should see TSV logged to the supplied stream for each
        barcode considered.
        """
        keys = sorted(self.samples.keys())
        barcodes_rev = [self.samples[s]["BarcodeRevAttrs"]["Seq"] for s in keys]
        stream = StringIO()
        record = self.load_input()
        assigned = assign_barcode_rev(record["I1"][0], barcodes_rev, send_stats=stream)
        for samp_name, sample in self.samples.items():
            if samp_name == self.expected["sample"]:
                barcode_expected = sample["BarcodeRevAttrs"]["Seq"]
        self.assertEqual(assigned, barcode_expected)
        self.assertEqual(stream.getvalue(), self.expected["stream_rev"])

    def test_demux(self):
        """Try demuxing a read pair using R1 and I1."""
        path_pat = str(DATA / "demux/demux_%s_%s.fastq.gz") % (self.record, "%s")
        fps = {key: path_pat % key for key in ("R1", "R2", "I1")}
        stream = StringIO()
        with TemporaryDirectory() as outdir:
            demux(self.samples, fps, outdir, send_stats=stream)
            files_observed = sorted([str(x.name) for x in Path(outdir).iterdir()])
            self.assertEqual(files_observed, self.expected["files"])
            # Grab a simple structure of the read data: for each sample, a list
            # of lists of ID/seq for R1, R2, I1.
            sample_data = self.load_trio(outdir)
        self.assertEqual(sample_data, self.expected["demuxed_data"])

    def load_input(self):
        """Load original R1/R2/I1 for a single case name into dict of SeqRecords."""
        path_pat = str(DATA / "demux/demux_%s_%s.fastq.gz") % (self.record, "%s")
        fps = {key: path_pat % key for key in ("R1", "R2", "I1")}
        gzp = lambda fp: gzip.open(fp, "rt")
        parse = lambda f: SeqIO.parse(f, "fastq")
        data = {"R1": [], "R2": [], "I1": []}
        with gzp(fps["R1"]) as r1_in, gzp(fps["R2"]) as r2_in, gzp(fps["I1"]) as i1_in:
            for trio in zip(parse(r1_in), parse(r2_in), parse(i1_in)):
                data["R1"].append(trio[0])
                data["R2"].append(trio[1])
                data["I1"].append(trio[2])
        return data

    def load_trio(self, outdir):
        """Load demuxed sample R1/R2/I1 from a directory"""
        gzp = lambda f: gzip.open(Path(outdir) / f, "rt")
        parse = lambda f: SeqIO.parse(f, "fastq")
        sample_data = {}
        for sample in self.samples.values():
            with gzp("%s.R1.fastq.gz" % sample["Sample"]) as r1_in, \
                    gzp("%s.R2.fastq.gz" % sample["Sample"]) as r2_in, \
                    gzp("%s.I1.fastq.gz" % sample["Sample"]) as i1_in:
                data = []
                for trio in zip(parse(r1_in), parse(r2_in), parse(i1_in)):
                    data.append([
                        trio[0].id, str(trio[0].seq),
                        trio[1].id, str(trio[1].seq),
                        trio[2].id, str(trio[2].seq)])
                sample_data[sample["Sample"]] = data
        return sample_data


class TestDemuxRevcmp(TestDemux):
    """Test demultiplexing with the ReverseComplement option.

    Here we have a single read pair which is a perfect match to sample1, once
    all input reads are reverse-complemented.
    """

    def setUp(self):
        self.record = "revcmp"
        self.samples = {
            "sample1": {
                "Sample": "sample1",
                "BarcodeFwdAttrs": {"Seq": "NNNNNAACTCTAA"},
                "BarcodeRevAttrs": {"Seq": "CGGGAAGA"}},
            "sample2": {
                "Sample": "sample2",
                "BarcodeFwdAttrs": {"Seq": "NNNNNNAAGGCCCT"},
                "BarcodeRevAttrs": {"Seq": "TGCAACAA"}},
            "sample3": {
                "Sample": "sample3",
                "BarcodeFwdAttrs": {"Seq": "NNNNNNNAATATGTC"},
                "BarcodeRevAttrs": {"Seq": "CATCTGTG"}}}
        prefix = ["F", "M00281:569:000000000-CN8J9:1:1101:16840:1816", "TCGGGAACTCTAA"]
        stream_fwd = "\n".join([
            "\t".join(prefix + ["NNNNNAACTCTAA", "8.0", "*"]),
            "\t".join(prefix + ["NNNNNNAAGGCCCT", "4.0", ""]),
            "\t".join(prefix + ["NNNNNNNAATATGTC", "3.0", ""])]) + "\n"
        prefix = ["R", "M00281:569:000000000-CN8J9:1:1101:16840:1816", "CGGGAAGA"]
        stream_rev = "\n".join([
            "\t".join(prefix + ["CGGGAAGA", "8.0", "*"]),
            "\t".join(prefix + ["TGCAACAA", "4.0", ""]),
            "\t".join(prefix + ["CATCTGTG", "2.0", ""])]) + "\n"
        names = sorted(self.samples.keys()) + ["unassigned"]
        self.expected = {
            "sample": "sample1",
            "demuxed_data":  {
                "sample1": [[
                    "M00281:569:000000000-CN8J9:1:1101:16840:1816",
                    "AAGCAGTGGTATCAACGCAGAGTACATGTGGGTCTCTGCAAATGAAT" # barcode removed
                    "CGCCTGAGAGTCGAGGGCACGGCCGTGTCCTCCTGTGTAAGTCTCACAGCTTGTGACGAC"
                    "AATGCCAGGGACTTAGCGGGTGCCAGTTCTTGCCTCGGTTCCTTGCATTGTCACCATTTT"
                    "ACAACTCATTTCTTGTCTTTTTCCTTTGTTTTCTTTTCACCTTCTCCTCTTCCTCCCCCC"
                    "CGTGCCCTTCTTTCTTCCCCCTCCCTCTTCCCTCTCTCTTCTTCCGTCTTCTCCTTTCTT"
                    "TTAAAAATC",
                    "M00281:569:000000000-CN8J9:1:1101:16840:1816",
                    "CGCTGAGGAGACGGTGACCCGAACTCCCCGGCCCCATACATCCAATTATTTGTCTCTTTT"
                    "TTTCTATTCCATTTCCCCCTCCAAGAACTTGCCCCCGCTCAGTCCCTGGCATTGTCGTCT"
                    "CACTCTGTTCTACTTACACCTGTTTTCTCTCCCTTTCCCTCTCCTCTCAGGCCCTTCATT"
                    "TTCTTCTTCCCCCCTCTTCTCTTCGTTTTTACCCCTTCTTTTAGATTTCCCTACTTTCGG"
                    "AAGTGCGTCCTTTCGGGTCATAGTTTATCTCTCCGTTTTCCCCGTTTCCTTTTTACTTCC"
                    "TCTTTCTTC",
                    "M00281:569:000000000-CN8J9:1:1101:16840:1816",
                    "TCTTCCCG"]],
                "sample2": [],
                "sample3": []},
            "files": sorted(
                ["%s.R1.fastq.gz" % sample for sample in names] +
                ["%s.R2.fastq.gz" % sample for sample in names] +
                ["%s.I1.fastq.gz" % sample for sample in names]),
            "stream_fwd": stream_fwd,
            "stream_rev": stream_rev
            }

    def test_assign_barcode_fwd(self):
        keys = sorted(self.samples.keys())
        barcodes_fwd = [self.samples[s]["BarcodeFwdAttrs"]["Seq"] for s in keys]
        stream = StringIO()
        record = self.load_input()
        r1rec = record["R1"][0].reverse_complement(id=True)
        assigned = assign_barcode_fwd(r1rec, barcodes_fwd, send_stats=stream)
        for samp_name, sample in self.samples.items():
            if samp_name == self.expected["sample"]:
                barcode_expected = sample["BarcodeFwdAttrs"]["Seq"]
        self.assertEqual(assigned, barcode_expected)
        self.assertEqual(stream.getvalue(), self.expected["stream_fwd"])

    def test_assign_barcode_rev(self):
        keys = sorted(self.samples.keys())
        barcodes_rev = [self.samples[s]["BarcodeRevAttrs"]["Seq"] for s in keys]
        stream = StringIO()
        record = self.load_input()
        i1rec = record["I1"][0].reverse_complement(id=True)
        assigned = assign_barcode_rev(i1rec, barcodes_rev, send_stats=stream)
        for samp_name, sample in self.samples.items():
            if samp_name == self.expected["sample"]:
                barcode_expected = sample["BarcodeRevAttrs"]["Seq"]
        self.assertEqual(assigned, barcode_expected)
        self.assertEqual(stream.getvalue(), self.expected["stream_rev"])

    def test_demux(self):
        path_pat = str(DATA / "demux/demux_%s_%s.fastq.gz") % (self.record, "%s")
        fps = {key: path_pat % key for key in ("R1", "R2", "I1")}
        stream = StringIO()
        with TemporaryDirectory() as outdir:
            demux(self.samples, fps, outdir, dorevcmp=True, send_stats=stream)
            files_observed = sorted([str(x.name) for x in Path(outdir).iterdir()])
            self.assertEqual(files_observed, self.expected["files"])
            sample_data = self.load_trio(outdir)
        self.assertEqual(sample_data, self.expected["demuxed_data"])


class TestDemuxIDMismatch(TestDemux):
    """Test demultiplexing with mismatched sequence IDs.

    This can come up for example if you inadvertantly mix bcl2fastq I1 reads
    with MiSeq-written R1/R2 reads because the read order can be different.
    Individual fwd/rev works fine but file-based should complain.
    """

    def setUp(self):
        self.record = "idmismatch"
        self.samples = {
            "sample1": {
                "Sample": "sample1",
                "BarcodeFwdAttrs": {"Seq": "NNNNNAACTCTAA"},
                "BarcodeRevAttrs": {"Seq": "CGGGAAGA"}},
            "sample2": {
                "Sample": "sample2",
                "BarcodeFwdAttrs": {"Seq": "NNNNNNAAGGCCCT"},
                "BarcodeRevAttrs": {"Seq": "TGCAACAA"}},
            "sample3": {
                "Sample": "sample3",
                "BarcodeFwdAttrs": {"Seq": "NNNNNNNAATATGTC"},
                "BarcodeRevAttrs": {"Seq": "CATCTGTG"}}}
        prefix = ["F", "M05588:232:000000000-CLL8M:1:1101:19922:1593", "CGGGCTAAGGCCCT"]
        stream_fwd = "\n".join([
            "\t".join(prefix + ["NNNNNNAAGGCCCT", "8.0", "*"]),
            "\t".join(prefix + ["NNNNNNNAATATGTC", "4.0", ""]),
            "\t".join(prefix + ["NNNNNAACTCTAA", "4.0", ""])]) + "\n"
        prefix = ["R", "M05588:232:000000000-CLL8M:1:1101:12845:1660", "TGCAACAA"]
        stream_rev = "\n".join([
            "\t".join(prefix + ["TGCAACAA", "8.0", "*"]),
            "\t".join(prefix + ["CGGGAAGA", "4.0", ""]),
            "\t".join(prefix + ["CATCTGTG", "3.0", ""])]) + "\n"
        self.expected = {
            "sample": "sample2",
            "stream_fwd": stream_fwd,
            "stream_rev": stream_rev
            }

    def test_demux(self):
        """Test demux for mismatched sequence IDs.

        Here we should get a DemuxError because the sequence ID isn't identical
        between R1/R2/I1.
        """
        path_pat = str(DATA / "demux/demux_%s_%s.fastq.gz") % (self.record, "%s")
        fps = {key: path_pat % key for key in ("R1", "R2", "I1")}
        stream = StringIO()
        with TemporaryDirectory() as outdir:
            with self.assertRaises(DemuxError):
                demux(self.samples, fps, outdir, dorevcmp=True, send_stats=stream)
