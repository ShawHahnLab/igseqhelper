"""
Misc utility functions.
"""
import re
import builtins
import logging
from math import ceil, log10
from pathlib import Path
from Bio import Phylo
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

LOGGER = logging.getLogger(__name__)

# https://www.bioinformatics.org/sms/iupac.html
# With added stubs to more easily handle regular bases and gaps
IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "T": {"T"},
    "G": {"G"},
    "-": {"-"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "T", "G"}}

def revcmp(record):
    """Reverse complement a SeqRcord, keeping ALL metadata.

    If a string is given, reverse-complement that.
    """
    try:
        return record.reverse_complement(
            id=True, name=True, description=True, features=True,
            annotations=True, letter_annotations=True, dbxrefs=True)
    except AttributeError:
        return Seq(record).reverse_complement()

def parse_seq_desc(record):
    """Parse key=val pairs from a SeqRecord description to a dict.

    If a sring is given, take that as the description.  Any entries without "="
    are skipped.
    """
    try:
        txt = record.description
    except AttributeError:
        txt = record
    attrs = [pair.split("=", 1) for pair in txt.split(" ")]
    attrs = {pair[0]: pair[1] for pair in attrs if len(pair) == 2}
    return attrs

def strdist(str1, str2, compare=lambda chr1, chr2: chr1 != chr2):
    """Naive string distance comparison.

    This just compares each pair of characters using the given function
    (not-equals by default).  The optional comare argument should be a function
    that returns 0 if a pair of characters match, 1 otherwise.
    """
    return sum([compare(x, y) for x, y in zip(str1, str2)])

def strdist_iupac(str1, str2):
    """String distance considering ambiguous bases as matches.

    For example, A <-> A gives 0, A <-> R gives 0, A <-> Y gives 1.
    Unrecognized characters count as a mismatch.
    """
    def compare(chr1, chr2):
        # It's a mismatch if all of these things are true:
        # 1) the characters are not the same
        # 2) the set of characters represented by the 1st doesn't have the 2nd
        # 3) the set of characters represented by the 2nd doesn't have the 1st
        return chr1 != chr2 and \
            chr2 not in IUPAC.get(chr1, []) and \
            chr1 not in IUPAC.get(chr2, [])
    return strdist(str1, str2, compare)

def strdist_iupac_squeezed(str1, str2):
    """String distance considering ambiguous bases as matches, with trimming.

    As in strdist_iupac, but excluding gap/nongap comparison at edges.
    For example this will give 2:
    str1:       CGCAC-TG------
    str2:       ---ACGTCGCACGC
    mismatches:      X X
    """
    def compare(chr1, chr2):
        return chr1 != chr2 and \
            chr2 not in IUPAC.get(chr1, []) and \
            chr1 not in IUPAC.get(chr2, [])
    match1 = re.match("(^-*)([^-].*[^-])(-*)$", str1)
    match2 = re.match("(^-*)([^-].*[^-])(-*)$", str2)
    start = max(match1.end(1), match2.end(1))
    stop = min(match1.start(3), match2.start(3))
    str1_squeezed = str1[start:stop]
    str2_squeezed = str2[start:stop]
    return strdist(str1_squeezed, str2_squeezed, compare)

def ambiguify_alignment(fasta_fp_in, fasta_fp_out, seqid="CombinedAlignment"):
    """Condense all variation in each position to IUPAC codes.

    This will take a FASTA alignment over a number of sequences, and output a
    single sequence that encompasses all variation at each position in the
    appropriate IUAPC ambiguity code.

    NOTE: As currently written this doesn't know how to account for a gap, so
    it will become "?" for any positions with gaps.  Sequences that vary in
    length will cause a ValueError.
    """
    seqs_in = []
    for record in SeqIO.parse(fasta_fp_in, "fasta"):
        seqs_in.append(str(record.seq))
    if not seqs_in:
        with open(fasta_fp_out, "w") as _:
            pass
        return
    lengths = {len(s) for s in seqs_in}
    if len(lengths) > 1:
        raise ValueError("Input sequences vary in length (not an alignment?)")
    # Set up a dictionary mapping sets of explicit bases to their IUPAC
    # ambiguity codes.
    # (I know because of how I set up IUPAC above that it's 1:1 with keys and
    # values, but be careful...)
    iupac_rev = {"".join(sorted(list(val))): key for key, val in IUPAC.items()}
    def getbase(pos):
        chars = "".join(sorted(list({seq[pos] for seq in seqs_in})))
        return iupac_rev.get(chars, "?")
    bases = [getbase(pos) for pos in range(len(seqs_in[0]))]
    seq_out = "".join(bases)
    record = SeqRecord(Seq(seq_out), id=seqid, description="")
    SeqIO.write(record, fasta_fp_out, "fasta")

def trim_alignment_to(fasta_fp_in, fasta_fp_out, seq_id):
    """Trim all sequences in a FASTA to match a given sequence ID.

    This removes start and end regions of other sequences so they line up with
    one particular sequence in the alignment.
    """
    for record in SeqIO.parse(fasta_fp_in, "fasta"):
        if record.id == seq_id:
            match = re.match("^-*(.*[^-])-*$", str(record.seq))
            start = match.start(1)
            end = match.end(1)
            break
    else:
        raise Exception("Seq ID %s not found" % seq_id)
    with open(fasta_fp_out, "w") as f_out:
        for record in SeqIO.parse(fasta_fp_in, "fasta"):
            SeqIO.write(record[start:end], f_out, "fasta")

def get_clade_seq_ids(tree_fp, query_ids, tree_fmt="nexus"):
    """Return all sequence IDs in the same clade as query IDs.

    tree_fp: path to tree on disk
    query_id: list of sequence IDs to use to define the clade
    tree_fmt: format of tree file
    """
    tree = Phylo.read(tree_fp, tree_fmt)
    clade = tree.common_ancestor(query_ids)
    ids = [node.name for node in clade.get_terminals()]
    ids = [re.sub("^'(.*)'$", "\\1", seq_id) for seq_id in ids]
    return ids

# Vaguely based on a few answers to https://stackoverflow.com/questions/952914
# Why did I need this again?
def flatten(iterable):
    """Take a nested iterable and retun a flat list with its contents."""
    flat = []
    for obj in iterable:
        try:
            flat.extend(flatten(obj))
        except TypeError:
            flat.append(obj)
    return flat

def normalize_read_files(fps):
    """Handle either a dict of R1/R2/I1 or a directory path, creating a list of dicts."""
    try:
        if Path(fps).is_dir():
            fps = Path(fps).glob("*.fastq.gz")
        else:
            raise ValueError("give either R1/R2/I1 dict or directory path for demux")
    except TypeError:
        return [fps]
    sets = {"R1": [], "R2": [], "I1": []}
    samples = {"R1": [], "R2": [], "I1": []}
    pattern = r"([^/]+)_S([0-9]+)_L[0-9]{3}_(R1|R2|I1|I2)_001*\.fastq\.gz"
    for path in fps:
        match = re.search(pattern, str(path))
        if match.group(3) in sets:
            sets[match.group(3)].append(path)
            samples[match.group(3)].append(match.group(1))
    if not samples["R1"] == samples["R2"] == samples["I1"]:
        raise ValueError("Sample name mismatch between R1/R2/I1")
    sets = zip(sets["R1"], sets["R2"], sets["I1"])
    return [dict(zip(["R1", "R2", "I1"], s)) for s in sets]

def make_chunk_str(count):
    """10 -> ["001", "002", ..., "010"]"""
    pad = max(3, ceil(log10(count)))
    return [str(num+1).zfill(pad) for num in range(count)]

class RoundRobinWriter:
    """Split write calls cyclically to a list of files.

    For example, if write is called 1000 times and three output filenames are
    given, the first output file will have 334 writes and the second and third
    will have 333 each.  If write is only called once, the second and third
    file will exist but will be left empty. RoundRobinWriter objects implement
    just three file-like methods: open, close, and write.  For this to work
    each call to write must contain exactly one record.
    """

    def __init__(self, file_paths_out, opener=builtins.open, **kwargs):
        """Set up output filenames and opener function.

        open() isn't called yet.  That happens either explicitly or when
        __enter__ is called as in "with RoundRobinWriter(...) as ...".  By
        default the builtin file open function is used, but another open
        function can supplied such as gzip.open.  Any other arguments are
        passed to the opener when it is called.
        """
        LOGGER.info("RoundRobinWriter: file_paths_out: %s", file_paths_out)
        LOGGER.info("RoundRobinWriter: opener: %s", opener)
        self.file_paths_out = file_paths_out
        self.files_out = None
        self.opener = opener
        self.idx = 0
        self.kwargs = kwargs
        self.kwargs["mode"] = self.kwargs.get("mode", "wt")

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type:
            LOGGER.error("RoundRobinWriter: exc_type: %s", exc_type)
            LOGGER.error("RoundRobinWriter: exc_value: %s", exc_value)
            LOGGER.error("RoundRobinWriter: traceback: %s", traceback)
        self.close()

    def open(self):
        """Open output files for writing."""
        LOGGER.info("RoundRobinWriter: open")
        self.files_out = [self.opener(fp, **self.kwargs) for fp in self.file_paths_out]

    def close(self):
        """Close all open file handles."""
        LOGGER.info("RoundRobinWriter: close")
        for obj in self.files_out:
            obj.close()

    def write(self, txt):
        """File-like write method.

        Bio.SeqIO.write calls _FormatToString[format] to get each txt for each
        record, and calls fp.write once for each record.
        """
        self.files_out[self.idx].write(txt)
        self.idx = (self.idx + 1) % len(self.file_paths_out)
