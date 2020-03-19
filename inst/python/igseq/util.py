"""
Misc utility functions.
"""
import re
from Bio import Phylo
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
