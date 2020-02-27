"""
Misc utility functions.
"""
import re
from Bio import Phylo
from Bio.Seq import Seq
from Bio import SeqIO

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
