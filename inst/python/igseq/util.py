"""
Misc utility functions.
"""
import re
from Bio import Phylo
from Bio.Seq import Seq

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
