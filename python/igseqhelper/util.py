"""
Misc utility functions.
"""
import logging
from Bio.Seq import Seq

LOGGER = logging.getLogger(__name__)

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
