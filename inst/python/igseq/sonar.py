"""
Helpers for running SONAR.
"""

import re
import logging
from Bio import SeqIO

LOGGER = logging.getLogger(__name__)

def gather_germline(input_imgt_fp, input_extra_v_fp, output_fp, segment):
    """
    Prepare germline references for SONAR.

    Note that it expects a particular sequence naming scheme for genes and
    module 2 will fail if we don't follow it.
    """
    with open(output_fp, "w") as f_out:
        seqids_seen = set()
        for record in SeqIO.parse(input_imgt_fp, "fasta"):
            new_id = munge_seqid_for_sonar(record.id, seqids_seen)
            record.id = new_id
            record.description = ""
            SeqIO.write(record, f_out, "fasta")
        if segment == "V":
            for record in SeqIO.parse(input_extra_v_fp, "fasta"):
                new_id = munge_seqid_for_sonar(record.id, seqids_seen)
                record.id = new_id
                record.description = ""
                SeqIO.write(record, f_out, "fasta")

def munge_seqid_for_sonar(seqid, seqids_seen):
    """Modify given sequence ID to work with SONAR.

    SONAR expects a particular format for the sequence ID entries for V(D)J
    genes and alleles.  This will try to parse out  the right portion and
    return it.  If it can't be found, an error is logged and the original ID is
    returned instead.  This works with IMGT IDs and at least one other
    reference database we're using.

    Also, the given set of already-observed sequence IDs will be used to record
    each new ID and check against the previous set.  An error is logged if a
    sequence ID has already been observed.  Note that the set is modified in place.
    """
    # See lineage/2.1-calculate_id-div.py
    # (IG[HKL]V[^*]+)[^,\s]+
    match = re.search(r"\|?(IG[HKL][VDJ][-*0-9A-Za-z_]+)\|?", seqid)
    if match:
        new_id = match.group(1)
    else:
        LOGGER.error("Sequence ID format not recognized, SONAR might get upset: %s", seqid)
        new_id = seqid
    if new_id in seqids_seen:
        LOGGER.error("Sequence ID already seen in this set: %s", new_id)
    seqids_seen.add(new_id)
    return new_id
