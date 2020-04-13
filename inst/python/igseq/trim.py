"""
Helpers for trimming raw reads.
"""

import logging
from igseq.util import revcmp
from igseq.data import MetadataError

LOGGER = logging.getLogger(__name__)

def adapter_fwd(sample, sequences):
    """Get the adapter sequence to trim off the end of R1.

    sample: dictionary of sample attributes
    sequences: dictionary of sequence attributes
    """
    # The PCR primer specific to the antibody type occurs just *after* the
    # beginning of R2 (in the 3' direction, that is), so we'll trim that off
    # the end of R1.
    # Each sample's chain type maps to a single type-specific PCR primer in the
    # sequences metadata.
    type_map = {
        "gamma": "IgG",
        "delta": "IgD",
        "kappa": "IgK",
        "lambda": "IgL",
        "mu": "RhIgM"}
    antibody_type = sample["Type"]
    try:
        pcr_seq_name = type_map[antibody_type]
    except KeyError:
        raise MetadataError("Unknown antibody type %s" % antibody_type)
    try:
        adapter = sequences[pcr_seq_name]["Seq"]
    except KeyError:
        raise MetadataError("Missing entry in sequence metadata: %s" % pcr_seq_name)
    return revcmp(adapter)

def adapter_rev(sample, sequences):
    """Get the adapter sequence to trim off the end of R2.

    sample: dictionary of sample attributes
    sequences: dictionary of sequence attributes
    """
    # Reverse is trickier than forward.  The 2nd round PCR Forward Primer 1 /
    # P5 Sequencing Primer Site occurs just *before* the beginning of R1, so we
    # could cut on that, but we don't want to leave a dangling barcode segment
    # to get paired back in when combining R1 and R2 later (it's preferable if
    # our paired sequences start exactly after the forward barcode).  So we'll
    # include the barcode for each specific sample as part of the adapter here.
    # Note that the "N" characters in the barcode sequences should be handled
    # just fine by cutadapt.
    try:
        p5seq = sequences["P5_Seq"]["Seq"]
    except KeyError:
        raise MetadataError("Missing entry in sequence metadata: P5_Seq")
    try:
        bcattrs = sample["BarcodeFwdAttrs"]
    except KeyError:
        LOGGER.error(
            "adapter_rev: No forward barcode for %s, trimming by P5_Seq only",
            sample["Sample"])
        return revcmp(p5seq)
    barcode = bcattrs["Seq"]
    # Reverse complement the barcode and prepend to the constant region.
    return revcmp(barcode) + revcmp(p5seq)
