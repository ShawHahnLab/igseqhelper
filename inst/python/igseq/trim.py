"""
Helpers for trimming raw reads.
"""

from igseq import SAMPLES
from igseq.util import revcmp

def adapter_fwd(sample_name, samples=None):
    """Get the adapter sequence to trim off the end of R1."""
    # The first round PCR primer occurs just *after* the beginning of R2, so
    # we'll trim that off the end of R1.
    if not samples:
        samples = SAMPLES
    return "TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC"

def adapter_rev(sample_name, samples=None):
    """Get the adapter sequence to trim off the end of R2."""
    if not samples:
        samples = SAMPLES
    # Reverse is trickier than forward.  The 2nd round PCR Forward Primer 1 /
    # P5 Sequencing Primer Site occurs just *before* the beginning of R1, so we
    # could cut on that, but we don't want to leave a dangling barcode segment
    # to get paired back in when combining R1 and R2 later (it's preferable if
    # our paired sequences start exactly after the forward barcode).  So we'll
    # include the barcode for each specific sample as part of the adapter here.
    # Note that the "N" characters in the barcode sequences should be handled
    # just fine by cutadapt.
    p5seq = "AGATCGGAAGAGCGTCGTGTAGGGAAAGA" # (reverse complemented)
    # Extract the (should be "the" but my wonky metadata handling means we
    # could run into multiple, so we'll check against that) forward barcode to
    # use as part of the adapter.
    bcs = []
    for runid in samples:
        for samp in samples[runid]:
            if samp["Sample"] == sample_name:
                bcs.append(samp["BarcodeFwd"])
    bcs = set(bcs)
    if len(bcs) > 1:
        raise ValueError(
            ("sample %s barcoded differently in different "
             "runs (%s), cleanup your metadata!") % (sample_name, str(bcs)))
    if len(bcs) == 0:
        raise ValueError("No barcode found for sample %s" % sample_name)
    # Reverse complement the barcode and prepend to the constant region.
    barcode = revcmp(bcs.pop())
    return barcode + p5seq
