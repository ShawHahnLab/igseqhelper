"""
Helpers for trimming raw reads.
"""

from igseq.util import revcmp

def adapter_fwd(sample_name, samples):
    """Get the adapter sequence to trim off the end of R1."""
    # The first round PCR primer occurs just *after* the beginning of R2, so
    # we'll trim that off the end of R1.
    # TODO this should be handled via sequences.csv
    return "TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC"

def adapter_rev(sample_name, samples):
    """Get the adapter sequence to trim off the end of R2."""
    # Reverse is trickier than forward.  The 2nd round PCR Forward Primer 1 /
    # P5 Sequencing Primer Site occurs just *before* the beginning of R1, so we
    # could cut on that, but we don't want to leave a dangling barcode segment
    # to get paired back in when combining R1 and R2 later (it's preferable if
    # our paired sequences start exactly after the forward barcode).  So we'll
    # include the barcode for each specific sample as part of the adapter here.
    # Note that the "N" characters in the barcode sequences should be handled
    # just fine by cutadapt.
    # TODO this should be handled via sequences.csv
    p5seq = "AGATCGGAAGAGCGTCGTGTAGGGAAAGA" # (reverse complemented)

    barcode = samples[sample_name]["BarcodeFwdAttrs"]["Seq"]
    # Reverse complement the barcode and prepend to the constant region.
    barcode = revcmp(barcode)
    return barcode + p5seq
