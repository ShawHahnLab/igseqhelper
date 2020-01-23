"""
Demultiplexing helpers.

This approach starts from raw reads and removes the forward barcode (including
varying-length prefix nucleotides) from the output.

Note that cutadapt has some support for demultiplexing, but from a quick read
through the user guide I don't think it's flexible enough to handle this
weirdness.
"""

import sys
import gzip
import logging
import os
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from igseq.util import revcmp

LOGGER = logging.getLogger(__name__)
PID = os.getpid()

DEMUX_ALIGNER = PairwiseAligner()
# target = first argument, query = second argument
DEMUX_ALIGNER.target_internal_gap_score = -1000
DEMUX_ALIGNER.query_internal_gap_score = -1000


class DemuxError(Exception):
    """A special exception for demultiplexing problems."""


def demux(samples, fps, outdir=".", output_files=None, revcmp=False, send_stats=sys.stdout):
    """Demultiplex one run based on dictionaries of sample and barcode data.

    samples: list of sample dicts for this run
    fps: dict of "R1", "R2", and "I1" keys pointing to file paths to fastq.gz.
    outdir: output directory to write demultiplexed fastq.gz files to
    output_files: expected output files; will touch empty ones if any aren't expected
    revcmp: reverse-complement R1 and R2 during demultiplexing?
    send_stats: file object to write a table of per-read demux results to
    """
    LOGGER.debug("%d demux: revcmp: %s", PID, str(revcmp))
    if not output_files:
        output_files = []
    else:
        output_files = [Path(fp) for fp in output_files]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for key in fps:
        LOGGER.debug("%d demux: %s input file: %s", PID, key, fps[key])
    # NOTE
    # with too many samples at once, this will cause an OS error due to too
    # many open files. In that case we'd have to open/close as needed.  It's
    # easy here to just open a bunch and store handles in a dictionary, though.
    fp_outs = {
        "R1": {s["Sample"]: Path(outdir) / (s["Sample"] + ".R1.fastq.gz") for s in samples},
        "R2": {s["Sample"]: Path(outdir) / (s["Sample"] + ".R2.fastq.gz") for s in samples},
        "I1": {s["Sample"]: Path(outdir) / (s["Sample"] + ".I1.fastq.gz") for s in samples}
        }
    fp_outs["R1"]["None"] = Path(outdir) / ("unassigned.R1.fastq.gz")
    fp_outs["R2"]["None"] = Path(outdir) / ("unassigned.R2.fastq.gz")
    fp_outs["I1"]["None"] = Path(outdir) / ("unassigned.I1.fastq.gz")
    for key in fp_outs:
        for sample in fp_outs[key]:
            LOGGER.debug("%d demux: %s %s output file: %s", PID, key, sample, fp_outs[key][sample])
    try:
        f_outs = {
            "R1": {key: gzip.open(fp_outs["R1"][key], "wt") for key in fp_outs["R1"].keys()},
            "R2": {key: gzip.open(fp_outs["R2"][key], "wt") for key in fp_outs["R2"].keys()},
            "I1": {key: gzip.open(fp_outs["I1"][key], "wt") for key in fp_outs["I1"].keys()}
            }
        _trio_demux(fps, f_outs, samples, revcmp, send_stats)
    finally:
        for f_out in f_outs["R1"].values():
            f_out.close()
        for f_out in f_outs["R2"].values():
            f_out.close()
        for f_out in f_outs["I1"].values():
            f_out.close()
    fp_outs_all = \
        list(fp_outs["R1"].values()) + \
        list(fp_outs["R2"].values()) + \
        list(fp_outs["I1"].values())
    for extra in set(output_files) - set(fp_outs_all):
        with gzip.open(extra, "wt") as _:
            pass

def gzp(fp_in):
    """Convenience wrapper for gzip opener."""
    return gzip.open(fp_in, "rt")

def fqparse(f_in):
    """Convenience wrapper for FASTQ parser."""
    return SeqIO.parse(f_in, "fastq")

def _trio_demux(fps, f_outs, samples, dorevcmp, send_stats):
    barcodes = {
        "F": [s["BarcodeFwd"] for s in samples],
        "R": [s["BarcodeRev"] for s in samples]
        }
    bc_map = {(s["BarcodeFwd"], s["BarcodeRev"]): s["Sample"] for s in samples}
    for key in bc_map:
        LOGGER.debug("%d barcode map: %s -> %s", PID, str(key), bc_map[key])
    counter = 0
    hits = defaultdict(int)
    with gzp(fps["R1"]) as r1_in, gzp(fps["R2"]) as r2_in, gzp(fps["I1"]) as i1_in:
        for trio in zip(fqparse(r1_in), fqparse(r2_in), fqparse(i1_in)):
            if not trio[0].id == trio[1].id == trio[2].id:
                raise DemuxError("Sequence ID mismatch between R1/R2/I1")
            trio = list(trio)
            if dorevcmp:
                trio[0] = revcmp(trio[0])
                trio[1] = revcmp(trio[1])
                trio[2] = revcmp(trio[2])
            assigned = (
                assign_barcode_fwd(trio[0], barcodes["F"], send_stats=send_stats),
                assign_barcode_rev(trio[2], barcodes["R"], send_stats=send_stats))
            samp = bc_map.get(assigned)
            hits[assigned] += 1
            # trim barcode from forward read
            if assigned[0]:
                trio[0] = trio[0][len(assigned[0]):]
            SeqIO.write(trio[0], f_outs["R1"][str(samp)], "fastq")
            SeqIO.write(trio[1], f_outs["R2"][str(samp)], "fastq")
            SeqIO.write(trio[2], f_outs["I1"][str(samp)], "fastq")
            # Log message every 10,000th read trio
            if counter % 10000 == 0:
                keys = sorted(hits, key=hits.__getitem__)[::-1]
                if len(keys) > 1:
                    LOGGER.debug(
                        "%d read trio %d: top hits so far: %s (%d reads), %s (%d reads)",
                        PID,
                        counter,
                        str(keys[0]), hits[keys[0]],
                        str(keys[1]), hits[keys[1]])
            counter += 1


def assign_barcode_fwd(record, barcodes,
                       max_mismatch=1, min_next_mismatch=1, send_stats=sys.stdout):
    """Assign a forward barcode to a single R1 record.

    The forward barcodes come with a varying (known per barcode but not before
    we've assigned barcodes) prefix of at least four random nucleotides.

    record: SeqRecord object for R1 with one read
    barcodes: list of known forward barcodes, N included
    max_mismatch: how many mismatches to a barcode should be allowed?
    min_next_mismatch: how far away does a match need to be from the next-closest barcode?
    send_stats: file object to write a table of per-read demux results to
    """
    # The number of unknown nucleotides varies depending on the barcode, and we
    # don't know which barcode we're dealing with.  We'll just take a prefix
    # from the read that  includes enough to handle any barcode in use here.
    num_ns = [len(bcode) - len(bcode.lstrip("N")) for bcode in barcodes]
    prefix = record.seq[0:(max(num_ns) + 8)]
    # Align the N-stripped barcodes to the read prefix, and get the scores.
    align = lambda bcode: DEMUX_ALIGNER.align(prefix, bcode.lstrip("N"))
    aligns = {bcode: align(bcode) for bcode in barcodes}
    scores = {bcode: aligns[bcode].score for bcode in aligns.keys()}
    match = _match_barcode(scores, max_mismatch, min_next_mismatch)
    if match:
        prefix = prefix[0:len(match)]
    _send_stats(send_stats, "F", record.id, prefix, scores, match)
    return match

def assign_barcode_rev(record, barcodes,
                       max_mismatch=1, min_next_mismatch=1, send_stats=sys.stdout):
    """Assign a reverse barcode to a single I1 record.

    The reverse barcode is stored in the I1 read as the reverse-complement of
    what's listed in the protocol tables.

    record: SeqRecord object for I1 with one read
    barcodes: list of known reverse barcodes
    max_mismatch: how many mismatches to a barcode should be allowed?
    min_next_mismatch: how far away does a match need to be from the next-closest barcode?
    send_stats: file object to write a table of per-read demux results to
    """
    # Reverse is considerably simpler than forward, since we just need to take
    # the reverse-complement of the entire I1 read.
    prefix = record.seq.reverse_complement()
    align = lambda bcode: DEMUX_ALIGNER.align(prefix, bcode)
    aligns = {bcode: align(bcode) for bcode in barcodes}
    scores = {bcode: aligns[bcode].score for bcode in aligns.keys()}
    match = _match_barcode(scores, max_mismatch, min_next_mismatch)
    _send_stats(send_stats, "R", record.id, prefix, scores, match)
    return match

def _match_barcode(scores, max_mismatch, min_next_mismatch):
    """Select the best-matching barcode from scores and thresholds."""
    # The first and next-best scores
    match_score = max(scores.values())
    next_scores = [val for val in scores.values() if val < match_score]
    # but there may not be a next best, so be careful
    if next_scores:
        match_score_second = max(next_scores)
    else:
        match_score_second = 0
    # Match the best score with the first barcode that produced it
    is_match = lambda b: scores[b] == match_score
    match = [barcode for barcode in scores.keys() if is_match(barcode)][0]
    # Count the number of mismatches for the best and next best
    # The N removal is only needed for fwd but is fine to leave the same for
    # rev
    mismatch = len(match.lstrip("N")) - match_score
    mismatch_next = len(match.lstrip("N")) - match_score_second
    # Invalidate the barcode match if it isn't definitive enough
    if not (mismatch <= max_mismatch and mismatch_next >= min_next_mismatch):
        match = None
    return match

def _send_stats(send_stats, read_dir, recid, prefix, scores, match):
    """Log barcode assignment details to TSV on a file handle, if defined.

    Columns are:
     1. "F" or "R" for forward or reverse read
     2. read ID
     3. the matched region of the input sequence
     4. the barcode compared with
     5. the score for that barcode
     6. "*" If that barcode was taken as the match
    """
    # pylint: disable=too-many-arguments
    if send_stats:
        for barcode in sorted(scores, key=scores.__getitem__)[::-1]:
            send_stats.write("\t".join([
                read_dir, # F or R
                recid, # sequence ID
                str(prefix), # sequence segment matched with
                barcode, # barcode under consideration
                str(scores[barcode]), # result for this barcode
                {True: "*", False: ""}[match == barcode]])) # barcode chosen?
            send_stats.write("\n")
