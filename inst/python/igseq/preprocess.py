"""
Various preprocessing tools (filtering, demultiplexing, etc.)
"""
from pathlib import Path
import gzip
import re
import sys
import logging
from Bio import SeqIO
from Bio.Align import PairwiseAligner

LOGGER = logging.getLogger(__name__)

# The 5PIIA constant segment used in all samples' forward primers
PRIMER_CONSTANT = "AAGCAGTGGTATCAACGCAGAGT"
PAT_UNDET = "Undetermined_S0_L001_%s_001.fastq.gz"
MSG_INTERVAL = 1000000

DEMUX_ALIGNER = PairwiseAligner()
# target = first argument, query = second argument
DEMUX_ALIGNER.target_internal_gap_score = -1000
DEMUX_ALIGNER.query_internal_gap_score = -1000

def primerfilt(runid, revcmp, dir_path_in, dir_path_out):
    """Filter reads to only those including the expected primer region in R1.

    This is a quick and simple way to separate out the material that was part
    of the protocol from everything else (PhiX or other samples).
    """
    runid = Path(runid).name
    counter = 0
    dir_in = Path(dir_path_in) / runid
    dir_out = Path(dir_path_out) / runid
    dir_out.mkdir(parents=True, exist_ok=True)
    LOGGER.debug("primerfilt: %s dir_path_in: %s", runid, dir_path_in)
    LOGGER.debug("primerfilt: %s dir_path_out: %s", runid, dir_path_out)
    LOGGER.debug("primerfilt: %s revcmp: %s", runid, revcmp)
    with \
            gzip.open(str(dir_in / PAT_UNDET) % "R1", "rt") as r1_in, \
            gzip.open(str(dir_in / PAT_UNDET) % "R2", "rt") as r2_in, \
            gzip.open(str(dir_in / PAT_UNDET) % "I1", "rt") as i1_in:
        with \
                gzip.open(str(dir_out / PAT_UNDET) % "R1", "wt") as r1_out, \
                gzip.open(str(dir_out / PAT_UNDET) % "R2", "wt") as r2_out, \
                gzip.open(str(dir_out / PAT_UNDET) % "I1", "wt") as i1_out:
            for trio in zip( \
                    SeqIO.parse(r1_in, "fastq"), \
                    SeqIO.parse(r2_in, "fastq"), \
                    SeqIO.parse(i1_in, "fastq")):
                trio = list(trio)
                if counter % MSG_INTERVAL == 0:
                    LOGGER.info("primerfilt: %s Reads: %d\n", runid, counter)
                counter += 1
                if revcmp:
                    flip = lambda r: r.reverse_complement(
                        id=True, name=True, description=True, features=True,
                        annotations=True, letter_annotations=True, dbxrefs=True)
                    trio[0] = flip(trio[0])
                    trio[1] = flip(trio[1])
                if PRIMER_CONSTANT in trio[0].seq:
                    SeqIO.write(trio[0], r1_out, "fastq")
                    SeqIO.write(trio[1], r2_out, "fastq")
                    # I1/I2 apparently do NOT have the same sequence IDs as
                    # R1/R2 which ruins my original plan to filter downstream
                    # on ID lists.  We'll force it here.
                    trio[2].id = trio[0].id
                    SeqIO.write(trio[2], i1_out, "fastq")

def demux_record_fwd(record, barcodes_fwd,
                     max_mismatch=1, min_next_mismatch=1, send_stats=sys.stdout):
    """Demultiplex a single sequence record by pairing to a forward barcode.

    This is the older, post-assembly demultiplexing.  See demux.py for the
    demultiplexer on raw reads.

    barcodes_fwd: dict mapping barcode sequences to sample names
    max_mismatch: how many mismatches to a barcode should be allowed?
    min_next_mismatch: how far away does a match need to be from the next-closest barcode?
    send_stats: file object to write a table of per-read demux results to
    """
    # forward barcode is 8 bases with a varying prefix of unknown bases (4 on up).
    # trim off the four we *know* aren't part of the barcode, and take enough
    # to account for the longest expected barcode.
    prefix = record.seq[0:max([len(barcodes_fwd[bc]) for bc in barcodes_fwd])]
    aligns_fwd = {barcode: DEMUX_ALIGNER.align(prefix, barcode) for barcode in barcodes_fwd.keys()}
    scores_fwd = {barcode: aligns_fwd[barcode].score for barcode in aligns_fwd.keys()}
    match_fwd_score = max(scores_fwd.values())
    next_scores = [val for val in scores_fwd.values() if val < match_fwd_score]
    # there may not be a next best, so be careful
    if next_scores:
        match_fwd_score_second = max(next_scores)
    else:
        match_fwd_score_second = 0
    is_match = lambda b: scores_fwd[b] == match_fwd_score
    match_fwd = [barcode for barcode in scores_fwd.keys() if is_match(barcode)][0]
    samp_fwd = None
    if len(match_fwd) - match_fwd_score <= max_mismatch and \
        len(match_fwd) - match_fwd_score_second >= min_next_mismatch:
        samp_fwd = barcodes_fwd[match_fwd]
    if send_stats:
        send_stats.write("\t".join(
            [record.id, str(prefix)] +
            [str(scores_fwd[barcode]) for barcode in scores_fwd.keys()] +
            [match_fwd, str(samp_fwd)]))
        send_stats.write("\n")
    return samp_fwd

def demux_fwd(samples, fp_paired, fp_i1, outdir=".", output_files=None, send_stats=sys.stdout):
    """Demultiplex one run based on dictionaries of sample and barcode data.

    This is the older, post-assembly demultiplexing.  See demux.py for the
    demultiplexer on raw reads.

    samples: list of sample dicts for this run
    fp_paired: path to input fastq.gz for paired reads
    fp_i1: path to input fastq.gz for I1 reads
    outdir: output directory to write demultiplexed fastq.gz files to
    output_files: expected output files; will touch empty ones if any aren't expected
    send_stats: file object to write a table of per-read demux results to
    """
    if not output_files:
        output_files = []
    else: output_files = [Path(fp) for fp in output_files]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    barcodes_fwd = {s["BarcodeFwd"]: s["Sample"] for s in samples}
    # NOTE
    # with too many samples at once, this will cause an OS error due to too
    # many open files. In that case we'd have to open/close as needed.  It's
    # easy here to just open a bunch and store handles in a dictionary, though.
    fp_outs = {s["Sample"]: Path(outdir) / (s["Sample"] + ".fastq.gz") for s in samples}
    fp_outs["None"] = Path(outdir) / ("unassigned.fastq.gz")
    try:
        f_outs = {key: gzip.open(fp_outs[key], "wt") for key in fp_outs.keys()}
        with gzip.open(fp_paired, "rt") as pair_in, gzip.open(fp_i1, "rt") as i1_in:
            for pair in zip(SeqIO.parse(pair_in, "fastq"), SeqIO.parse(i1_in, "fastq")):
                samp = demux_record_fwd(pair[0], barcodes_fwd, send_stats=send_stats)
                SeqIO.write(pair[0], f_outs[str(samp)], "fastq")
    finally:
        for f_out in f_outs.values():
            f_out.close()
    # touch any remaining filenames
    for extra in set(output_files) - set(fp_outs.values()):
        with gzip.open(extra, "wt") as _:
            pass

def primertrim(fqgz_in_fp, fqgz_out_fp, primer_fwd):
    """Trim off matches to a given sequence at the start of each reads.

    This is a simple helper to trim off everything up to and including the
    given forward primer, with only exact matching.
    """
    with gzip.open(fqgz_in_fp, "rt") as f_in, gzip.open(fqgz_out_fp, "wt") as f_out:
        for record in SeqIO.parse(f_in, "fastq"):
            match = re.match(".*(" + primer_fwd + ").*", str(record.seq))
            if match:
                idx = match.start(1) + len(primer_fwd)
                record = record[idx:]
                SeqIO.write(record, f_out, "fastq")
