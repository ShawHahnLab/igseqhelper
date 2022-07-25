"""
This ambiguous alignment stuff is currently a mess.  proceed with caution.
"""

import csv
import logging
from warnings import warn
from Bio import SeqIO, AlignIO
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import Array
from Bio.Seq import Seq

warn(DeprecationWarning("Don't use this junk!"), stacklevel=2)

LOGGER = logging.getLogger(__name__)

# https://www.bioinformatics.org/sms/iupac.html
# sorted alphabetically
IUPAC_DNA = {
    "A": ("A", ),
    "C": ("C", ),
    "G": ("G", ),
    "T": ("T", ),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("C", "G"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N":  ("A", "C", "G", "T")
    }

def make_iupac_scores(s_match=1, s_mismatch=-1):
    """Create lookup table of IUPAC code pairs to match/mismatch scores.

    Return value is a dictionary containg tuple pairs of single-letter codes,
    where values are scores for matches or mismatches.  Code pairs containing
    overlapping nucleotides will be considered a match, anything else, a
    mismatch.

    For example:

      (A, A): match
      (A, N): match (N contains A)
      (A, T): mismatch
      (C, S): match (S contains C)
      (S, Y): match (both contain C)
      (W, S): mismatch (no overlap)
    """
    scores = {}
    for ambig in IUPAC_DNA:
        bases1 = IUPAC_DNA[ambig]
        for ambig2 in IUPAC_DNA:
            bases2 = IUPAC_DNA[ambig2]
            for base1 in bases1:
                if base1 in bases2:
                    match = True
                    break
            else:
                match = False
            scores[(ambig, ambig2)] = {True: s_match, False: s_mismatch}[match]
            # special case: N implies neutral
            if ambig == "N" or ambig2 == "N":
                scores[(ambig, ambig2)] = 0
    return scores

def ambiguify_alignment(fasta_fp):
    """Turn an aligned FASTA into a single ambiguous sequence."""
    aln = AlignIO.read(fasta_fp, "fasta")
    bases = []
    for idx in range(aln.get_alignment_length()):
        bases_here = [rec.seq[idx] for rec in aln]
        bases_here = [b for b in bases_here if b != "-"]
        if bases_here:
            bases_here = tuple(sorted(set(bases_here)))
            for key in IUPAC_DNA:
                if bases_here == IUPAC_DNA[key]:
                    break
            else:
                raise Exception("base combination %s not in IUPAC list" % str(bases_here))
        else:
            key = "-"
        bases.append(key)
    sequence = Seq("".join(bases))
    return sequence

def score_ambiguous_matches(data_fp, ref_fp, out_fp, min_length=None):
    """Make a CSV table of ambiguous alignment scores from data to reference.

    data_fp: path fo FASTA with sequences to score.
    ref_fp: path to FASTA with aligned reference sequences to refer to as a
            single consensus.
    out_fp: path to FASTA to write output to.  The named columns are SeqID,
            Score, and Length.
    """
    LOGGER.info("score_ambiguous_matches: ambiguifying alignment from %s", ref_fp)
    ref = ambiguify_alignment(ref_fp)
    LOGGER.info("score_ambiguous_matches: generating IUPAC score LUT")
    scoring_dict = make_iupac_scores()
    LOGGER.info("score_ambiguous_matches: scoring seqs and writing output to %s", data_fp)

    aligner = PairwiseAligner()
    # https://biopython.org/docs/latest/api/Bio.Align.substitution_matrices.html
    array = Array(dims=2)
    array.update(scoring_dict)
    aligner.substitution_matrix = array
    aligner.gap_score = -5
    aligner.query_left_gap_score = -100
    aligner.query_right_gap_score = -100

    with open(out_fp, "w") as f_out:
        writer = csv.writer(f_out)
        #writer.writerow(["SeqID", "Score", "AlignmentStart", "AlignmentEnd", "SeqLength"])
        writer.writerow(["SeqID", "Score", "SeqLength"])
        for record in SeqIO.parse(data_fp, "fasta"):
            if min_length and len(record) < min_length:
                LOGGER.debug(
                    "score_ambiguous_matches: dropping %s, len %d < %d",
                    record.id,
                    len(record),
                    min_length)
                continue

            # PairwiseAligner
            alignments = aligner.align(ref, record.seq)
            print(alignments[0])
            score = alignments.score

            LOGGER.debug(
                "score_ambiguous_matches: score %s: %d",
                record.id,
                score)
            writer.writerow([record.id, score, len(record)])
