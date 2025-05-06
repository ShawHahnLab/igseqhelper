#!/usr/bin/env python

"""
Make a quality-aware consensus sequence from FASTQ reads
"""

from pathlib import Path
from argparse import ArgumentParser
from collections import defaultdict
from Bio import Align
from igseq.record import RecordReader, RecordWriter, RecordHandler

ALIGNER = Align.PairwiseAligner()
ALIGNER.gap_score = -10

def _load_fastq(fqgz_in):
    ref = None
    # index in ref -> base call -> list of Q scores
    scores = defaultdict(lambda: defaultdict(list))
    with RecordReader(fqgz_in) as reader:
        for row in reader:
            if not ref:
                ref = row["sequence"]
            aln = ALIGNER.align(ref, row["sequence"])[0]
            quals = RecordHandler.decode_phred(row["sequence_quality"])
            indices = aln.indices[:]
            for idx in range(aln.length):
                idx_target = indices[0][idx] # index in original target (ref) seq
                idx_query = indices[1][idx] # index in query seq
                # -1 for gaps in target, which we don't care about, gaps in
                # read, which won't affect things anyway
                if idx_target != -1 and idx_query != -1:
                    basecall = aln[1][idx]
                    scores[idx_target][basecall].append(quals[idx_query])
    return scores

def _make_consensus(scores):
    consensus = ""
    quals = []
    check = 0
    for idx, base_scores in scores.items():
        assert check == idx # TODO rearrange this stupid dictionary approach
        # numbers of observed base calls per base, for use in sorting by
        # decreasing abundance
        tally = {base: len(scores_here) for base, scores_here in base_scores.items()}
        # reformatting into pairs of base and quality score
        pairs = []
        for base, scores_here in base_scores.items():
            for score in scores_here:
                pairs.append((base, score))
        # sort by quality first, then by abundance of base call
        pairs.sort(key=lambda pair, t=tally: (pair[1], t[pair[0]]), reverse=True)
        # top one is winner
        base, score = pairs[0]
        consensus += base
        quals.append(score)
        check += 1
    return consensus, quals

def fastq_consensus(fqgz_in, path_out):
    """Make a quality-aware consensus sequence from FASTQ reads"""
    # first, gather up base calls and associated Q scores at each position in
    # the reference (which is presumed to be the first read in the file)
    # Sort all entries by quality and then by decreasing abundance (across all
    # qualities) per base call, and just take the top one.  Then we can also
    # just note the Q score for each selected read+position, too.
    scores = _load_fastq(fqgz_in)
    consensus, quals = _make_consensus(scores)
    path_out = Path(path_out)
    path_out.parent.mkdir(parents=True, exist_ok=True)
    with RecordWriter(path_out) as writer:
        writer.write({
            "sequence_id": "consensus",
            "sequence": consensus,
            "sequence_quality": RecordHandler.encode_phred(quals)})

def main():
    """CLI for fastq_consensus"""
    parser = ArgumentParser()
    addarg = parser.add_argument
    addarg("input", help="FASTQ input file (.fastq or .fastq.gz)")
    addarg("output", help="Output file for consensus (FASTA or FASTQ or tabular)")
    args = parser.parse_args()
    fastq_consensus(args.input, args.output)

if __name__ == "__main__":
    main()
