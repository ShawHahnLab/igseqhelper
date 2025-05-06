#!/usr/bin/env python

"""
Replace SONAR cluster centroid sequences in an AIRR TSV with updated versions.

This uses a directory of sequence files, presumed to be named
{centroid}.{suffix} with one record each, to supply the new centroid sequence
for each centroid row in the given SONAR rearrangements TSV file.  The sequence
alignment field (the full sequenced trimmed to the V(D)J region) is determined
by pairwise alignment of the original trimmed sequence to the new full-length
sequence, followed by trimming the new sequence to match.
"""

import sys
from pathlib import Path
from Bio import SeqIO, Align
from igseq.record import RecordReader, RecordWriter

ALIGNER = Align.PairwiseAligner()
# For the purposes of the alignment we're doing here, gaps are expensive,
# except for end gaps on the query (the original V(D)J-trimmed sequence) when
# aligning to the target (the new full-length sequence).
ALIGNER.gap_score = -10
ALIGNER.query_end_gap_score = 0

def sonar_recluster(rearr_in, rearr_out, fasta_dir):
    """Replace SONAR cluster centroid sequences with updated versions"""
    seqs = {}
    for fasta in Path(fasta_dir).glob("*.fa"):
        rec = SeqIO.read(fasta, "fasta")
        seqs[fasta.stem] = str(rec.seq)
    with RecordReader(rearr_in) as reader, RecordWriter(rearr_out) as writer:
        for row in reader:
            if row["cluster_count"]:
                row["sequence"] = seqs[row["sequence_id"]]
                # target is the new full-length sequence, and query the
                # previous trimmed sequence
                aln = ALIGNER.align(row["sequence"], row["sequence_alignment"])[0]
                # edges of query sequence in alignment
                idx_query_start = list(aln.indices[1]).index(0)
                idx_query_end = list(aln.indices[1]).index(len(row["sequence_alignment"])-1) + 1
                # trimming target (new full seq) and query (prev. trimmed seq)
                # according to those edges
                target_trim = aln[0][idx_query_start:idx_query_end]
                query_trim = aln[1][idx_query_start:idx_query_end]
                assert query_trim == row["sequence_alignment"] # it better be!
                row["sequence_aligment"] = target_trim
            writer.write(row)

def main():
    """CLI for sonar_recluster"""
    sonar_recluster(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
