#!/usr/bin/env python

"""
Get V(D)J trimmed FASTA for updated SONAR cluster sequences
"""

import sys
from Bio import Align
from igseq.record import RecordReader, RecordWriter

ALIGNER = Align.PairwiseAligner()
# For the purposes of the alignment we're doing here, gaps are expensive,
# except for end gaps on the query (the original V(D)J-trimmed sequence) when
# aligning to the target (the new full-length sequence).
ALIGNER.gap_score = -10
ALIGNER.query_end_gap_score = 0

def sonar_recluster_fasta(csv_in, rearr_in, fasta_out):
    """Get V(D)J trimmed FASTA for updated SONAR cluster sequences"""
    with RecordReader(rearr_in) as reader:
        airr = {row["sequence_id"]: row for row in reader}
    with RecordReader(csv_in) as reader, RecordWriter(fasta_out) as writer:
        for row in reader:
            seq_old_vdj = airr[row["sequence_id"]]["sequence_alignment"]
            # target is the new full-length sequence, and query the
            # previous trimmed sequence
            aln = ALIGNER.align(row["sequence"], seq_old_vdj)[0]
            # edges of query sequence in alignment
            idx_query_start = list(aln.indices[1]).index(0)
            idx_query_end = list(aln.indices[1]).index(len(seq_old_vdj)-1) + 1
            # trimming target (new full seq) and query (prev. trimmed seq)
            # according to those edges
            seq_new_vdj = aln[0][idx_query_start:idx_query_end]
            query_trim = aln[1][idx_query_start:idx_query_end]
            assert query_trim == seq_old_vdj # it better be!
            fields = {
                "duplicate_count": row["duplicate_count"],
                "cluster_count": row["cluster_count"]}
            desc = " ".join(f"{key}={val}" for key, val in fields.items())
            writer.write({
                "sequence_id": row["sequence_id"],
                "sequence_description": desc,
                "sequence": seq_new_vdj})

def main():
    """CLI for sonar_recluster_fasta"""
    sonar_recluster_fasta(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()
