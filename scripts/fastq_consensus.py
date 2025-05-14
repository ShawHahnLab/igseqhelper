#!/usr/bin/env python

"""
Make quality-aware consensus sequences from sets of FASTQ reads
"""

from pathlib import Path
from argparse import ArgumentParser
from collections import defaultdict
from Bio import Align
from igseq.record import RecordReader, RecordWriter, RecordHandler

ALIGNER = Align.PairwiseAligner()
ALIGNER.gap_score = -10

def load_cluster_fastqs(fqgz_dir_in):
    """Load .fastq.gz files into dictionary of files->reads and list of reads"""
    # assuming files are named {centroid}.fastq.gz
    clustermap = {} # cluster ID -> list of read IDs
    reads = [] # list of read dictionaries
    for fqgz_in in Path(fqgz_dir_in).glob("*.fastq.gz"):
        key = fqgz_in.name.removesuffix(".fastq.gz")
        clustermap[key] = []
        with RecordReader(fqgz_in) as reader:
            for row in reader:
                row["sequence_quality"] = RecordHandler.decode_phred(row["sequence_quality"])
                clustermap[key].append(row["sequence_id"])
                reads.append(row)
    return clustermap, reads

def tally_scores(reads):
    """List base call quality score details per position across FASTQ reads

    Each returned list item is a dictionary giving lists of quality scores per
    observed base call, where the list indexes correspond to positions in an
    inferred reference sequence across the set of reads.

    The reference is chosen from among the supplied reads by selecting the top
    match after sorting for:

     1. Most abundant sequence
     2. Most abundant read length
     3. Longest read
     4. Sequence content
    """
    # take most abundant sequence as the reference, with secondary sorting by
    # most common length
    totals = defaultdict(int)
    lentotals = defaultdict(int)
    for row in reads:
        totals[row["sequence"]] += 1
        lentotals[len(row["sequence"])] += 1
    stats = [(totals[seq], lentotals[len(seq)], len(seq), seq) for seq in totals]
    stats.sort(reverse=True)
    ref = stats[0][3]
    # align the deduplicated sequences
    alns = {}
    for attrs in stats:
        # I don't know what's going on under the hood in PairwiseAligner, but
        # accessing some of these properties is evidently quite expensive for
        # some reason, so we'll just stash the bits we want before getting to
        # the loop below (rather than storing the whole object)
        query = attrs[3]
        aln = ALIGNER.align(ref, query)[0]
        alns[query] = {
            "aln_query": aln[1][:],
            "length": aln.length,
            "indices": [[int(num) for num in items] for items in aln.indices]}
    # index in ref -> base call -> list of Q scores
    scores = defaultdict(lambda: defaultdict(list))
    for row in reads:
        aln = alns[row["sequence"]]
        for idx in range(aln["length"]):
            idx_target = aln["indices"][0][idx] # index in original target (ref) seq
            idx_query = aln["indices"][1][idx] # index in query seq
            # -1 for gaps in target, which we don't care about, gaps in
            # read, which won't affect things anyway
            if idx_target != -1 and idx_query != -1:
                basecall = aln["aln_query"][idx]
                scores[idx_target][basecall].append(row["sequence_quality"][idx_query])
    # Dictionary keyed by index is stupid.  Switch to a plain list... but make
    # sure it all makes sense first.
    scores = sorted(scores.items())
    assert list(range(1+max(pair[0] for pair in scores))) == [pair[0] for pair in scores]
    scores = [pair[1] for pair in scores]
    return scores

def make_consensus(scores):
    """Make consensus from list of base call and Q info for a set of reads"""
    consensus = ""
    quals = []
    check = 0
    for base_scores in scores:
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

def __centroids_by_consensus(reads, clustermap, consensuses):
    # First, make a mapping of each read sequence to the cluster that contains
    # it
    # (just a helper for matching things up below)
    readmap = {row["sequence_id"]: row for row in reads}
    # (each read seq -> centroid)
    # In theory this could have clashing entries, but in practice it won't
    # since identical reads will have always been placed in the same cluster.
    read_seqs_to_centroids = {}
    for centroid, read_ids in clustermap.items():
        for read_id in read_ids:
            seq = readmap[read_id]["sequence"]
            read_seqs_to_centroids[seq] = centroid
    # Now, with the existing mapping of centroids to consensus seqs and this
    # mapping of read seqs to centroids, gather all relevant centroids for each
    # unique consensus sequence
    # (each consensus seq -> set of matching centroids via eithercons or reads)
    cons_back = defaultdict(set)
    for centroid in clustermap:
        # (the second item of each tuple is the quality scores but we don't
        # care about that here)
        cons_seq = consensuses[centroid][0]
        # This will gather up each seq -> centroid match including those for
        # identical consensus sequences between clusters
        cons_back[cons_seq].add(centroid)
        # also, if this consensus matches a read exactly, note that read's
        # centroid too
        read_centroid = read_seqs_to_centroids.get(cons_seq)
        if read_centroid:
            cons_back[cons_seq].add(read_centroid)
    return cons_back

# Making that happen hurt my brain way more than I expected.  Thankfully we
# have Stack Overflow.  Below adapted from:
# https://stackoverflow.com/a/16306257
def __group_overlaps(sets):
    result = []
    for item in sets:
        for members in result:
            if members.intersection(item):
                members.update(item)
                break
        else:
            result.append(set(item))
    # switch sorted lists for everything, so it's all deterministic
    result = [sorted(items) for items in result]
    result.sort()
    return result

def merge_clusters(reads, clustermap, consensuses):
    """Given clusters of reads and consensuses per cluster, merge overlapping clusters

    reads: full, static list of reads
    clustermap: full set of centroid IDs -> lists of read IDs
    consensuses: consensus for each cluster in clustermap

    This will return an updated version of clustermap, where any clusters whose
    consensus sequences perfectly match a read or consensus sequence of another
    cluster are merged together.
    """
    #
    # the set of cluster centroids relevant for each consensus sequence, either
    # in terms of what cluster has that consensus or what single read
    # (belonging to whatever cluster) is that sequence
    cons_back = __centroids_by_consensus(reads, clustermap, consensuses)
    # But we're not done yet!  The same centroid may show up in the context of
    # multiple distinct consensus seqs.  Next, assign a new representative for
    # each centroid based on all those it touches.
    # We don't actually care what those sequences are at this point; it's just
    # about the (potentially overlapping) sets of centroids.
    centroid_groups = __group_overlaps(cons_back.values())
    # Now we have the centroids of the previous clustered grouped together and
    # we can take the first one of each group as the ID of the new cluster (for
    # those that are grouped; singletons will remain the same) and group all
    # applicable reads.
    clustermap_new = {}
    for centroids in centroid_groups:
        clustermap_new[centroids[0]] = []
        for centroid2 in centroids:
            clustermap_new[centroids[0]] += clustermap[centroid2]
    return clustermap_new

def _write_iter_details(dir_out_iter, clustermap, consensuses, reads):
    dir_out_iter.mkdir(parents=True, exist_ok=True)
    _write_consensus_table(dir_out_iter/"out.csv", clustermap, consensuses, reads)
    for centroid, seq_ids in clustermap.items():
        with RecordWriter(dir_out_iter/f"{centroid}.consensus.fastq") as writer:
            writer.write({
                "sequence_id": centroid,
                "sequence": consensuses[centroid][0],
                "sequence_quality": RecordHandler.encode_phred(consensuses[centroid][1])})
        reads_tmp = [rec.copy() for rec in reads if rec["sequence_id"] in seq_ids]
        with RecordWriter(dir_out_iter/f"{centroid}.fastq.gz") as writer:
            for rec in reads_tmp:
                rec["sequence_quality"] = RecordHandler.encode_phred(
                    rec["sequence_quality"])
                writer.write(rec)

def _write_consensus_table(csv_out, clustermap, consensuses, reads):
    # my original idea, with # one row per original centroid:
    #
    # sequence_id_orig      SONAR original sequence ID
    # sequence_orig         SONAR original representative sequence
    # cluster_count_orig    SONAR original cluster_count
    # duplicate_count_orig  SONAR original duplicate_count
    # sequence_cons         initial consensus for the original cluster (worth having this?)
    # sequence_id           final centroid ID (choose from among IDs by abundance?)
    # sequence              final consensus sequence
    # duplicate_count       observed duplicate count for final consensus
    # cluster_count         observed read count in final cluster
    #
    # ...but time is limited.  Instead:
    #
    # sequence_id           final centroid ID
    # sequence              final consensus sequence
    # sequence_quality      final consensus sequence: inferred quality scores
    # sequence_quality_min  min of sequence_quality (as integer Q score)
    # duplicate_count       observed duplicate count for final consensus
    # cluster_count         observed read count in final cluster
    rows = []
    for centroid, seq_ids in clustermap.items():
        cons_seq, cons_qual = consensuses[centroid]
        cons_qual_txt = RecordHandler.encode_phred(cons_qual)
        reads_here = [rec for rec in reads if rec["sequence_id"] in seq_ids]
        rows.append({
            "sequence_id": centroid,
            "sequence": cons_seq,
            "sequence_quality": cons_qual_txt,
            "sequence_quality_min": min(qual for qual in cons_qual),
            "duplicate_count": sum(cons_seq == rec["sequence"] for rec in reads_here),
            "cluster_count": len(reads_here)})
    rows.sort(key=lambda row: (
        -row["cluster_count"],
        -row["sequence_quality_min"],
        -row["duplicate_count"],
        row["sequence_id"]))
    with RecordWriter(csv_out) as writer:
        for row in rows:
            writer.write(row)

def fastq_consensus_recluster(fqgz_dir_in, csv_out, dir_out=None, max_iter=10):
    """Iteratively make consensus sequences and update read clustering"""
    # 1. Load all read data from the per-cluster .fastq.gz files
    # 2. Take consensus using the most abundant read as the reference (to mirror
    #    SONAR's approach)
    # 3. Compare each consensus sequence with all reads and all other consensus
    #    sequences.  For any matches, merge the relevant clusters together for the
    #    next iteration
    # 4. Steps 2 through 4 again, but this time with the combined sets of reads for
    #    applicable clusters; for any clusters that have been merged, calculated a
    #    new consensus sequence.  Repeat until there's no more reclustering, up
    #    to some limit
    # 5. Final output is a table of per-cluster info
    if dir_out:
        dir_out = Path(dir_out)
    # Load initial per-cluster .fastq.gz files
    # it's easeiest for our purposes here to actually just store a mapping of
    # reads to clusters on the one hand and a big list of all reads on the
    # other.
    clustermap, reads = load_cluster_fastqs(fqgz_dir_in)
    clustermap_diff = clustermap.copy()
    consensuses = {}
    for cluster_idx in range(max_iter):
        print(f"{fqgz_dir_in}: iteration {cluster_idx}")
        # Make consensus sequences for any clusters that have been updated
        # since the last iteration.  On the first iteration, that'll just be
        # everything.  On later iterations the "diff" dictionary will only have
        # entries for what needs updating.
        consensuses_new = {}
        for centroid, seq_ids in clustermap_diff.items():
            consensuses_new[centroid] = make_consensus(tally_scores(
                [rec for rec in reads if rec["sequence_id"] in seq_ids]))
        print(f"{fqgz_dir_in}: iteration {cluster_idx}: "
                f"{len(consensuses_new)} newly-created consensuses")
        # Even if we've brought more reads in for a particular cluster the
        # consensus might be the same as before.
        actually_new = sum(item != consensuses.get(c) for c, item in consensuses_new.items())
        print(f"{fqgz_dir_in}: iteration {cluster_idx}: ({actually_new} actually different)")
        consensuses.update(consensuses_new)
        print(f"{fqgz_dir_in}: iteration {cluster_idx}: {len(consensuses)} total consensuses")
        # Find any matches between the latest consensus sequences on the one
        # hand, and all the reads and consensus sequences on the other.  Merge
        # any clusters with matches.
        print(f"{fqgz_dir_in}: iteration {cluster_idx}: {len(clustermap)} previous clusters")
        clustermap_new = merge_clusters(reads, clustermap, consensuses)
        print(f"{fqgz_dir_in}: iteration {cluster_idx}: {len(clustermap_new)} total clusters now")
        # at this point we have an updated mapping of what reads go with what
        # cluster, and an updated set of consensus sequences for each cluster.
        # If an output directory was specified, write the final (updated
        # compared to the starting point) cluster info for this iteration.
        if dir_out:
            _write_iter_details(
                dir_out/f"iter{cluster_idx:03}", clustermap_new, consensuses, reads)
        # figure out which clusters were updated since the last iteration.  If
        # none, we're done.
        clustermap_diff = {}
        for centroid, reads_new in clustermap_new.items():
            if clustermap[centroid] != reads_new:
                clustermap_diff[centroid] = reads_new[:]
        print(f"{fqgz_dir_in}: iteration {cluster_idx}: {len(clustermap_diff)} updated clusters")
        clustermap = clustermap_new
        # (we might have consensus sequences left over from clusters that are
        # merged into other clusters)
        consensuses = {key: consensuses[key] for key in clustermap}
        if not clustermap_diff:
            break
    # at this point clustermap and consensuses should have the final info
    _write_consensus_table(csv_out, clustermap, consensuses, reads)

def main():
    """CLI for fastq_consensus"""
    parser = ArgumentParser()
    addarg = parser.add_argument
    addarg("input", help="FASTQ input file (.fastq or .fastq.gz)")
    addarg("output", help="Output CSV with final consensus info")
    addarg("-O", "--output-dir", help="Optional output directory for detailed intermediate info")
    addarg("-M", "--max-iterations", default=10, help="Maximum number of cluster-merging cycles")
    args = parser.parse_args()
    fastq_consensus_recluster(args.input, args.output, args.output_dir)

if __name__ == "__main__":
    main()
