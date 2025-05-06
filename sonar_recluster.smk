"""
Refine SONAR cluster centroid sequences using raw FASTQ data.

SONAR's method to define a sequence for each cluster centroid is to take the
most abundant unique sequence per cluster.  That keeps it reasonably quick and
simple for whole-repertoire clustering, but we can make more accurate consensus
sequences for subsets of a repertoire (like relevant sequences for members of a
lineage) with the original input FASTQ given to SONAR.
"""

# TODO merge this into sonar.smk and put directly in line between module 2 and
# module 3 if all works as intended

rule sonar_recluster_get_rearrangements:
    """filter the full SONAR AIRR TSV to just those rows for this island"""
    output:
        rearr=WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/rearrangements.tsv"
    input:
        seqids=WD_SONAR/"output/tables/islandSeqs_{antibody_lineage}.txt",
        rearr=WD_SONAR/"output/tables/{specimen}_rearrangements.tsv"
    run:
        with open(input.seqids) as f_in:
            seqids = [line.strip() for line in f_in]
        with open(input.rearr) as f_in, open(output.rearr, "w") as f_out:
            reader = csv.DictReader(f_in, delimiter="\t")
            writer = csv.DictWriter(f_out, reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                if row["centroid"] in seqids:
                    writer.writerow(row)

rule sonar_recluster_get_fastqs:
    output: directory(WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/fastq")
    input:
        fqgzdir="analysis/samples-by-specimen/{specimen}.{chain_type}",
        rearr_island=WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/rearrangements.tsv"
    shell:
        """
            sonar_cluster_fastqs.py "{input.rearr_island}" "{output}" "{input.fqgzdir}"
        """

rule sonar_recluster_fastq_consensus:
    output: directory(WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/cons")
    input: WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/fastq"
    shell:
        """
            for fqgz in {input}/*.fastq.gz; do
                fasta={output}/$(basename "$fqgz" | sed s/fastq.gz/fa/)
                fastq_consensus.py "$fqgz" "$fasta"
            done
        """

rule sonar_recluster:
    """Update each SONAR cluster centroid sequence with quality-aware consensus"""
    output: WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/rearrangements_reclustered.tsv"
    input:
        cons=WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/cons",
        rearr_island=WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/rearrangements.tsv"
    shell: "sonar_recluster.py {input.rearr_island} {output} {input.cons}"

rule sonar_recluster_collapse:
    """Collapse duplicated centroid sequences after updating"""
    output: WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/rearrangements_reclustered_collapsed.tsv"
    input: WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/rearrangements_reclustered.tsv"
    shell: "sonar_recluster_collapse.py {input} {output}"

rule sonar_recluster_fasta:
    """Refined version of SONAR module 2's FASTA output"""
    output: WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}.fa"
    input: WD_SONAR/"output/sequences/nucleotide/islandSeqs_recluster_{antibody_lineage}/rearrangements_reclustered_collapsed.tsv"
    shell: "igseq convert --col-seq sequence_alignment {input} {output}"
