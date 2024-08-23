"""
Rules for running MINING-D using deduplicated CDRH3 sequences from SONAR's clustered IgM+ reads.
"""

import gzip
from collections import defaultdict

# The MINING-D authors actually clustered the CDR3 sequences *themselves*
# whereas I'm just deduplicating the CDR3s from SONAR's antibody cluster
# centroids.  Possibly clustering their way would be even better; not sure.
def cdrh3s_from_airr(tsvgzs_in, fasta_out):
    seqs = defaultdict(int)
    for path in tsvgzs_in:
        opener = open
        if str(path).endswith(".gz"):
            opener = gzip.open
        with opener(path, "rt") as f_in:
            reader = DictReader(f_in, delimiter="\t")
            for row in reader:
                # If this is SONAR AIRR only use the cluster centroids
                if "cluster_count" in row and not row["cluster_count"]:
                    continue
                # IgBLAST puts both cdr3 and junction cols in its AIRR TSV
                # output while SONAR only provides junction.  We'll work
                # with either.
                if row["locus"] == "IGH" and (row.get("cdr3") or row.get("junction")):
                    try:
                        seq = row["cdr3"]
                    except KeyError:
                        seq = row["junction"][3:-3]
                    seqs[seq] += 1
    # Write deduplicated sequences sorted by abundance
    seqs = list(seqs.items())
    seqs.sort(key=lambda pair: -pair[1])
    with open(fasta_out, "wt") as f_out:
        for idx, pair in enumerate(seqs):
            seqid = f"cdr3_{idx}"
            desc = f"count={pair[1]}"
            seq = pair[0]
            f_out.write(f">{seqid} {desc}\n{seq}\n")

def miningd_sonar_input(w):
    specs = []
    for specimen in SPECIMENS.values():
        if specimen["Subject"] == w.subject and "IgM" in specimen["CellType"]:
            specs.append(specimen["Specimen"])
    return expand(
        "analysis/sonar/{subject}.mu/{specimen}/output/tables/{specimen}_rearrangements.tsv",
        subject=w.subject, specimen=specs)

rule miningd_get_cdrh3s:
    output: "analysis/mining-d/{subject}.cdr3.fasta"
    input: miningd_sonar_input
    run:
        cdrh3s_from_airr(input, output[0])

# Note that the order and naming of sequences changes from run to run, but the
# sequences themselves do appear to typically be the same for any given input.
# The miningd_final rule below should give consistent output.
# Note that the sequences are not *always* the same, though-- there's very
# slight randomness in the output in some cases.
rule miningd_run:
    output: "analysis/mining-d/{subject}.output.{pval}.fasta"
    input: "analysis/mining-d/{subject}.cdr3.fasta"
    log:
        main="analysis/mining-d/{subject}.{pval}.log.txt",
        conda="analysis/mining-d/{subject}.{pval}.conda_build.txt"
    conda: "envs/mining-d.yaml"
    threads: 8
    params:
        # The authors thought 600 was a good choice for rhesus macaque and
        # human
        num_k_mers=600,
        # 4.5e-36 is The default setting.  Will tend to leave off a few
        # nucleotides at the edges, but will avoid false positives.
        # (If a wildcard value is given that's not "default" or "sensitive",
        # take that text itself as the P value threshold.)
        p_val_th=lambda w: {"default": "4.5e-36", "sensitive": "4.5e-16"}.get(w.pval, w.pval)
    shell:
        """
            conda list --explicit > {log.conda}
            numseqs_in=$(grep -c "^>" {input} || echo 0)
            echo "start: $(date +'%Y-%m-%d %H:%M:%S %Z')" >> {log.main}
            echo "threads: {threads}" >> {log.main}
            echo "input: {input}" >> {log.main}
            echo "input seqs: $numseqs_in" >> {log.main}
            echo "num_k_mers: {params.num_k_mers}" >> {log.main}
            echo "p_val_th: {params.p_val_th}" >> {log.main}
            python MINING-D/MINING_D.py -t {threads} \
                -n {params.num_k_mers} \
                -p {params.p_val_th} \
                -i {input} -o {output}
            numseqs_out=$(grep -c "^>" {output} || echo 0)
            echo "output seqs: $numseqs_out" >> {log.main}
            echo "end: $(date +'%Y-%m-%d %H:%M:%S %Z')" >> {log.main}
        """

# Sort the sequences and write as plain text, one sequence per line, to get the
# same output for the same input.
rule miningd_final:
    output: "analysis/mining-d/{subject}.output.{pval}.txt"
    input: "analysis/mining-d/{subject}.output.{pval}.fasta"
    run:
        with open(input[0]) as f_in, open(output[0], "w") as f_out:
            seqs = [str(record.seq) for record in SeqIO.parse(f_in, "fasta")]
            seqs.sort()
            for seq in seqs:
                f_out.write(f"{seq}\n")
