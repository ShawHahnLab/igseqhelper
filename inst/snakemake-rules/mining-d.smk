"""
Rules for running MINING-D using CDRH3 sequences from IgBLAST output for our
demultiplexed, adapter-trimmed, merged IgM+ reads.
"""

import gzip

def cdrh3s_from_airr(tsvgzs_in, fasta_out):
    with open(fasta_out, "wt") as f_out:
        for path in tsvgzs_in:
            with gzip.open(path, "rt") as f_in:
                reader = DictReader(f_in, delimiter="\t")
                for row in reader:
                    if row["locus"] == "IGH" and row["cdr3"]:
                        seq = row["cdr3"]
                        seqid = row["sequence_id"]
                        desc = " productive=" + row["productive"]
                        f_out.write(f">{seqid} {desc}\n{seq}\n")

rule miningd_get_cdrh3s:
    # using this pedantic level of detail so I can use grouped_samples_input
    # from the by-run rules, even though it's always going to be by-subject and
    # mu sequences from IgM+ cells.
    output: "analysis/mining-d/{thing}_{name}.{celltype}.{type}.cdr3.fasta"
    input: lambda w: grouped_samples_input(w, "analysis/igblast/{runid}/{samp}.tsv.gz")
    run:
        cdrh3s_from_airr(input, output[0])

rule miningd_run:
    output: "analysis/mining-d/{subject}.output.fasta" 
    input: "analysis/mining-d/subject_{subject}.igm.mu.cdr3.fasta"
    conda: str(BASEDIR/"mining-d.yml")
    threads: 8
    params:
        # The authors thought 600 was a good choice for rhesus macaque and
        # human
        num_k_mers=600,
        # 4.5e-36 is The default setting.  Will tend to leave off a few
        # nucleotides at the edges, but will avoid false positives.
        p_val_th="4.5e-36"
    shell: "python MINING-D/MINING_D.py -t {threads} -n {params.num_k_mers} -p {params.p_val_th} -i {input} -o {output}"
