"""
Rules for reporting on MINING-D output.
"""

rule report_miningd_combo_tree:
    """Make a tree of all subjects' MINING-D output compared with all known D sequences."""
    output: "analysis/reporting/mining-d/all.nex"
    input: "analysis/reporting/mining-d/all.msa.fasta"
    log:
        conda="analysis/reporting/mining-d/all.nex.conda_build.txt"
    conda: "envs/igseq.yaml"
    shell:
        """
            conda list --explicit > {log.conda}
            igseq tree {input} {output}
        """

rule report_miningd_combo_msa:
    """Align all subjects' MINING-D output with all known D sequences."""
    output: "analysis/reporting/mining-d/all.msa.fasta"
    input:
        subject="analysis/reporting/mining-d/all.fasta",
        refs="analysis/reporting/mining-d/refs.fa"
    log:
        conda="analysis/reporting/mining-d/all.msa.fasta.conda_build.txt"
    conda: "envs/igseq.yaml"
    shell:
        """
            conda list --explicit > {log.conda}
            cat {input.refs} {input.subject} | igseq msa --input-format fa - {output}
        """

rule report_miningd_combo_fasta:
    """Combine all subjects' MINING-D outputs into one FASTA."""
    output: "analysis/reporting/mining-d/all.fasta"
    input: TARGET_REPORT_MININGD_FASTAS
    params:
        subjects = AVAILABLE_MININGD["subject"]
    run:
        with open(output[0], "wt") as f_out:
            for subject, path in zip(params.subjects, input):
                with open(path) as f_in:
                    for rec in SeqIO.parse(f_in, "fasta"):
                        rec.id = f"{subject}_{rec.id}"
                        SeqIO.write(rec, f_out, "fasta-2line")

rule report_miningd_tree:
    """Make a tree of each subject's MINING-D output compared with all known D sequences."""
    output: "analysis/reporting/mining-d/{subject}/{subject}.nex"
    input:
        msa="analysis/reporting/mining-d/{subject}/{subject}.msa.fasta",
        subject="analysis/reporting/mining-d/{subject}/{subject}.fasta",
        refs="analysis/reporting/mining-d/refs.fa"
    log:
        conda="analysis/reporting/mining-d/{subject}/{subject}.nex.conda_build.txt"
    conda: "envs/igseq.yaml"
    shell:
        """
            conda list --explicit > {log.conda}
            igseq tree \
                -C subject=#880000 -C refs=#000000 \
                -L subject=<(grep '^>' {input.subject} | cut -c 2- | cut -f 1 -d ' ') \
                -L refs=<(grep '^>' {input.refs} | cut -c 2- | cut -f 1 -d ' ') \
                {input.msa} {output}
        """

rule report_miningd_msa:
    """Align a subject's MINING-D output with all known D sequences."""
    output: "analysis/reporting/mining-d/{subject}/{subject}.msa.fasta"
    input:
        subject="analysis/reporting/mining-d/{subject}/{subject}.fasta",
        refs="analysis/reporting/mining-d/refs.fa"
    log:
        conda="analysis/reporting/mining-d/{subject}/{subject}.msa.fasta.conda_build.txt"
    conda: "envs/igseq.yaml"
    shell:
        """
            conda list --explicit > {log.conda}
            cat {input.refs} {input.subject} | igseq msa --input-format fa - {output}
        """

rule report_miningd_refs:
    """
    Make a combined D gene FASTA across rhesus refs.

    Output has one line per unqiue sequence across all references, with
    automated seq IDs, and IDs within individual references in the
    descriptions.
    """
    output: "analysis/reporting/mining-d/refs.fasta"
    run:
        from igseq.util import FILES, DATA
        from base64 import b32encode
        from hashlib import sha1
        refs = {
            "germ/rhesus/imgt/IGH/IGHD.fasta": "IMGT",
            "germ/rhesus/sonarramesh/IGH/IGHD.fasta": "Ramesh",
            "germ/rhesus/kimdb/IGH/IGHD.fasta": "KIMDB"}
        ref_seqs = defaultdict(list)
        for path in FILES:
            key = str(path.relative_to(DATA))
            if key not in refs:
                continue
            with open(path) as f_in:
                for rec in SeqIO.parse(f_in, "fasta"):
                    ref_seqs[str(rec.seq)].append((refs[key], rec.id))
        with open(output[0], "wt") as f_out:
            for seq, id_list in ref_seqs.items():
                # ID syntax inspired by 10.1016/j.immuno.2023.100025
                # except I'm making it deterministic (why make it random when you
                # can make it reproducible?)
                seqhash = sha1()
                seqhash.update(seq.encode("UTF8"))
                seqhash = b32encode(seqhash.digest()).decode("ASCII")
                seqid = "IGHD0-" + seqhash[:4] + "*00"
                seqdesc = " ".join([f"{ref}:{seqid}" for ref, seqid in id_list])
                f_out.write(f">{seqid} {seqdesc}\n{seq}\n")

rule report_miningd_output:
    output: "analysis/reporting/mining-d/{subject}/{subject}.fasta"
    input: "analysis/mining-d/{subject}.output.default.fasta"
    shell: "cp {input} {output}"
