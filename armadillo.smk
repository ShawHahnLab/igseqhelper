rule armadillo:
    """General ARMADiLLO run"""
    output: expand("analysis/armadillo/{{category}}.{{locus_chain}}/{{case}}/{{case}}.{suffix}", suffix=["tiles.html", "ARMADiLLO.html", "freq_table.html", "ARMADiLLO.Detailed.text", "freq_table.txt", "ARMADiLLO.fasta"])
    input:
        seq_from="analysis/armadillo/{category}.{locus_chain}/{case}/from.fasta",
        seq_to="analysis/armadillo/{category}.{locus_chain}/{case}/{case}.fasta"
    log:
        main="analysis/armadillo/{category}.{locus_chain}/{case}/armadillo.log.txt",
        conda="analysis/armadillo/{category}.{locus_chain}/{case}/armadillo.conda_build.txt"
    threads: 8
    params:
        max_iter=10000,
        random_seed=123,
        species="rhesus",
        outdir="analysis/armadillo/{category}.{locus_chain}/{case}"
    conda: "armadillo.yml"
    shadow: "shallow"
    shell:
        """
            seqid_uca=$(grep "^>" {input.seq_from} | cut -c 2- | cut -f 1 -d ' ')
            seqid_seq=$(grep "^>" {input.seq_to} | cut -c 2- | cut -f 1 -d ' ')
            conda list --explicit > {log.conda}
            (
              date +'%Y-%m-%d %H:%M:%S %Z'
              echo "Threads: {threads}"
              echo "Max iter: {params.max_iter}"
              echo "Random seed: {params.random_seed}"
              echo "Species: {params.species}"
              echo "Chain: {wildcards.locus_chain}"
              echo "UCA Seq ID: $seqid_uca"
              echo "Antibody Seq ID: $seqid_seq"
              ARMADiLLO -species rhesus -chain {wildcards.locus_chain} \
                  -threads {threads} \
                  -fulloutput \
                  -max_iter {params.max_iter} -random_seed {params.random_seed} \
                  -uca {input.seq_from} -seq {input.seq_to} 2>&1
              date +'%Y-%m-%d %H:%M:%S %Z'
            ) | tee {log.main}
            if grep -q -x "0 mutations found" {log.main}; then
                touch {output}
            else
                for suffix in tiles.html ARMADiLLO.html freq_table.html ARMADiLLO.Detailed.text freq_table.txt ARMADiLLO.fasta; do
                    mv "$seqid_seq.$suffix" "{params.outdir}/{wildcards.case}.$suffix"
                done
                mv *.css "{params.outdir}/"
            fi
        """

def input_for_armadillo_input_ancs(w):
    # I use keywords gamma/kappa/lambda (implicit IgG+ source material for
    # this workflow) while ARMADiLLO uses heavy/kappa/lambda (so, the locus, no
    # isotype assumed)
    chain_type = "gamma" if w.locus_chain == "heavy" else w.locus_chain
    subject = ANTIBODY_LINEAGES[w.antibody_lineage]["Subject"]
    return (f"analysis/sonar/{subject}.{chain_type}/longitudinal-custom-{w.antibody_lineage}"
        f"/output/sequences/nucleotide/longitudinal-custom-{w.antibody_lineage}_inferredAncestors.fa")

rule armadillo_input_ancs:
    """ARMADiLLO input: pair of inferred ancestors leading toward a particular isolate in the tree"""
    output:
        seq_from="analysis/armadillo/{antibody_lineage}.{antibody_isolate}.{locus_chain}/{anc1}.{anc2}/from.fasta",
        seq_to="analysis/armadillo/{antibody_lineage}.{antibody_isolate}.{locus_chain}/{anc1}.{anc2}/{anc1}.{anc2}.fasta"
    input: input_for_armadillo_input_ancs
    run:
        with open(output.seq_from, "w") as f_from, open(output.seq_to, "w") as f_to:
            for rec in SeqIO.parse(input[0], "fasta"):
                anc, tips, length = rec.id.split(";", 3)
                tips = tips.split(",")
                if wildcards.antibody_isolate in tips:
                    anc = f"anc{anc}"
                    seqid = ".".join([
                        wildcards.antibody_lineage,
                        wildcards.antibody_isolate,
                        wildcards.locus_chain,
                        anc])
                    if anc == wildcards.anc1:
                        f_from.write(f">{seqid}\n{rec.seq}\n")
                    elif anc == wildcards.anc2:
                        f_to.write(f">{seqid}\n{rec.seq}\n")
