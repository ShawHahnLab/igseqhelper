"""
Rules for running ARMADiLLO to infer mutation probabilities, mainly between
inferred ancestors.

https://armadillo.dhvi.duke.edu/info
"""

from functools import partial

def parse_anc_id(txt):
    """Parse SONAR/IgPhyML ancestor sequence ID"""
    anc, tips, length = txt.split(";", 3)
    tips = tips.split(",")
    return anc, tips, float(length)

def input_for_armadillo_by_ancs(w, antibody_lineage, antibody_isolate, locus_chain):
    subject = ANTIBODY_LINEAGES[antibody_lineage]["Subject"]
    chain_type = "gamma" if locus_chain == "heavy" else locus_chain
    path_ancs = checkpoints.sonar_module_3_igphyml_custom.get(
        subject=subject,
        chain_type=chain_type,
        antibody_lineage=antibody_lineage).output.inferred_nucl
    ancs = []
    for rec in SeqIO.parse(path_ancs, "fasta"):
        anc, tips, _ = parse_anc_id(rec.id)
        if antibody_isolate in tips:
            ancs.append(anc)
    ancs = ancs[::-1]
    cases = []
    for idx, anc_from in enumerate(ancs):
        for anc_to in ancs[(idx+1):]:
            cases.append(f"anc{anc_from}.anc{anc_to}")
    return expand(
        "analysis/armadillo/ancs.{antibody_lineage}.{antibody_isolate}.{locus_chain}/{case}/{case}.tiles.html",
        antibody_lineage=antibody_lineage,
        antibody_isolate=antibody_isolate,
        locus_chain=locus_chain,
        case=cases)

def input_for_armadillo_by_isolate(w, antibody_lineage, antibody_isolate, locus_chain):
    subject = ANTIBODY_LINEAGES[antibody_lineage]["Subject"]
    chain_type = "gamma" if locus_chain == "heavy" else locus_chain
    cases = f"anc1.{antibody_isolate}"
    return expand(
        "analysis/armadillo/isolate.{antibody_lineage}.{antibody_isolate}.{locus_chain}/{case}/{case}.tiles.html",
        antibody_lineage=antibody_lineage,
        antibody_isolate=antibody_isolate,
        locus_chain=locus_chain,
        case=cases)

def armadillo_setup_helper_rules():
    for antibody_isolate, attrs in ANTIBODY_ISOLATES.items():
        for chain in ["heavy", "light"]:
            seq_col = chain.capitalize() + "Seq"
            if attrs["IncludeInTracing"] == "Y" and attrs[seq_col]:
                antibody_lineage = attrs["AntibodyLineage"]
                locus_chain = "heavy"
                if chain == "light":
                    locus_chain = {"IGL": "lambda", "IGK": "kappa"}[
                        attrs["AntibodyLineageAttrs"]["LightLocus"]]
                chain_suffix = "" if chain == "heavy" else f" ({locus_chain})"
                # whole bunch of ancX -> ancY pairs
                rule:
                    f"ARMADiLLO between inferred ancestors for lineage {antibody_lineage} heading to {antibody_isolate}, {locus_chain} chain{chain_suffix}"
                    name: f"armadillo_ancs_lineage_{antibody_lineage}_isolate_{antibody_isolate}_{chain}"
                    input: partial(input_for_armadillo_by_ancs, antibody_lineage=antibody_lineage, antibody_isolate=antibody_isolate, locus_chain=locus_chain)
                # just anc1 (UCA) to each specific isolate
                rule:
                    f"ARMADiLLO between inferred UCA for lineage {antibody_lineage} and isolate {antibody_isolate}, {chain} chain{chain_suffix}"
                    name: f"armadillo_isolate_lineage_{antibody_lineage}_isolate_{antibody_isolate}_{chain}"
                    input: partial(input_for_armadillo_by_isolate, antibody_lineage=antibody_lineage, antibody_isolate=antibody_isolate, locus_chain=locus_chain)

armadillo_setup_helper_rules()

### ARMADiLLO itself

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
    conda: "envs/armadillo.yaml"
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

### ARMADiLLO - layout for an earlier inferred ancestor versus a later

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
        seq_from="analysis/armadillo/ancs.{antibody_lineage}.{antibody_isolate}.{locus_chain}/{anc1}.{anc2}/from.fasta",
        seq_to="analysis/armadillo/ancs.{antibody_lineage}.{antibody_isolate}.{locus_chain}/{anc1}.{anc2}/{anc1}.{anc2}.fasta"
    input: input_for_armadillo_input_ancs
    run:
        with open(output.seq_from, "w") as f_from, open(output.seq_to, "w") as f_to:
            for rec in SeqIO.parse(input[0], "fasta"):
                anc, tips, _ = parse_anc_id(rec.id)
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

### ARMADiLLO - layout for a UCA (ancestor 1) to a particular isolate

def input_for_armadillo_input_isolate(w):
    chain_type = "gamma" if w.locus_chain == "heavy" else w.locus_chain
    subject = ANTIBODY_LINEAGES[w.antibody_lineage]["Subject"]
    return (f"analysis/sonar/{subject}.{chain_type}/longitudinal-custom-{w.antibody_lineage}"
        f"/output/sequences/nucleotide/longitudinal-custom-{w.antibody_lineage}_inferredAncestors.fa")

rule armadillo_input_isolate:
    """ARMADiLLO input: earliest inferred ancestor and a specific mature isolate"""
    output:
        seq_from="analysis/armadillo/isolate.{antibody_lineage}.{antibody_isolate}.{locus_chain}/anc1.{antibody_isolate}/from.fasta",
        seq_to="analysis/armadillo/isolate.{antibody_lineage}.{antibody_isolate}.{locus_chain}/anc1.{antibody_isolate}/anc1.{antibody_isolate}.fasta"
    input: input_for_armadillo_input_isolate
    run:
        with open(output.seq_from, "w") as f_from:
            for rec in SeqIO.parse(input[0], "fasta"):
                anc, tips, _ = parse_anc_id(rec.id)
                if anc == "1" and wildcards.antibody_isolate in tips:
                    anc = f"anc{anc}"
                    seqid = ".".join([
                        wildcards.antibody_lineage,
                        wildcards.antibody_isolate,
                        wildcards.locus_chain,
                        anc])
                    f_from.write(f">{seqid}\n{rec.seq}\n")
        isol = ANTIBODY_ISOLATES[wildcards.antibody_isolate]
        seqid = isol["AntibodyIsolate"]
        seq = isol["HeavySeq" if wildcards.locus_chain == "heavy" else "LightSeq"]
        with open(output.seq_to, "w") as f_to:
            f_to.write(f">{seqid}\n{seq}\n")
