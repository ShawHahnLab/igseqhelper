"""
Rules for select final output files per-subject-per-lineage.

Everything here is just copying over lightweight files created in other rules.
"""

def set_chain_type(w):
    """Set the chain_type wildcard based on chain and lineage and return as dict"""
    w = dict(w)
    for lineage, attrs in ANTIBODY_LINEAGES.items():
        if w["antibody_lineage"] == lineage:
            break
    else:
        raise ValueError(f"Can't find lineage {antibody_lineage}")
    w["chain_type"] = "gamma"
    if w["chain"] != "heavy":
        w["chain_type"] = {"IGK": "kappa", "IGL": "lambda"}[attrs["LightLocus"]]
    return w

def summary_setup_lineage_helper_rules():
    """Define per-lineage summary rules dynamically"""
    for lineage, attrs in ANTIBODY_LINEAGES.items():
        subject = attrs["Subject"]
        targets = []
        for chain in ["heavy", "light"]:
            targets_chain = []
            w = set_chain_type({"subject": attrs["Subject"], "antibody_lineage": lineage, "chain": chain})
            # will default to just igphyml-based outputs for automatic alignment,
            # and only expect the custom alignment outputs if a custom alignment is
            # already in place
            things = [""]
            custom_aln = Path(expand(
                "analysis/sonar/{subject}.{chain_type}/alignment.{antibody_lineage}.fa", **w)[0])
            if custom_aln.exists():
                things += [".custom"]
            targets_chain += expand(
                "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_igphyml{thing}.{ext}",
                thing=things, ext=["tree", "tree.pdf"], **w)
            targets_chain += expand(
                "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_aligned{thing}.afa",
                thing=things, **w)
            targets_chain += expand(
                "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_collected.fa", **w)
            targets_chain += expand(
                "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_inferredAncestors{thing}.fa",
                thing=things, **w)
            targets_chain += expand(
                "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_inferredAncestors{thing}.common.fa",
                thing=things, **w)
            targets += targets_chain
            rule:
                f"Summary outputs for subject {subject} lineage {lineage} chain {chain}"
                name: f"summary_lineage_{lineage}_{chain}"
                input: targets_chain
        rule:
            f"Summary outputs for subject {subject} lineage {lineage}"
            name: f"summary_lineage_{lineage}"
            input: targets + [f"summary/{subject}/{lineage}/{lineage}_divergence.pdf"]

def summary_setup_subject_helper_rules():
    """Define per-subject summary rules dynamically"""
    # note which subjects have IgM reads sequenced for which chain
    m_subjects = {"mu": set(), "kappa": set(), "lambda": set()}
    for sample, attrs in SAMPLES.items():
        if attrs["Run"] and "IgM" in attrs["SpecimenAttrs"]["CellType"]:
            m_subjects[attrs["Type"]].add(attrs["SpecimenAttrs"]["Subject"])
    # organize by locus rather than chain type
    m_subjects["IGH"] = m_subjects.pop("mu")
    m_subjects["IGK"] = m_subjects.pop("kappa")
    m_subjects["IGL"] = m_subjects.pop("lambda")
    for subject in {attrs["Subject"] for attrs in SPECIMENS.values()}:
        targets = []
        # MINING-D output from heavy chain IgM
        if subject in m_subjects["IGH"]:
            targets += expand(
                "summary/{subject}/{subject}.miningd.{pval}.txt",
                subject=subject, pval=["default", "sensitive"])
        # IgDiscover output from IgM for whichever loci were sequenced
        loci = {locus for locus in m_subjects.keys() if subject in m_subjects[locus]}
        for locus in loci:
            segments = ["V", "J"]
            if locus == "IGH":
                segments += "D"
            targets += expand(
                "summary/{subject}/{subject}.germline.{locus}{segment}.fasta",
                subject=subject, locus=locus, segment=segments)
        rule:
            f"Summary outputs for subject {subject}"
            name: f"summary_subject_{subject}"
            input: targets

summary_setup_lineage_helper_rules()
summary_setup_subject_helper_rules()

### Antibody Lineage summary rules

rule summary_tree:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_igphyml.{ext}"
    input: lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/longitudinal-{antibody_lineage}_igphyml.{ext}", **set_chain_type(w))
    shell: "cp {input} {output}"

rule summary_tree_custom:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_igphyml.custom.{ext}"
    input: lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal-custom-{antibody_lineage}/output/longitudinal-custom-{antibody_lineage}_igphyml.{ext}", **set_chain_type(w))
    shell: "cp {input} {output}"

# Using `igseq convert` for FASTA files to ensure we get unwrapped versions
rule summary_aln_auto:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_aligned.afa"
    input: lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/work/phylo/longitudinal-{antibody_lineage}_aligned.afa", **set_chain_type(w))
    shell: "igseq convert {input} {output}"
    
rule summary_aln_custom:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_aligned.custom.afa"
    input: lambda w: expand("analysis/sonar/{subject}.{chain_type}/alignment.{antibody_lineage}.fa", **set_chain_type(w))
    shell: "igseq convert {input} {output}"

rule summary_collected:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_collected.fa"
    input: lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/sequences/nucleotide/longitudinal-{antibody_lineage}-collected.fa", **set_chain_type(w))
    shell: "igseq convert {input} {output}"

rule summary_ancestors_auto:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_inferredAncestors.fa"
    input: lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal-{antibody_lineage}/output/sequences/nucleotide/longitudinal-{antibody_lineage}_inferredAncestors.fa", **set_chain_type(w))
    shell: "igseq convert {input} {output}"

rule summary_ancestors_custom:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_inferredAncestors.custom.fa"
    input: lambda w: expand("analysis/sonar/{subject}.{chain_type}/longitudinal-custom-{antibody_lineage}/output/sequences/nucleotide/longitudinal-custom-{antibody_lineage}_inferredAncestors.fa", **set_chain_type(w))
    shell: "igseq convert {input} {output}"

rule summary_ancestors_common:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_inferredAncestors.common.fa"
    input: lambda w: expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.common.fa", **set_chain_type(w))
    shell: "igseq convert {input} {output}"

rule summary_ancestors_custom_common:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_{chain}_inferredAncestors.custom.common.fa"
    input: lambda w: expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.custom.common.fa", **set_chain_type(w))
    shell: "igseq convert {input} {output}"

rule summary_germline_divergence_plot:
    output: "summary/{subject}/{antibody_lineage}/{antibody_lineage}_divergence.pdf"
    input: "analysis/reporting/by-lineage/{antibody_lineage}.divergence.pdf"
    shell: "cp {input} {output}"

### Subject summary rules

rule summary_miningd:
    output: "summary/{subject}/{subject}.miningd.{thing}.txt"
    input: "analysis/mining-d/{subject}.output.{thing}.txt"
    shell: "cp {input} {output}"

rule summary_germline:
    output: "summary/{subject}/{subject}.germline.{locus}{segment}.fasta"
    input: "analysis/germline/{subject}.{locus}/{segment}.fasta"
    shell: "igseq convert {input} {output}"
