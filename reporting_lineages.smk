def report_lineages_divergence_input(w):
    loci = ["IGH"]
    for attrs in ANTIBODY_LINEAGES.values():
        if attrs["Lineage"] == w.antibody_lineage:
            loci.append(attrs["LightLocus"])
            subject = attrs["Subject"]
            break
    else:
        raise ValueError(f"Can't find lineage {w.antibody_lineage}")
    targets = {
        "members": expand(
            "analysis/reporting/sonar/{{antibody_lineage}}.{chain_type}/igphyml_collected.csv",
            chain_type = [{"IGH": "gamma", "IGK": "kappa", "IGL": "lambda"}[loc] for loc in loci]),
        "mabs": "analysis/reporting/by-lineage/{antibody_lineage}.mabs.csv",
        "refs": expand("analysis/germline/{subject}.{locus}/{segment}.fasta",
            subject=subject, locus=loci, segment=["V", "D", "J"])}
    return targets

def report_lineages_divergence_param_refs(w, input):
    # to condense V/D/J to just parent dirs
    return sorted(list({Path(path).parent for path in input.refs}))

rule report_lineages_divergence_plot:
    output: "analysis/reporting/by-lineage/{antibody_lineage}.divergence.pdf"
    input: "analysis/reporting/by-lineage/{antibody_lineage}.divergence.csv"
    shell: "germ_div.py -Q {input} -o {output}"

rule report_lineages_divergence:
    output: "analysis/reporting/by-lineage/{antibody_lineage}.divergence.csv"
    input: unpack(report_lineages_divergence_input)
    params:
        refs=report_lineages_divergence_param_refs
    shell: "germ_div.py -S rhesus -r {params.refs} -Q {input.members} isol={input.mabs} -G Bulk isol=Isolate -o {output}"

rule report_lineages_mabs:
    output: temp("analysis/reporting/by-lineage/{antibody_lineage}.mabs.csv")
    run:
        with open(output[0], "w") as f_out:
            writer = DictWriter(
                f_out,
                fieldnames=["timepoint", "chain", "sequence_id", "sequence"],
                lineterminator="\n")
            writer.writeheader()
            for attrs in ANTIBODY_ISOLATES.values():
                if attrs["Lineage"] == wildcards.antibody_lineage:
                    for chain in ["heavy", "light"]:
                        seq = attrs[chain.capitalize() + "Seq"]
                        if seq:
                            writer.writerow({
                                "timepoint": attrs["Timepoint"],
                                "chain": chain,
                                "sequence_id": attrs["Isolate"] + f"_{chain}",
                                "sequence": seq})
