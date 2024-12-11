"""
An alternate demux/trim/merge workflow for duplicated barcodes across loci.
"""

def input_for_merge_split_by_locus_final(w):
    samples_grouped = make_grouped_samples(w.run)
    for row in samples_grouped:
        if w.sample in row["SamplesMap"]:
            break
    else:
        raise ValueError
    if len(row["SamplesMap"]) > 1:
        return expand(
            "analysis/merge/{run}.grouped.split_by_locus/{sample}.{locus}.fastq.gz",
            run=w.run, sample=row["Sample"], locus=row["SamplesMap"][w.sample])
    # If there was actually only one sample for this particular case, we can
    # skip over the slow by-locus demux and just take the regular demux+merge
    # output
    return expand(
        "analysis/merge/{run}.grouped/{sample}.fastq.gz",
        run=w.run, sample=row["Sample"])

rule merge_split_by_locus_link:
    """Symlink merge file from locus-based workflow into main merge path"""
    output: "analysis/merge/{run}/{sample}.fastq.gz"
    input: "analysis/merge/{run}.grouped.final/{sample}.fastq.gz"
    shell: "ln -sr {input} {output}"

rule merge_split_by_locus_final:
    """Symlink real sample file to locus-demux'd file"""
    output: "analysis/merge/{run}.grouped.final/{sample}.fastq.gz"
    input: input_for_merge_split_by_locus_final
    shell: "ln -sr {input} {output}"

rule merge_split_by_locus:
    """Split files from a cross-locus demux using IgBLAST"""
    output: expand("analysis/merge/{{run}}.grouped.split_by_locus/{{sample}}.{locus}.fastq.gz", locus=["IGH", "IGK", "IGL", "other"])
    input: "analysis/merge/{run}.grouped/{sample}.fastq.gz"
    threads: 16
    params:
        output_pattern=lambda w: f"analysis/merge/{w.run}.grouped.split_by_locus/{w.sample}.{{locus}}.fastq.gz"
    shell: "split_by_locus.py -t {threads} {input} --output-pattern '{params.output_pattern}'"

def make_grouped_samples(runid):
    """Make placeholder sample info for samples with identical barcoding on one run"""
    out = []
    for attrs in SAMPLES.values():
        if attrs["Run"] == runid:
            out.append(attrs.copy())
    out_grouped = defaultdict(dict)
    for row in out:
        key = (row["BarcodeFwd"], row["BarcodeRev"])
        locus = {"lambda": "IGL", "kappa": "IGK"}.get(row["Type"], "IGH")
        if locus in out_grouped[key]:
            raise ValueError(
                f"Duplicate {locus} entry for Run/BCFwd/BCRev "
                f"{row['Run']}/{row['BarcodeFwd']}/{row['BarcodeRev']}")
        row["Locus"] = locus
        out_grouped[key][locus] = row
    out_flat = []
    for group in out_grouped.values():
        # make combo sample name sorted by locus, plus some extra entries for
        # keeping track of things
        locus_samples = sorted((key, row["Sample"]) for key, row in group.items())
        sample = "_".join(pair[1] for pair in locus_samples)
        samples = "/".join(pair[0] + ":" + pair[1] for pair in locus_samples)
        loci = "/".join(pair[0] for pair in locus_samples)
        rows = list(group.values())
        # make a single row out for each group using that sample name
        row_out = {key: rows[0][key] for key in ["Sample", "Run", "BarcodeFwd", "BarcodeRev"]}
        row_out["Sample"] = sample
        row_out["Loci"] = loci
        row_out["Samples"] = samples
        row_out["SamplesMap"] = {pair[1]: pair[0] for pair in locus_samples}
        out_flat.append(row_out)
    return out_flat

rule demux_locus_samples:
    """Make a placeholder samples.csv for cases that need demuxing by locus"""
    output: "analysis/demux/{run}.samples_grouped.csv"
    input: ancient("metadata/samples.csv")
    run:
        fields = ["Sample", "Run", "BarcodeFwd", "BarcodeRev", "Loci", "Samples"]
        sample_attrs = make_grouped_samples(wildcards.run)
        sample_attrs = [{k: v for k, v in attrs.items() if k in fields} for attrs in sample_attrs]
        with open(output[0], "w") as f_out:
            writer = csv.DictWriter(
                f_out, fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(sample_attrs)

def make_run_rules_locus_demux():
    """Generate dynamic per-sequencing-run rules where applicable (alternate workflow)

    This is equivalent to the per-run rules in by-run.smk but for the alternate
    workflow where we need to split files by antibody locus."""
    by_run = defaultdict(list)
    for samp_attrs in SAMPLES.values():
        if samp_attrs["Run"]:
            by_run[samp_attrs["Run"]].append(samp_attrs)
    for runid, samp_attrs_list in by_run.items():
        samp_names = [samp["Sample"] for samp in samp_attrs_list]
        bc_pairs = [(samp["BarcodeFwd"], samp["BarcodeRev"]) for samp in samp_attrs_list]
        if len(bc_pairs) > len(set(bc_pairs)):
            # special case where we have duplicated barcode pairs and need to
            # demux by locus too.
            samples_grouped = make_grouped_samples(runid)
            samp_names_grpd = [sample["Sample"] for sample in samples_grouped]
            rule:
                name: f"demux_{runid}"
                output:
                    outdir=directory(f"analysis/demux/{runid}.grouped"),
                    fqgz=expand("analysis/demux/{run}.grouped/{sample}.{rp}.fastq.gz", run=runid, sample=samp_names_grpd, rp=["R1", "R2"])
                input:
                    reads=f"analysis/reads/{runid}",
                    samples=f"analysis/demux/{runid}.samples_grouped.csv"
                log:
                    conda=f"analysis/demux/{runid}.grouped/conda_build.txt"
                conda: "envs/igseq.yaml"
                shell:
                    """
                        conda list --explicit > {log.conda}
                        igseq demux --samples {input.samples} --details {output.outdir}/details.csv.gz {input.reads} --outdir {output.outdir}
                    """
            # These are just helpers to group outputs from other rules by run
            # ID, but with some extra steps for this case
            rule:
                name: f"merge_{runid}"
                input: expand("analysis/merge/{run}/{sample}.fastq.gz", run=runid, sample=samp_names)
            rule:
                name: f"merge_{runid}_split_final"
                input: expand("analysis/merge/{run}.grouped.final/{sample}.fastq.gz", run=runid, sample=samp_names)
            rule:
                name: f"merge_{runid}_split"
                input: expand("analysis/merge/{run}.grouped.split_by_locus/{sample}.{locus}.fastq.gz", run=runid, sample=samp_names_grpd, locus=["IGH", "IGK", "IGL", "other"])
            rule:
                name: f"merge_{runid}_grouped"
                input: expand("analysis/merge/{run}.grouped/{sample}.fastq.gz", run=runid, sample=samp_names_grpd)
            rule:
                name: f"trim_{runid}"
                input: expand("analysis/trim/{run}.grouped/{sample}.{rp}.fastq.gz", run=runid, sample=samp_names_grpd, rp=["R1", "R2"])

make_run_rules_locus_demux()
