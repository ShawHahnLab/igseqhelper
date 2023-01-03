"""
Rules for dealing with data aggregated to the specimen-level, grouped by chain.

This includes aggregation of samples into per-specimen files, pRESTO QC, and
IgDiscover (for IgM+/IgD+ specimens).
"""

def reads_by_specimen_input(wildcards):
    """Get trimmed sequencer samples matching specimen+antibody type.

    We may have prepped and sequenced the same specimen multiple times so this
    will give a list of one or more files.
    """
    pattern = "analysis/reads-by-sample/{sample}.{rp}.fastq.gz"
    filepaths = []
    for samp_name, samp_attrs in SAMPLES.items():
        if \
            samp_attrs["Specimen"] == wildcards.specimen and \
            samp_attrs["Chain"] == wildcards.chain and \
            samp_attrs["Type"] == wildcards.chain_type:
            filepaths.append(pattern.format(
                sample=samp_name,
                rp=wildcards.rp))
    if not filepaths:
        raise ValueError("No matching samples found for %s & %s & %s" %
            (wildcards.specimen, wildcards.chain, wildcards.chain_type))
    return filepaths

TARGET_READS_BY_SPECIMEN = expand(
        "analysis/reads-by-specimen/{chain}.{chain_type}/{specimen}.{rp}.fastq.gz",
        zip,
        chain=[samp["Chain"] for samp in SAMPLES.values()],
        chain_type=[samp["Type"] for samp in SAMPLES.values()],
        specimen=[samp["Specimen"] for samp in SAMPLES.values()],
        rp=["R1", "R2"])

rule all_reads_by_specimen:
    input: TARGET_READS_BY_SPECIMEN

rule reads_by_specimen:
    """Aggregate trimmed sequencer samples into per-specimen FASTQ files."""
    output: "analysis/reads-by-specimen/{chain}.{chain_type}/{specimen}.{rp}.fastq.gz"
    input: reads_by_specimen_input
    run:
        if len(input) == 1:
            path_orig = Path(input[0]).resolve()
            shell("ln -s {path_orig} {output}")
        else:
            shell("cat {input} > {output}")

include: "presto.smk"
