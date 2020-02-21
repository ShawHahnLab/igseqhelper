"""
Various QC checks.
"""

from igseq.qc import make_qualtrim_csv, get_rearr_centroids_by_raw_reads, rarefy_sonar_clusters

TARGET_QUALTRIM_GRID = expand(
    outputs_per_run("qc/{run}/qualtrim.{sample}.{{rp}}.csv", SAMPLES),
    rp=["R1", "R2", "I1"])

rule all_qualtrim_grid:
    input: TARGET_QUALTRIM_GRID

rule qualtrim_grid:
    """Make a CSV table summarizing cutadapt trim cutoffs vs output length."""
    output: "qc/{run}/qualtrim.{sample}.{rp}.csv"
    input: "demux/{run}/{sample}.{rp}.fastq.gz"
    run: make_qualtrim_csv(input[0], output[0])

def input_sonar_clusters_by_read(w):
    targets = {}
    subject = SPECIMENS[w.specimen]["Subject"]
    targets["rearr"] = "sonar-analysis/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv".format(
        subject=subject, chain=w.chain, chain_type=w.chain_type, specimen=w.specimen)
    for samp_name, samp_attrs in SAMPLES.items():
        if samp_attrs["Specimen"] == w.specimen and \
            samp_attrs["Chain"] == w.chain and \
            samp_attrs["Type"] == w.chain_type:
            targets[samp_name] = "demux/{run}/{sample}.I1.fastq.gz".format(
                run=samp_attrs["Run"],
                sample=samp_name)
    return targets

rule sonar_clusters_by_read:
    """Make a CSV pairing raw sequence read IDs with the SONAR cluster they map to."""
    output: "qc/{specimen}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    input: unpack(input_sonar_clusters_by_read)
    run:
        fp_rearr = input.rearr
        fps_fqgz = [val for key, val in input.items() if key is not "rearr"]
        get_rearr_centroids_by_raw_reads(fp_rearr, fps_fqgz, output[0])

rule rarefy_sonar_clusters:
    output: "qc/{specimen}.{chain}.{chain_type}/sonar_clusters_rarefaction.csv"
    input: "qc/{specimen}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    run: rarefy_sonar_clusters(input[0], output[0])
