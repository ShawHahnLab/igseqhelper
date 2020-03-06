### Reporting

import igseq.reporting
from igseq.data import amplicon_files

IGM_CHAINS = [entry["Chain"] for entry in SAMPLES.values() if "IgM+" in entry["SpecimenAttrs"]["CellType"]]
IGM_CHAINTYPES = [entry["Type"] for entry in SAMPLES.values() if "IgM+" in entry["SpecimenAttrs"]["CellType"]]
IGM_SUBJECTS = [entry["SpecimenAttrs"]["Subject"] for entry in SAMPLES.values() if "IgM+" in entry["SpecimenAttrs"]["CellType"]]
IGM_SPECIMENS = [entry["Specimen"] for entry in SAMPLES.values() if "IgM+" in entry["SpecimenAttrs"]["CellType"]]

TARGET_QUALTRIM_GRID = expand(
    outputs_per_run("reporting/{run}/qualtrim.{sample}.{{rp}}.csv", SAMPLES),
    rp=["R1", "R2", "I1"])

TARGET_IGDISCOVER_CLUSTERPLOTS = expand(
    "reporting/igdiscover/{chain}.{chain_type}/{specimen}/clusterplots.png",
    zip, chain=IGM_CHAINS, chain_type=IGM_CHAINTYPES, specimen=IGM_SPECIMENS)

rule all_qualtrim_grid:
    input: TARGET_QUALTRIM_GRID

rule all_igdiscover_clusterplots:
    input: TARGET_IGDISCOVER_CLUSTERPLOTS

TARGET_REPORT_INPUTS = expand(
    "reporting/{thing}.csv",
    thing=["counts_by_sample", "counts_by_run", "counts_amplicon_summary",
           "counts_assembly_summary", "counts_presto_qual_summary"]) + \
           TARGET_QUALTRIM_GRID + \
           TARGET_IGDISCOVER_CLUSTERPLOTS

TARGET_REPORT_ALL = ["report.pdf"]

TARGET_REPORT_COUNTS = expand(
    outputs_per_run("counts/demux/{run}/{sample}.{{rp}}.fastq.gz.counts", SAMPLES),
    rp=["R1", "R2", "I1"]) + expand(
    "counts/demux/{run}/unassigned.{rp}.fastq.gz.counts", run=RUNS.keys(), rp=["R1", "R2", "I1"])

TARGET_AMPLICON_COUNTS = amplicon_files(
    "counts/presto/data/{chain}.{chain_type}/{specimen}.R1.fastq.counts", SAMPLES)

TARGET_ASSEMBLY_COUNTS = amplicon_files(
    "counts/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq.counts", SAMPLES)

TARGET_PRESTO_QUAL_COUNTS = amplicon_files(
    "counts/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq.counts", SAMPLES)

rule render_report:
    output: TARGET_REPORT_ALL
    input: TARGET_REPORT_INPUTS
    script: "igseq/inst/report.Rmd"

rule qualtrim_grid:
    """Make a CSV table summarizing cutadapt trim cutoffs vs output length."""
    output: "reporting/{run}/qualtrim.{sample}.{rp}.csv"
    input: "demux/{run}/{sample}.{rp}.fastq.gz"
    run: igseq.reporting.make_qualtrim_csv(input[0], output[0])

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
    output: "reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    input: unpack(input_sonar_clusters_by_read)
    run:
        fp_rearr = input.rearr
        fps_fqgz = [val for key, val in input.items() if key is not "rearr"]
        igseq.reporting.get_rearr_centroids_by_raw_reads(fp_rearr, fps_fqgz, output[0])

rule rarefy_sonar_clusters:
    output: "reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_rarefaction.csv"
    input: "reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    run: igseq.reporting.rarefy_sonar_clusters(input[0], output[0])

# Sample-based

rule counts_table:
    """Just a big ol' list of file paths and read counts."""
    output: "reporting/counts.csv"
    input: TARGET_REPORT_COUNTS
    shell:
        """
            echo "Filename,NumSequences" > {output}
            paste -d , <(echo "{input}" | tr ' ' '\\n') <(cat {input}) >> {output}
        """

rule counts_sample_summary:
    """A per-sample summary of raw read counts."""
    output: "reporting/counts_by_sample.csv"
    input: "reporting/counts.csv"
    run: igseq.reporting.counts_sample_summary(input[0], output[0], SAMPLES)

rule counts_run_summary:
    """A per-run summary of raw read counts."""
    output: "reporting/counts_by_run.csv"
    input: "reporting/counts_by_sample.csv"
    run: igseq.reporting.counts_run_summary(input[0], output[0])

# Specimen+chain -based

rule counts_presto_amplicon_summary:
    """A per-amplicon summary of read counts."""
    output: "reporting/counts_amplicon_summary.csv"
    input: TARGET_AMPLICON_COUNTS
    run: igseq.reporting.counts_specimen_summary(input, output[0], SPECIMENS)

rule counts_presto_assembly_summary:
    """A per-specimen summary of paired read counts."""
    output: "reporting/counts_assembly_summary.csv"
    input: TARGET_ASSEMBLY_COUNTS
    run: igseq.reporting.counts_assembly_summary(input, output[0], SPECIMENS)

rule counts_presto_qual_summary:
    """A per-specimen summary of paired read counts."""
    output: "reporting/counts_presto_qual_summary.csv"
    input: TARGET_PRESTO_QUAL_COUNTS
    run: igseq.reporting.counts_presto_qual_summary(input, output[0], SPECIMENS)

# TODO next: also tally after pRESTO's QC and primer checking.  Maybe do a
# summary table in the report of counts following assembly, qc, and primer
# checking showing attrition across steps.

rule igdiscover_clusterplot_grid:
    """A summary grid PNG of all the various clusterplot PNGs for a single IgDiscover dir.

    This will make an 8xN grid of tiled thumbnails for quick viewing.
    """
    output: "reporting/igdiscover/{chain}.{chain_type}/{specimen}/clusterplots.png"
    input: "igdiscover/{chain}.{chain_type}/{specimen}/stats/stats.json"
    params:
        width=100,
        height=100,
        cols=8
    shell:
        """
            plotdir=$(dirname {input})/../final/clusterplots
            num=$(ls $plotdir/*.png | wc -l)
            num=$((num / {params.cols} + 1))
            # My install always complains about fonts and returns with exit
            # code 1.  Not ideal but I'm just ignoring it for now.
            montage $plotdir/*.png -geometry {params.width}x{params.height}+0+0 -tile {params.cols}x$num {output} || exitval=$?
        """
