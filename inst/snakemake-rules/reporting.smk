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

TARGET_IGDISCOVER_ALLELE_ALIGNMENTS = expand(
    "reporting/{specimen}.{chain}.{chain_type}/germline.VDJ.aligned.csv",
    zip, chain=IGM_CHAINS, chain_type=IGM_CHAINTYPES, specimen=IGM_SPECIMENS)

TARGET_SONAR_RAREFACTION = expand(
    "reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_rarefaction.csv",
    zip, chain=IGG_CHAINS, chain_type=IGG_CHAINTYPES, specimen=IGG_SPECIMENS)

rule all_qualtrim_grid:
    input: TARGET_QUALTRIM_GRID

rule all_igdiscover_clusterplots:
    input: TARGET_IGDISCOVER_CLUSTERPLOTS

rule all_igdiscover_allele_alignments:
    input: TARGET_IGDISCOVER_ALLELE_ALIGNMENTS

rule all_sonar_rarefaction:
    input: TARGET_SONAR_RAREFACTION

TARGET_REPORT_ALL = [Path("report.pdf").resolve()]

TARGET_REPORT_INPUTS = expand(
    "reporting/{thing}.csv",
    thing=["counts_by_sample", "counts_by_run", "counts_amplicon_summary",
           "counts_assembly_summary", "counts_presto_qual_summary"]) + \
           TARGET_QUALTRIM_GRID + \
           TARGET_IGDISCOVER_CLUSTERPLOTS + \
           TARGET_IGDISCOVER_ALLELE_ALIGNMENTS

TARGET_REPORT_COUNTS = expand(
    outputs_per_run("counts/demux/{run}/{sample}.{{rp}}.fastq.gz.counts", SAMPLES),
    rp=["R1", "R2", "I1"]) + expand(
        "counts/demux/{run}/unassigned.{rp}.fastq.gz.counts",
        run=set([entry["Run"] for entry in SAMPLES.values()]),
        rp=["R1", "R2", "I1"])

TARGET_AMPLICON_COUNTS = amplicon_files(
    "counts/presto/data/{chain}.{chain_type}/{specimen}.R1.fastq.counts", SAMPLES)

TARGET_ASSEMBLY_COUNTS = amplicon_files(
    "counts/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq.counts", SAMPLES)

TARGET_PRESTO_QUAL_COUNTS = amplicon_files(
    "counts/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq.counts", SAMPLES)

rule render_report:
    output: TARGET_REPORT_ALL
    input: TARGET_REPORT_INPUTS
    params: report=RINST/"report.Rmd"
    shell:
        """
            Rscript --vanilla -e 'rmarkdown::render("{params.report}", output_file="{output}", quiet=TRUE)'
        """

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

### IgDiscover allele handling

rule igdiscover_allele_alignments_table:
    """Convert combined alignment into CSV of relevant attributes per sequence."""
    output: csv="reporting/{specimen}.{chain}.{chain_type}/germline.VDJ.aligned.csv"
    input: fasta="reporting/{specimen}.{chain}.{chain_type}/germline.VDJ.aligned.fasta"
    run:
        igseq.reporting.convert_combined_alignment(
            input.fasta, output.csv, SPECIMENS, ANTIBODY_ISOLATES, wildcards)

rule igdiscover_allele_alignments_align_vdj:
    """Combine V(D)J alignments into one combined alignment."""
    output: fasta="reporting/{specimen}.{chain}.{chain_type}/germline.VDJ.aligned.fasta"
    input:
        antibodies="reporting/{specimen}.{chain}.{chain_type}/antibodies.fasta",
        with_v="reporting/{specimen}.{chain}.{chain_type}/germline.V.aligned.fasta",
        with_d="reporting/{specimen}.{chain}.{chain_type}/germline.D.aligned.fasta",
        with_j="reporting/{specimen}.{chain}.{chain_type}/germline.J.aligned.fasta"
    run:
        if wildcards.chain == "heavy":
            with_d = input.with_d
        else:
            with_d = None
        igseq.reporting.combine_aligned_segments(\
            input.antibodies, input.with_v, with_d, input.with_j, output.fasta)

rule igdiscover_allele_alignments_align_j:
    """Align discovered J alleles to antibody lineages based on V/J alignment."""
    output: fasta="reporting/{specimen}.{chain}.{chain_type}/germline.J.aligned.fasta"
    input:
        antibodies="reporting/{specimen}.{chain}.{chain_type}/antibodies.fasta",
        with_d="reporting/{specimen}.{chain}.{chain_type}/germline.D.aligned.fasta",
        j="igdiscover/{chain}.{chain_type}/{specimen}/final/database/J.fasta"
    run: igseq.reporting.align_next_segment(input.antibodies, input.with_d, input.j, output.fasta)

rule igdiscover_allele_alignments_align_d:
    """Align discovered D alleles to antibody lineages based on V alignment."""
    output: fasta="reporting/{specimen}.{chain}.{chain_type}/germline.D.aligned.fasta"
    input:
        antibodies="reporting/{specimen}.{chain}.{chain_type}/antibodies.fasta",
        with_v="reporting/{specimen}.{chain}.{chain_type}/germline.V.aligned.fasta",
        d="igdiscover/{chain}.{chain_type}/{specimen}/final/database/D.fasta"
    run:
        if wildcards.chain == "heavy":
            igseq.reporting.align_next_segment(
                input.antibodies, input.with_v, input.d, output.fasta)
        else:
            shell("cp {input.with_v} {output.fasta}")

rule igdiscover_allele_alignments_align_v:
    """Align discovered V alleles to antibody lineages."""
    output: fasta="reporting/{specimen}.{chain}.{chain_type}/germline.V.aligned.fasta"
    input:
        antibodies="reporting/{specimen}.{chain}.{chain_type}/antibodies.fasta",
        v="igdiscover/{chain}.{chain_type}/{specimen}/final/database/V.fasta",
    shell:
        """
            if [[ -s {input.antibodies} ]]; then
                clustalw -align -profile1={input.antibodies} -profile2={input.v} -sequences -output=fasta -outfile={output}
            else
                touch {output}
            fi
        """

rule gather_antibody_sequences:
    """Gather mature antibody sequences to start off the alignment against."""
    output: "reporting/{specimen}.{chain}.{chain_type}/antibodies.fasta"
    run:
        igseq.reporting.gather_antibodies(
            SPECIMENS[wildcards.specimen]["Subject"],
            wildcards.chain, ANTIBODY_ISOLATES, output[0])
