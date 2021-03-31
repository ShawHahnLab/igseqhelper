"""
Report generation and helper rules.

These prepare summary spreadsheets and plots from various steps in the
analysis.  It's a bit of a mess but the idea is everything funnels toward
analysis/report.pdf, with as much as possible handled in a modular fashion here
before going in the Rmarkdown for the report.
"""

from igseq.reporting.demux import make_barcode_summary
from igseq.reporting.trim import make_qualtrim_csv
from igseq.reporting.sonar import (
    get_rearr_centroids_by_raw_reads,
    rarefy_sonar_clusters,
    sonar_island_summary,
    sonar_island_stats)
from igseq.reporting.counts import (
    counts_sample_summary,
    counts_run_summary,
    counts_specimen_summary,
    counts_assembly_summary,
    counts_presto_qual_summary,
    counts_sonar_module1_summary)
from igseq.reporting.igdiscover import (
    convert_combined_alignment,
    combine_aligned_segments,
    align_next_segment,
    gather_antibodies)
from igseq.data import amplicon_files, transpose_sample_md

SAMPLE_MD_IGM = transpose_sample_md(SAMPLES, "IgM+")

# Only for samples where we have the run info
TARGET_QUALTRIM_GRID = expand(
    outputs_per_run(
        "analysis/reporting/by-run/{run}/qualtrim.{sample}.{{rp}}.csv",
        {key: SAMPLES[key] for key in SAMPLES if SAMPLES[key]["Run"]}),
    rp=["R1", "R2", "I1"])

TARGET_BARCODE_SUMMARY = expand(
    "analysis/reporting/by-run/{run}/barcode_summary.csv", run=RUNS.keys())

TARGET_IGDISCOVER_CLUSTERPLOTS = expand(
    "analysis/reporting/igdiscover/{chain}.{chain_type}/{specimen}/clusterplots.png",
    zip, chain=SAMPLE_MD_IGM["chains"], chain_type=SAMPLE_MD_IGM["chaintypes"], specimen=SAMPLE_MD_IGM["specimens"])

TARGET_SONAR_RAREFACTION = expand(
    "analysis/reporting/by-specimen/{specimen}.{antibody_lineage}.{chain}.{chain_type}/sonar_clusters_rarefaction.csv",
    zip, **igseq.sonar.setup_sonar_combos(
        transpose_sample_md({key: SAMPLES[key] for key in SAMPLES if SAMPLES[key]["Run"]}, "IgG+"),
        ANTIBODY_LINEAGES))

TARGET_SONAR_ISLAND_STATS = expand(
    "analysis/reporting/by-lineage/{antibody_lineage}/{specimen}.{chain}.{chain_type}/island_stats.csv",
    zip, **igseq.sonar.setup_sonar_combos(
        transpose_sample_md({key: SAMPLES[key] for key in SAMPLES if SAMPLES[key]["Run"]}, "IgG+"),
        ANTIBODY_LINEAGES))

rule all_qualtrim_grid:
    input: TARGET_QUALTRIM_GRID

rule all_barcode_summary:
    input: TARGET_BARCODE_SUMMARY

rule all_igdiscover_clusterplots:
    input: TARGET_IGDISCOVER_CLUSTERPLOTS

rule all_sonar_rarefaction:
    input: TARGET_SONAR_RAREFACTION

rule all_sonar_island_stats:
    input: TARGET_SONAR_ISLAND_STATS

TARGET_REPORT_ALL = [Path("analysis/report.pdf").resolve()]

TARGET_REPORT_INPUTS = expand(
    "analysis/reporting/counts/counts_{thing}.csv",
    thing=["by_sample", "by_run", "amplicon_summary",
           "assembly_summary", "presto_qual_summary", "sonar_module1_summary"]) + \
           TARGET_QUALTRIM_GRID + \
           TARGET_BARCODE_SUMMARY + \
           TARGET_IGDISCOVER_CLUSTERPLOTS + \
           TARGET_SONAR_RAREFACTION

TARGET_REPORT_COUNTS = expand(
    outputs_per_run("analysis/counts/demux/{run}/{{chunk}}/{sample}.{{rp}}.fastq.gz.counts", {key: SAMPLES[key] for key in SAMPLES if SAMPLES[key]["Run"]}),
    chunk=CHUNKS,
    rp=["R1", "R2", "I1"]) + expand(
        "analysis/counts/demux/{run}/{chunk}/unassigned.{rp}.fastq.gz.counts",
        run=set([entry["Run"] for entry in SAMPLES.values() if entry["Run"]]),
        chunk=CHUNKS,
        rp=["R1", "R2", "I1"])

TARGET_AMPLICON_COUNTS = amplicon_files(
    "analysis/counts/presto/data/{chain}.{chain_type}/{specimen}.R1.fastq.counts", SAMPLES)

TARGET_ASSEMBLY_COUNTS = amplicon_files(
    "analysis/counts/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq.counts", SAMPLES)

TARGET_PRESTO_QUAL_COUNTS = amplicon_files(
    "analysis/counts/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq.counts", SAMPLES)

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
    output: "analysis/reporting/by-run/{run}/qualtrim.{sample}.{rp}.csv"
    input: expand("analysis/demux/{{run}}/{chunk}/{{sample}}.{{rp}}.fastq.gz", chunk=CHUNKS)
    run: make_qualtrim_csv(input, output[0])

rule barcode_summary:
    """Make a CSV table summarizing barcodes identified for one run."""
    output: "analysis/reporting/by-run/{run}/barcode_summary.csv"
    input: expand("analysis/demux/{{run}}/{chunk}.barcodes.csv", chunk=CHUNKS)
    run:
        samples = {k: v for k, v in SAMPLES.items() if v["Run"] == wildcards.run}
        make_barcode_summary(output[0], input, SEQUENCES, samples)

def input_sonar_clusters_by_read(w):
    targets = {}
    subject = SPECIMENS[w.specimen]["Subject"]
    targets["rearr"] = "analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/{specimen}/output/tables/{specimen}_rearrangements.tsv".format(
        subject=subject, chain=w.chain, chain_type=w.chain_type, antibody_lineage=w.antibody_lineage, specimen=w.specimen)
    for samp_name, samp_attrs in SAMPLES.items():
        if samp_attrs["Specimen"] == w.specimen and \
            samp_attrs["Chain"] == w.chain and \
            samp_attrs["Type"] == w.chain_type:
            if not samp_name in targets:
                targets[samp_name] = []
            targets[samp_name] += expand("analysis/demux/{run}/{chunk}/{sample}.I1.fastq.gz",
                run=samp_attrs["Run"],
                chunk=CHUNKS,
                sample=samp_name)
    return targets

rule sonar_clusters_by_read:
    """Make a CSV pairing raw sequence read IDs with the SONAR cluster they map to."""
    output: "analysis/reporting/by-specimen/{specimen}.{antibody_lineage}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    input: unpack(input_sonar_clusters_by_read)
    run:
        fp_rearr = input.rearr
        fps_fqgz = []
        for key, val in input.items():
            if key is not "rearr":
                fps_fqgz.extend(val)
        get_rearr_centroids_by_raw_reads(fp_rearr, fps_fqgz, output[0])

rule rarefy_sonar_clusters:
    output: "analysis/reporting/by-specimen/{specimen}.{antibody_lineage}.{chain}.{chain_type}/sonar_clusters_rarefaction.csv"
    input: "analysis/reporting/by-specimen/{specimen}.{antibody_lineage}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    run: rarefy_sonar_clusters(input[0], output[0])

# Summaries of SONAR's ID/DIV information for the selected "islands" of sequences

def sonar_island_stats_input(w):
    subject = SPECIMENS[w.specimen]["Subject"]
    prefix = "analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/{specimen}/output/"
    targets = {
        "island": expand(
            prefix + "tables/islandSeqs.txt",
            subject=subject, chain=w.chain, chain_type=w.chain_type, antibody_lineage=w.antibody_lineage, specimen=w.specimen)[0],
        "iddiv": expand(
            prefix + "tables/{specimen}_goodVJ_unique_id-div.tab",
            subject=subject, chain=w.chain, chain_type=w.chain_type, antibody_lineage=w.antibody_lineage, specimen=w.specimen)[0],
        "fasta": expand(
            prefix + "sequences/nucleotide/{specimen}_islandSeqs.fa",
            subject=subject, chain=w.chain, chain_type=w.chain_type, antibody_lineage=w.antibody_lineage, specimen=w.specimen)[0]}
    return targets

rule sonar_island_stats:
    """Condense the full ID/DIV stats to just those for one island and sumamrize across antibodies."""
    output: "analysis/reporting/by-lineage/{antibody_lineage}/{specimen}.{chain}.{chain_type}/island_stats.csv"
    input: unpack(sonar_island_stats_input)
    run:
        timepoint = [attrs["Timepoint"] for spec, attrs in SPECIMENS.items() if spec == wildcards.specimen][0]
        sonar_island_stats(
            output[0], input.island, input.iddiv, input.fasta,
            extras={
                "specimen": wildcards.specimen,
                "timepoint": timepoint})

def sonar_island_summary_input(w):
    keep = lambda spec: spec["Subject"] == ANTIBODY_LINEAGES[w.antibody_lineage]["Subject"] and "IgG" in spec["CellType"]
    specimens = [specname for specname, specattrs in SPECIMENS.items() if keep(specattrs)]
    return expand(
        "analysis/reporting/by-lineage/{antibody_lineage}/{specimen}.{chain}.{chain_type}/island_stats.csv",
        antibody_lineage=w.antibody_lineage,
        specimen=specimens,
        chain=w.chain,
        chain_type=w.chain_type)

rule sonar_island_summary:
    """Further condense ID/DIV stats to one file per lineage."""
    output: "analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/island_stats_summary.csv"
    input: sonar_island_summary_input
    run: sonar_island_summary(output[0], input)

# Alignments of antibody isolate sequences with the germline V(D)J alleles

rule lineage_allele_alignments_table:
    """Convert combined alignment into CSV of relevant attributes per sequence."""
    output: csv="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/VDJ.aligned.csv"
    input: fasta="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/VDJ.aligned.fasta"
    run:
        convert_combined_alignment(
            input.fasta, output.csv, ANTIBODY_LINEAGES, ANTIBODY_ISOLATES, wildcards)

rule lineage_align_vdj:
    """Combine antibody/germline V(D)J alignments into one combined alignment."""
    output: fasta="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/VDJ.aligned.fasta"
    input:
        target="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies.aln.fa",
        with_v="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/V.aligned.fasta",
        with_d="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/D.aligned.fasta",
        with_j="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/J.aligned.fasta"
    run:
        if wildcards.chain == "heavy":
            with_d = input.with_d
        else:
            with_d = None
        combine_aligned_segments(\
            input.target, input.with_v, with_d, input.with_j, output.fasta)

rule lineage_align_j:
    """Align J alleles to target sequences based on V/J alignment."""
    output: fasta="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/J.aligned.fasta"
    input:
        target="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies.aln.fa",
        with_d="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/D.aligned.fasta",
        j="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/germline.J.fasta"
    run: align_next_segment(input.target, input.with_d, input.j, output.fasta)

rule lineage_align_d:
    """Align D alleles to target sequences based on V alignment."""
    output: fasta="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/D.aligned.fasta"
    input:
        target="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies.aln.fa",
        with_v="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/V.aligned.fasta",
        d="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/germline.D.fasta"
    run:
        if wildcards.chain == "heavy":
            align_next_segment(
                input.target, input.with_v, input.d, output.fasta)
        else:
            shell("cp {input.with_v} {output.fasta}")

rule lineage_align_v:
    """Align V alleles to target sequences."""
    output: fasta="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies/V.aligned.fasta"
    input:
        target="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies.aln.fa",
        v="analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/germline.V.fasta"
    shell:
        """
            if [[ -s {input.target} ]]; then
                clustalw -align -profile1={input.target} -profile2={input.v} -sequences -output=fasta -outfile={output}
            else
                touch {output}
            fi
        """

rule lineage_align_abs:
    """Align target sequences to one another before bringing discovered alleles into the mix."""
    output: "analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies.aln.fa"
    input: "analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies.fasta"
    shell:
        """
            if [[ -s {input} ]]; then
                clustalw -align -infile={input} -output=fasta -outfile={output}
            else
                touch {output}
            fi
        """

rule lineage_gather_antibody_sequences:
    """Gather mature antibody sequences to start off the alignment against."""
    output: "analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/antibodies.fasta"
    run:
        gather_antibodies(
            wildcards.antibody_lineage, wildcards.chain,
            ANTIBODY_ISOLATES, ANTIBODY_LINEAGES, output[0])

rule lineage_gather_germline:
    """Get the germline sequences for one segment for a given lineage's subject."""
    # We can just use the existing logic from the SONAR rules to match lineage
    # to subject to relevant germline specimens and results.
    output: "analysis/reporting/by-lineage/{antibody_lineage}/{chain}.{chain_type}/germline.{segment}.fasta"
    input: lambda w: expand("analysis/sonar/{subject}/{{chain}}.{{chain_type}}/germline.{{segment}}.fasta", subject = ANTIBODY_LINEAGES[w.antibody_lineage]["Subject"])
    shell: "cp {input} {output}"

# Sample-based

# All the different samples and runs for samples that list their runs,
# AND the two unassigned count files for each define run,
# AND expand it across all the chunked files.
# Whew.
INPUT_COUNTS_SAMPLE_SUMMARY = expand(
    expand(
        "analysis/counts/demux/{run}/{{chunk}}/{item}.I1.fastq.gz.counts",
        zip,
        run=[s["Run"] for s in SAMPLES.values() if s["Run"]],
        item=[s["Sample"] for s in SAMPLES.values() if s["Run"]]) +
    expand(
        "analysis/counts/demux/{run}/{{chunk}}/{item}.I1.fastq.gz.counts",
        run={s["Run"] for s in SAMPLES.values() if s["Run"]},
        item=["unassigned", "unassigned.phix"]),
    chunk=CHUNKS)

rule counts_sample_summary:
    """A per-sample summary of raw read counts."""
    output: "analysis/reporting/counts/counts_by_sample.csv"
    input: INPUT_COUNTS_SAMPLE_SUMMARY
    run: counts_sample_summary(input, output[0], SAMPLES)

rule counts_run_summary:
    """A per-run summary of raw read counts."""
    output: "analysis/reporting/counts/counts_by_run.csv"
    input: "analysis/reporting/counts/counts_by_sample.csv"
    run: counts_run_summary(input[0], output[0])

# Specimen+chain -based

rule counts_presto_amplicon_summary:
    """A per-amplicon summary of read counts."""
    output: "analysis/reporting/counts/counts_amplicon_summary.csv"
    input: TARGET_AMPLICON_COUNTS
    run: counts_specimen_summary(input, output[0], SPECIMENS)

rule counts_presto_assembly_summary:
    """A per-specimen summary of paired read counts."""
    output: "analysis/reporting/counts/counts_assembly_summary.csv"
    input: TARGET_ASSEMBLY_COUNTS
    run: counts_assembly_summary(input, output[0], SPECIMENS)

rule counts_presto_qual_summary:
    """A per-specimen summary of paired read counts."""
    output: "analysis/reporting/counts/counts_presto_qual_summary.csv"
    input: TARGET_PRESTO_QUAL_COUNTS
    run: counts_presto_qual_summary(input, output[0], SPECIMENS)

rule counts_sonar_module1_summary:
    """A per-specimen summary of clustered and categorized read counts."""
    output: "analysis/reporting/counts/counts_sonar_module1_summary.csv"
    input: expand("analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/{specimen}/output/tables/{specimen}_rearrangements.tsv", zip, **igseq.sonar.setup_sonar_combos(SAMPLE_MD_IGG, ANTIBODY_LINEAGES))
    run:
        counts_sonar_module1_summary(
            input,
            output[0],
            igseq.sonar.setup_sonar_combos(SAMPLE_MD_IGG, ANTIBODY_LINEAGES))

# TODO next: also tally after pRESTO's QC and primer checking.  Maybe do a
# summary table in the report of counts following assembly, qc, and primer
# checking showing attrition across steps.

rule igdiscover_clusterplot_grid:
    """A summary grid PNG of all the various clusterplot PNGs for a single IgDiscover dir.

    This will make an 8xN grid of tiled thumbnails for quick viewing.
    """
    output: "analysis/reporting/igdiscover/{chain}.{chain_type}/{specimen}/clusterplots.png"
    input: "analysis/igdiscover/{chain}.{chain_type}/{specimen}/stats/stats.json"
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
