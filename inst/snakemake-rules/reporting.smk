"""
Report generation and helper rules.

These prepare summary spreadsheets and plots from various steps in the
analysis.  It's a bit of a mess but the idea is everything funnels toward
analysis/report.pdf, with as much as possible handled in a modular fashion here
before going in the Rmarkdown for the report.
"""

from igseq.reporting.trim import make_qualtrim_csv
from igseq.reporting.sonar import (
    get_rearr_centroids_by_raw_reads,
    rarefy_sonar_clusters)
from igseq.reporting.counts import (
    counts_sample_summary,
    counts_run_summary,
    counts_specimen_summary,
    counts_assembly_summary,
    counts_presto_qual_summary)
from igseq.reporting.igdiscover import (
    convert_combined_alignment,
    combine_aligned_segments,
    align_next_segment,
    gather_antibodies)
from igseq.data import amplicon_files, transpose_sample_md

SAMPLE_MD_IGM = transpose_sample_md(SAMPLES, "IgM+")

TARGET_QUALTRIM_GRID = expand(
    outputs_per_run("analysis/reporting/{run}/qualtrim.{sample}.{{rp}}.csv", SAMPLES),
    rp=["R1", "R2", "I1"])

TARGET_IGDISCOVER_CLUSTERPLOTS = expand(
    "analysis/reporting/igdiscover/{chain}.{chain_type}/{specimen}/clusterplots.png",
    zip, chain=SAMPLE_MD_IGM["chains"], chain_type=SAMPLE_MD_IGM["chaintypes"], specimen=SAMPLE_MD_IGM["specimens"])

def _get_allele_alignment_targets():
    targets = []
    for sample, sample_attrs in SAMPLES.items():
        specimen = sample_attrs["Specimen"]
        if not specimen in SAMPLE_MD_IGM["specimens"]:
            continue
        chain = sample_attrs["Chain"]
        chain_type = sample_attrs["Type"]
        subject = sample_attrs["SpecimenAttrs"]["Subject"]
        for lineage, lineage_attrs in ANTIBODY_LINEAGES.items():
            if lineage_attrs["Subject"] == subject:
                targets.append(
                    ("analysis/reporting/{specimen}.{chain}.{chain_type}"
                    "/antibodies.{lineage}/VDJ.aligned.csv").format(
                        specimen=specimen,
                        chain=chain,
                        chain_type=chain_type,
                        lineage=lineage))
    return targets

TARGET_IGDISCOVER_ALLELE_ALIGNMENTS = _get_allele_alignment_targets()

TARGET_SONAR_RAREFACTION = expand(
    "analysis/reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_rarefaction.csv",
    zip, chain=SAMPLE_MD_IGG["chains"], chain_type=SAMPLE_MD_IGG["chaintypes"], specimen=SAMPLE_MD_IGG["specimens"])

rule all_qualtrim_grid:
    input: TARGET_QUALTRIM_GRID

rule all_igdiscover_clusterplots:
    input: TARGET_IGDISCOVER_CLUSTERPLOTS

rule all_igdiscover_allele_alignments:
    input: TARGET_IGDISCOVER_ALLELE_ALIGNMENTS

rule all_sonar_rarefaction:
    input: TARGET_SONAR_RAREFACTION

TARGET_REPORT_ALL = [Path("analysis/report.pdf").resolve()]

TARGET_REPORT_INPUTS = expand(
    "analysis/reporting/{thing}.csv",
    thing=["counts_by_sample", "counts_by_run", "counts_amplicon_summary",
           "counts_assembly_summary", "counts_presto_qual_summary"]) + \
           TARGET_QUALTRIM_GRID + \
           TARGET_IGDISCOVER_CLUSTERPLOTS + \
           TARGET_IGDISCOVER_ALLELE_ALIGNMENTS

TARGET_REPORT_COUNTS = expand(
    outputs_per_run("analysis/counts/demux/{run}/{{chunk}}/{sample}.{{rp}}.fastq.gz.counts", SAMPLES),
    chunk=CHUNKS,
    rp=["R1", "R2", "I1"]) + expand(
        "analysis/counts/demux/{run}/{chunk}/unassigned.{rp}.fastq.gz.counts",
        run=set([entry["Run"] for entry in SAMPLES.values()]),
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
    output: "analysis/reporting/{run}/qualtrim.{sample}.{rp}.csv"
    input: expand("analysis/demux/{{run}}/{chunk}/{{sample}}.{{rp}}.fastq.gz", chunk=CHUNKS)
    run: make_qualtrim_csv(input, output[0])

def input_sonar_clusters_by_read(w):
    targets = {}
    subject = SPECIMENS[w.specimen]["Subject"]
    targets["rearr"] = "analysis/sonar/{subject}/{chain}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv".format(
        subject=subject, chain=w.chain, chain_type=w.chain_type, specimen=w.specimen)
    for samp_name, samp_attrs in SAMPLES.items():
        if samp_attrs["Specimen"] == w.specimen and \
            samp_attrs["Chain"] == w.chain and \
            samp_attrs["Type"] == w.chain_type:
            targets[samp_name] = "analysis/demux/{run}/{sample}.I1.fastq.gz".format(
                run=samp_attrs["Run"],
                sample=samp_name)
    return targets

rule sonar_clusters_by_read:
    """Make a CSV pairing raw sequence read IDs with the SONAR cluster they map to."""
    output: "analysis/reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    input: unpack(input_sonar_clusters_by_read)
    run:
        fp_rearr = input.rearr
        fps_fqgz = [val for key, val in input.items() if key is not "rearr"]
        get_rearr_centroids_by_raw_reads(fp_rearr, fps_fqgz, output[0])

rule rarefy_sonar_clusters:
    output: "analysis/reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_rarefaction.csv"
    input: "analysis/reporting/{specimen}.{chain}.{chain_type}/sonar_clusters_by_read.csv"
    run: rarefy_sonar_clusters(input[0], output[0])

# Sample-based

rule counts_table:
    """Just a big ol' list of file paths and read counts."""
    output: "analysis/reporting/counts.csv"
    input: TARGET_REPORT_COUNTS
    run:
        import csv
        with open(output[0], "wt") as f_out:
            writer = csv.DictWriter(
                f_out, fieldnames=["Filename", "NumSequences"], lineterminator="\n")
            writer.writeheader()
            for countsfile in input:
                with open(countsfile) as f_in:
                    counts = f_in.read().strip()
                writer.writerow({"Filename": countsfile, "NumSequences": counts})

rule counts_sample_summary:
    """A per-sample summary of raw read counts."""
    output: "analysis/reporting/counts_by_sample.csv"
    input: "analysis/reporting/counts.csv"
    run: counts_sample_summary(input[0], output[0], SAMPLES)

rule counts_run_summary:
    """A per-run summary of raw read counts."""
    output: "analysis/reporting/counts_by_run.csv"
    input: "analysis/reporting/counts_by_sample.csv"
    run: counts_run_summary(input[0], output[0])

# Specimen+chain -based

rule counts_presto_amplicon_summary:
    """A per-amplicon summary of read counts."""
    output: "analysis/reporting/counts_amplicon_summary.csv"
    input: TARGET_AMPLICON_COUNTS
    run: counts_specimen_summary(input, output[0], SPECIMENS)

rule counts_presto_assembly_summary:
    """A per-specimen summary of paired read counts."""
    output: "analysis/reporting/counts_assembly_summary.csv"
    input: TARGET_ASSEMBLY_COUNTS
    run: counts_assembly_summary(input, output[0], SPECIMENS)

rule counts_presto_qual_summary:
    """A per-specimen summary of paired read counts."""
    output: "analysis/reporting/counts_presto_qual_summary.csv"
    input: TARGET_PRESTO_QUAL_COUNTS
    run: counts_presto_qual_summary(input, output[0], SPECIMENS)

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

### IgDiscover allele handling

rule igdiscover_allele_alignments_table:
    """Convert combined alignment into CSV of relevant attributes per sequence."""
    output: csv="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/VDJ.aligned.csv"
    input: fasta="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/VDJ.aligned.fasta"
    run:
        convert_combined_alignment(
            input.fasta, output.csv, SPECIMENS, ANTIBODY_ISOLATES, wildcards)

rule igdiscover_allele_alignments_align_vdj:
    """Combine V(D)J alignments into one combined alignment."""
    output: fasta="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/VDJ.aligned.fasta"
    input:
        target="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}.aln.fa",
        with_v="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/V.aligned.fasta",
        with_d="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/D.aligned.fasta",
        with_j="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/J.aligned.fasta"
    run:
        if wildcards.chain == "heavy":
            with_d = input.with_d
        else:
            with_d = None
        combine_aligned_segments(\
            input.target, input.with_v, with_d, input.with_j, output.fasta)

rule igdiscover_allele_alignments_align_j:
    """Align discovered J alleles to target sequences based on V/J alignment."""
    output: fasta="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/J.aligned.fasta"
    input:
        target="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}.aln.fa",
        with_d="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/D.aligned.fasta",
        j="analysis/igdiscover/{chain}.{chain_type}/{specimen}/final/database/J.fasta"
    run: align_next_segment(input.target, input.with_d, input.j, output.fasta)

rule igdiscover_allele_alignments_align_d:
    """Align discovered D alleles to target sequences based on V alignment."""
    output: fasta="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/D.aligned.fasta"
    input:
        target="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}.aln.fa",
        with_v="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/V.aligned.fasta",
        d="analysis/igdiscover/{chain}.{chain_type}/{specimen}/final/database/D.fasta"
    run:
        if wildcards.chain == "heavy":
            align_next_segment(
                input.target, input.with_v, input.d, output.fasta)
        else:
            shell("cp {input.with_v} {output.fasta}")

rule igdiscover_allele_alignments_align_v:
    """Align discovered V alleles to target sequences.

    This will generally be the antibody sequences provided by
    gather_antibody_sequences, but if another FASTA is provided manually it'll
    work with that instead.
    """
    output: fasta="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}/V.aligned.fasta"
    input:
        target="analysis/reporting/{specimen}.{chain}.{chain_type}/{target}.aln.fa",
        v="analysis/igdiscover/{chain}.{chain_type}/{specimen}/final/database/V.fasta",
    shell:
        """
            if [[ -s {input.target} ]]; then
                clustalw -align -profile1={input.target} -profile2={input.v} -sequences -output=fasta -outfile={output}
            else
                touch {output}
            fi
        """

rule align_targets:
    """Align target sequences to one another before bringing discovered alleles into the mix."""
    output: "analysis/reporting/{specimen}.{chain}.{chain_type}/antibodies.{antibody_lineage}.aln.fa"
    input: "analysis/reporting/{specimen}.{chain}.{chain_type}/antibodies.{antibody_lineage}.fasta"
    shell:
        """
            if [[ -s {input} ]]; then
                clustalw -align -infile={input} -output=fasta -outfile={output}
            else
                touch {output}
            fi
        """

rule gather_antibody_sequences:
    """Gather mature antibody sequences to start off the alignment against."""
    output: "analysis/reporting/{specimen}.{chain}.{chain_type}/antibodies.{antibody_lineage}.fasta"
    run:
        gather_antibodies(
            wildcards.antibody_lineage, wildcards.chain, ANTIBODY_ISOLATES, output[0])







### Another approach: by lineage

rule lineage_gather_antibody_sequences:
    """Gather mature antibody sequences to start off the alignment against.

    This is the same as gather_antibody_sequences but handles the output differently.
    """
    output: "analysis/reporting/lineages/{subject}.{lineage}.{chain}/antibodies.fasta"
    run:
        gather_antibodies(
            wildcards.lineage, wildcards.chain, ANTIBODY_ISOLATES, output[0])

def lineage_gather_germline_input(w):
    """Input for lineage_gather_germline rule.

    This uses info from antibody lineage metadata to find the right germline
    allele FASTA files for the appropriate locus and segments for a given
    subject and chain."""
    # Pretty sure I us this logic in a dozen different places.  Should clean up.
    lineage_attrs = ANTIBODY_LINEAGES[w.lineage]
    if w.chain == "heavy":
        chain_type = "gamma"
        segments = ["V", "D", "J"]
    elif w.chain == "light":
        if "IGLV" in lineage_attrs["VL"]:
            chain_type = "lambda"
        elif "IGKV" in lineage_attrs["VL"]:
            chain_type = "kappa"
        else:
            raise ValueError("VL of %s not recognized" % lineage_attrs["VL"])
        segments = ["V", "J"]
    else:
        raise ValueError("heavy or light, not %s" % w.chain)
    fmt = "analysis/sonar/{subject}/{chain}.{chain_type}/germline.{segment}.fasta"
    mktarget = lambda s: fmt.format(
        subject=w.subject, chain=w.chain, chain_type=chain_type, segment=s)
    targets = {segment: mktarget(segment) for segment in segments}
    return targets

from Bio import SeqIO

rule lineage_gather_germline:
    """Get the germline alleles (across segments) matched to a given linage.

    This will create a FASTA file with at most three (V/D/J, for heavy) or two
    (V/J, for light) sequences.
    """
    output: "analysis/reporting/lineages/{subject}.{lineage}.{chain}/germline.fasta"
    input: unpack(lineage_gather_germline_input)
    run:
        lineage_attrs = ANTIBODY_LINEAGES[wildcards.lineage]
        with open(output[0], "w") as f_out:
            for key, val in input.items():
                if wildcards.chain == "heavy":
                    seqid = lineage_attrs[key + "H"]
                else:
                    seqid = lineage_attrs[key + "L"]
                with open(val) as f_in:
                    for record in SeqIO.parse(f_in, "fasta"):
                        if record.id == seqid:
                            SeqIO.write(record, f_out, "fasta")
                            break
