"""
Report generation and helper rules.

These prepare summary spreadsheets and plots from various steps in the
analysis.  It's a bit of a mess but the idea is everything funnels toward
analysis/report.pdf, with as much as possible handled in a modular fashion here
before going in the Rmarkdown for the report.
"""

from igseqhelper.reporting.sonar import (
    sonar_island_summary,
    sonar_island_stats)
from igseqhelper.reporting.counts import (
    counts_sample_summary,
    counts_run_summary,
    counts_specimen_summary,
    counts_presto_specimen_summary,
    counts_presto_assembly_summary,
    counts_presto_qual_summary,
    counts_sonar_module1_summary)
from igseqhelper.data import amplicon_files, transpose_sample_md

SAMPLE_MD_IGM = transpose_sample_md(SAMPLES, "IgM+")

# Only for samples where we have the run info
TARGET_QUALTRIM_GRID = expand(
    outputs_per_run(
        "analysis/reporting/by-run/{run}/qualtrim.{sample}.{{rp}}.csv",
        {key: SAMPLES[key] for key in SAMPLES if SAMPLES[key]["Run"]}),
    rp=["R1", "R2", "I1"])

TARGET_BARCODE_SUMMARY = expand(
    "analysis/reporting/by-run/{run}/barcode_summary.csv", run=[runattrs["Run"] for runattrs in RUNS.values() if runattrs["Protocol"] == "IgSeq"])

TARGET_SONAR_ISLAND_STATS = expand(
    "analysis/reporting/by-lineage/{antibody_lineage}/{specimen}.{chain}.{chain_type}/island_stats.csv",
    zip, **igseqhelper.sonar.setup_sonar_combos(
        transpose_sample_md({key: SAMPLES[key] for key in SAMPLES if SAMPLES[key]["Run"]}, "IgG+"),
        ANTIBODY_LINEAGES))

rule all_qualtrim_grid:
    input: TARGET_QUALTRIM_GRID

rule all_barcode_summary:
    input: TARGET_BARCODE_SUMMARY

rule all_sonar_island_stats:
    input: TARGET_SONAR_ISLAND_STATS

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
