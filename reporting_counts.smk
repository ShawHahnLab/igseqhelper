"""
Rules for reporting on totals of reads and such at different steps in the workflow.
"""

import lzma

rule report_counts_by_sample:
    """Sample read count summary CSV (one row per sample)"""
    output: "analysis/reporting/counts/counts_by_sample.csv"
    input:
        samples=ancient("metadata/samples.csv"),
        specimens=ancient("metadata/specimens.csv")
    shell: "report_counts_by_sample.py {input.samples} {input.specimens} {output}"

rule report_counts_by_run:
    """Run read count summary CSV (one row per run)"""
    output: "analysis/reporting/counts/counts_by_run.csv"
    input: "analysis/reporting/counts/counts_by_sample.csv"
    shell: "report_counts_by_run.py {input} {output}"

rule report_counts_by_specimen:
    """Specimen read and SONAR cluster summary CSV (one row per specimen per cell+chain type)"""
    output: "analysis/reporting/counts/counts_by_specimen.csv"
    input: "analysis/reporting/counts/counts_by_sample.csv"
    shell: "report_counts_by_specimen.py {input} {output}"

rule report_counts_by_subject:
    """Subject read and SONAR cluster summary CSV (one row per subject per cell+chain type)"""
    output: "analysis/reporting/counts/counts_by_subject.csv"
    input: "analysis/reporting/counts/counts_by_specimen.csv"
    shell: "report_counts_by_subject.py {input} {output}"

def sonar_airr_counts(input_tsv, output_csv, wcards):
    fmt_islands = ("analysis/sonar/{subject}.{chain_type}/"
        "{specimen}/output/tables/islandSeqs_{{antibody_lineage}}.txt").format(**wcards)
    lineages = [lin for lin, attrs in ANTIBODY_LINEAGES.items() if \
        attrs["Subject"] == wcards["subject"]]
    # Reads: grand total in the file.  Some don't make it this far if they
    #        didn't even look like antibody reads (e.g. ferritin)
    # GoodReads: The good (no stops, all parts intact) antibody reads
    # ClusteredReads: Good ones that made it into clusters (not quite all do)
    # ClusteredUnique The *number* of clusters
    # LineageMembers: sum of lineage members across any lineages for this
    #                 specimen (only included if files are available)
    counts = {"SONARReads": 0, "SONARGoodReads": 0, "SONARClusteredReads": 0, "SONARClusteredUnique": 0}
    opener = {".tsv": open, ".xz": lzma.open}.get(Path(input_tsv).suffix, open)
    with opener(input_tsv, "rt", encoding="ASCII") as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        for row in reader:
            counts["SONARReads"] += int(row["duplicate_count"])
            if row["status"] == "good":
                counts["SONARGoodReads"] += int(row["duplicate_count"])
                if row["cluster_count"]:
                    # counting the number of reads counted in clusters.  Ony good
                    # reads get clustered, and only those with the minimum dupliate
                    # count (by default singletons are not assigned clusters, good
                    # or not).
                    counts["SONARClusteredReads"] += int(row["cluster_count"])
                    # counting each cluster once for this one
                    counts["SONARClusteredUnique"] += 1
    # If there are islandSeqs files for the expected lineages, tally those up
    # too (but just if available)
    if lineages and fmt_islands:
        members = None
        for antibody_lineage in lineages:
            path = Path(fmt_islands.format(antibody_lineage = antibody_lineage))
            if path.exists():
                members = members or 0
                with open(path) as f_in:
                    for line in f_in:
                        members += 1
        if members is not None:
            counts["LineageMembers"] = members
    with open(output_csv, "wt") as f_out:
        writer = DictWriter(f_out, fieldnames=counts.keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerow(counts)

# Typically I wouldn't want to base rules on what happens to be on disk, but
# reporting is a special case.  (I'm not tying these targets into downstream
# summary rules, though, because I don't want those to ever trigger rebuilding
# actual analysis outputs.)
# TODO figure out how to enforce my global wildcard constraints here
# https://stackoverflow.com/questions/79744378
rule report_available_sonar_airr_counts:
    """"Tabulate SONAR AIRR counts CSV for whatever rearrangements.tsv files are on hand"""
    input: expand("analysis/reporting/counts/sonar/{subject}.{chain_type}.{specimen}.csv", zip, **glob_wildcards("analysis/sonar/{subject}.{chain_type,[a-z]+}/{specimen,[A-Za-z0-9]+}/output/tables/{specimen2}_rearrangements.tsv{suffix,[^/]*}")._asdict())

ruleorder: report_sonar_airr_counts_xz > report_sonar_airr_counts

rule report_sonar_airr_counts:
    output: csv="analysis/reporting/counts/sonar/{subject}.{chain_type}.{specimen}.csv"
    input: tsv="analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv"
    run:
        sonar_airr_counts(input.tsv, output.csv, vars(wildcards))

rule report_sonar_airr_counts_xz:
    output: csv="analysis/reporting/counts/sonar/{subject}.{chain_type}.{specimen}.csv"
    input: tsv="analysis/sonar/{subject}.{chain_type}/{specimen}/output/tables/{specimen}_rearrangements.tsv.xz"
    run:
        sonar_airr_counts(input.tsv, output.csv, vars(wildcards))

rule report_sonar_airr_counts_by_subject_and_type:
    output: touch("analysis/reporting/counts/sonar/{subject}.{chain_type}.done")
    input: lambda w: expand("analysis/reporting/counts/sonar/{{subject}}.{{chain_type}}.{specimen}.csv", specimen={samp["Specimen"] for samp in SAMPLES.values() if samp["SpecimenAttrs"]["Subject"] == w.subject and samp["Type"] == w.chain_type})
