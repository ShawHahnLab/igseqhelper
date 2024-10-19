from csv import DictWriter, DictReader
from Bio import SeqIO

def parse_seq_desc(record):
    """Parse key=val pairs from a SeqRecord description to a dict.

    If a sring is given, take that as the description.  Any entries without "="
    are skipped.
    """
    try:
        txt = record.description
    except AttributeError:
        txt = record
    attrs = [pair.split("=", 1) for pair in txt.split(" ")]
    attrs = {pair[0]: pair[1] for pair in attrs if len(pair) == 2}
    return attrs

# Just a helper for the below rules.
# Filter out any entries that don't match the metadata since we rely on the
# metadata when prepping the tables.  (This could be things like
# "subject-(datestamp)" or what have you)
def filter_wildcards_by_metadata(wildcards_dict):
    removes = set()
    for key, vec in wildcards_dict.items():
        for idx, item in enumerate(vec):
            if (key == "subject" and item not in {row["Subject"] for row in SPECIMENS.values()}) or \
                (key == "chain_type" and item not in {row["Type"] for row in SAMPLES.values()}) or \
                (key == "specimen" and item not in SPECIMENS.keys()) or \
                (key == "antibody_lineage" and item not in ANTIBODY_LINEAGES.keys()):
                removes.add(idx)
    # (careful, need to remove consistently across all wildcards)
    wildcards_dict_out = {}
    for key, vec in wildcards_dict.items():
        wildcards_dict_out[key] = [item for idx, item in enumerate(vec) if idx not in removes]
    return wildcards_dict_out

AVAILABLE_SONAR_ISLANDS = filter_wildcards_by_metadata(dict(zip(
    ["subject", "chain_type", "specimen", "antibody_lineage", "dummy"],
    glob_wildcards("analysis/sonar/{subject}.{chain_type}/{specimen}/"
        "output/sequences/nucleotide/islandSeqs_{dummy}.fa"))))

AVAILABLE_SONAR_ANCESTORS = filter_wildcards_by_metadata(dict(zip(
    ["subject", "chain_type", "antibody_lineage", "dummy"],
    glob_wildcards("analysis/sonar/{subject}.{chain_type}/longitudinal.auto.{antibody_lineage}/"
        "output/sequences/nucleotide/longitudinal-{dummy}_inferredAncestors.fa"))))

AVAILABLE_IGDISCOVER = dict(zip(
    ["ref", "chain_type", "subject", "segment"],
    glob_wildcards("analysis/igdiscover/{ref}/{chain_type}/{subject}/final/database/{segment}.fasta")))

AVAILABLE_MININGD = filter_wildcards_by_metadata(dict(zip(
    ["subject"], glob_wildcards("analysis/mining-d/{subject}.output.default.fasta"))))

TARGET_REPORT_COUNTS = expand("analysis/reporting/counts/counts_by_{thing}.csv", thing=["sample", "run", "specimen"])
TARGET_REPORT_SONAR_ISLAND_SUMMARIES = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/island_stats_summary.csv", zip, **AVAILABLE_SONAR_ISLANDS)
TARGET_REPORT_SONAR_MEMBERS_TABLES = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/members.csv", zip, **AVAILABLE_SONAR_ISLANDS)
TARGET_REPORT_SONAR_IGPHYML_COLLECTED = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_collected.csv", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_ANCESTORS = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.csv", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_ancestors.common.csv", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_ALIGNED = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_aligned.fa", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_TREES = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_tree.tree", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_SONAR_IGPHYML_TREEPDFS = expand("analysis/reporting/sonar/{antibody_lineage}.{chain_type}/igphyml_tree.pdf", zip, **AVAILABLE_SONAR_ANCESTORS)
TARGET_REPORT_IGDISCOVER_TREES = expand("analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.nex", zip, **AVAILABLE_IGDISCOVER)
TARGET_REPORT_IGDISCOVER_FINAL_DBS = expand("analysis/reporting/igdiscover/{ref}/{chain_type}/{subject}/{segment}.fasta", zip, **AVAILABLE_IGDISCOVER)
TARGET_REPORT_MININGD_FASTAS = expand("analysis/reporting/mining-d/{subject}/{subject}.fasta", zip, **AVAILABLE_MININGD)
TARGET_REPORT_MININGD_TREES = expand("analysis/reporting/mining-d/{subject}/{subject}.nex", zip, **AVAILABLE_MININGD)

rule all_report:
    input: TARGET_REPORT_COUNTS +
        TARGET_REPORT_SONAR_ISLAND_SUMMARIES +
        TARGET_REPORT_SONAR_MEMBERS_TABLES +
        TARGET_REPORT_SONAR_IGPHYML_COLLECTED +
        TARGET_REPORT_SONAR_IGPHYML_ANCESTORS +
        TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON +
        TARGET_REPORT_SONAR_IGPHYML_ALIGNED +
        TARGET_REPORT_SONAR_IGPHYML_TREES +
        TARGET_REPORT_SONAR_IGPHYML_TREEPDFS +
        TARGET_REPORT_IGDISCOVER_TREES +
        TARGET_REPORT_IGDISCOVER_FINAL_DBS +
        TARGET_REPORT_MININGD_TREES +
        ["analysis/reporting/mining-d/all.nex"]

rule report_available_sonar:
    input: TARGET_REPORT_SONAR_ISLAND_SUMMARIES +
        TARGET_REPORT_SONAR_MEMBERS_TABLES +
        TARGET_REPORT_SONAR_IGPHYML_COLLECTED +
        TARGET_REPORT_SONAR_IGPHYML_ANCESTORS +
        TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON +
        TARGET_REPORT_SONAR_IGPHYML_ALIGNED +
        TARGET_REPORT_SONAR_IGPHYML_TREES +
        TARGET_REPORT_SONAR_IGPHYML_TREEPDFS

rule report_available_igdiscover:
    input: TARGET_REPORT_IGDISCOVER_TREES +
        TARGET_REPORT_IGDISCOVER_FINAL_DBS

rule report_counts:
    input: TARGET_REPORT_COUNTS

rule report_available_sonar_island_summaries:
    input: TARGET_REPORT_SONAR_ISLAND_SUMMARIES

rule report_available_sonar_members_tables:
    input: TARGET_REPORT_SONAR_MEMBERS_TABLES

rule report_available_sonar_igphyml_collected:
    input: TARGET_REPORT_SONAR_IGPHYML_COLLECTED

rule report_available_sonar_igphyml_ancestors:
    input: TARGET_REPORT_SONAR_IGPHYML_ANCESTORS

rule report_available_sonar_igphyml_ancs_common:
    input: TARGET_REPORT_SONAR_IGPHYML_ANCS_COMMON

rule report_available_sonar_igphyml_aligned:
    input: TARGET_REPORT_SONAR_IGPHYML_ALIGNED

rule report_available_sonar_igphyml_trees:
    input: TARGET_REPORT_SONAR_IGPHYML_TREES

rule report_available_sonar_igphyml_treepdfs:
    input: TARGET_REPORT_SONAR_IGPHYML_TREEPDFS

rule report_available_igdiscover_trees:
    input: TARGET_REPORT_IGDISCOVER_TREES

rule report_available_igdiscover_final_dbs:
    input: TARGET_REPORT_IGDISCOVER_FINAL_DBS

rule report_available_miningd_fastas:
    input: TARGET_REPORT_MININGD_FASTAS

rule report_available_miningd_trees:
    input: TARGET_REPORT_MININGD_TREES
