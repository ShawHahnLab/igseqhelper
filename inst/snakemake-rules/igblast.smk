rule igblast_run:
    """A generic rule to igblast any query against any reference.
    
    Recognized names for references (like igdiscover results) and queries (like
    SONAR's filtered, clustered sequences) will let this find the reference and
    query files automatically.  Otherwise, manually place reference and/or
    query FASTA files in the right place first:

    analysis/igblast/{ref}/{querytype}/{query}.fasta
    analysis/igblast/{ref}/db/V.fasta
    analysis/igblast/{ref}/db/D.fasta
    analysis/igblast/{ref}/db/J.fasta

    example targets:
       analysis/igblast/igdiscover.RMT681WK72.heavy.mu/sonar.T681.heavy.gamma/RMT681WK24IGG.tsv
       analysis/igblast/sonarramesh.H/sonar.T681.heavy.gamma/RMT681WK24IGG.tsv
       analysis/igblast/igdiscover.RMT681WK72.heavy.mu/lineage.heavy/T681.tsv

    There's a bit of redundancy in things like the chain getting repeated, but
    it's easier that way with the generalized rules here.
    """
    output: "analysis/igblast/{ref}/{querytype}/{query}.tsv"
    input:
        query="analysis/igblast/{ref}/{querytype}/{query}.fasta",
        db_v="analysis/igblast/{ref}/db/V.nhr",
        db_d="analysis/igblast/{ref}/db/D.nhr",
        db_j="analysis/igblast/{ref}/db/J.nhr"
    params:
        outfmt=19, # AIRR TSV format
        organism="rhesus_monkey",
        db_v="analysis/igblast/{ref}/db/V",
        db_d="analysis/igblast/{ref}/db/D",
        db_j="analysis/igblast/{ref}/db/J"
    shell:
        """
            igblastn \
                -germline_db_V {params.db_v} \
                -germline_db_D {params.db_d} \
                -germline_db_J {params.db_j} \
                -outfmt {params.outfmt} \
                -organism {params.organism} \
                -ig_seqtype Ig \
                -query {input.query} \
                -out {output}
        """

### Queries

def input_for_igblast_query_fasta_sonar_inferred(w):
    return expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/"
        "longitudinal/output/sequences/nucleotide/longitudinal_inferredAncestors.fa",
        subject = ANTIBODY_LINEAGES[w.antibody_lineage]["Subject"],
        chain = w.chain,
        chain_type = w.chain_type,
        antibody_lineage = w.antibody_lineage)

rule igblast_query_fasta_sonar_inferred:
    """Set up query FASTA with SONAR's inferred ancestors for one lineage+chain+type."""
    output: "analysis/igblast/{ref}/sonar.{antibody_lineage}.{chain}.{chain_type}/longitudinal.inferredAncestors.fasta"
    input: input_for_igblast_query_fasta_sonar_inferred
    shell: "cp {input} {output}"

def input_for_igblast_query_fasta_sonar_island(w):
    return expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/"
        "{specimen}/output/sequences/nucleotide/{specimen}_islandSeqs.fa",
        subject = SPECIMENS[w.specimen]["Subject"],
        chain = w.chain,
        chain_type = w.chain_type,
        antibody_lineage = w.antibody_lineage,
        specimen = w.specimen)

rule igblast_query_fasta_sonar_island:
    """Set up query FASTA with SONAR's goodVJ unique for one lineage+chain+type+specimen."""
    output: "analysis/igblast/{ref}/sonar.{antibody_lineage}.{chain}.{chain_type}/{specimen}.island.fasta"
    input: input_for_igblast_query_fasta_sonar_island
    shell: "cp {input} {output}"

def input_for_igblast_query_fasta_sonar(w):
    return expand(
        "analysis/sonar/{subject}/{chain}.{chain_type}/{antibody_lineage}/"
        "{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        subject = SPECIMENS[w.specimen]["Subject"],
        chain = w.chain,
        chain_type = w.chain_type,
        antibody_lineage = w.antibody_lineage,
        specimen = w.specimen)

rule igblast_query_fasta_sonar:
    """Set up query FASTA with SONAR's goodVJ unique for one lineage+chain+type+specimen."""
    output: "analysis/igblast/{ref}/sonar.{antibody_lineage}.{chain}.{chain_type}/{specimen}.fasta"
    input: input_for_igblast_query_fasta_sonar
    shell: "cp {input} {output}"

rule igblast_query_fasta_antibodies:
    """Set up query FASTA with known mAbs for a given lineage and chain."""
    output: "analysis/igblast/{ref}/lineage.{chain}/{antibody_lineage}.fasta"
    run:
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        if wildcards.chain == "heavy":
            key = "HeavySeq"
        elif wildcards.chain  == "light":
            key = "LightSeq"
        else:
            raise ValueError
        with open(output[0], "wt") as f_out:
            for mab_attrs in ANTIBODY_ISOLATES.values():
                # include any that match this lineage and have a non-empty
                # sequence (sometimes we might only have one chain or the other
                # for any particular isolate)
                if mab_attrs["AntibodyLineage"] == wildcards.antibody_lineage and mab_attrs[key]:
                    SeqIO.write(SeqRecord(
                            seq=Seq(mab_attrs[key]),
                            id=mab_attrs["AntibodyIsolate"],
                            description=""),
                        f_out, "fasta-2line")

### References

rule igblast_ref:
    """Set up VDJ BLAST DB from reference FASTA."""
    output: "analysis/igblast/{ref}/db/{segment}.nhr"
    input: "analysis/igblast/{ref}/db/{segment}.fasta"
    params: out="analysis/igblast/{ref}/db/{segment}"
    shell: "makeblastdb -dbtype nucl -parse_seqids -in {input} -out {params.out}"

rule igblast_ref_fasta_igdiscover:
    """Set up reference FASTA from an individualized database."""
    # actually makes segment.nhr/nin/nog/nsd/nsi/nsq but we'll just use .nhr for the rules
    output: "analysis/igblast/igdiscover.{specimen}.{chain}.{chain_type}/db/{segment}.fasta"
    input: "analysis/igdiscover/{chain}.{chain_type}/{specimen}/final/database/{segment}.fasta"
    shell: "cp {input} {output}"

def input_for_igblast_ref_fasta_sonarramesh(w):
    locus = w.locus
    if w.segment == "D":
        locus = "H"
    return expand("SONAR/germDB/Ig{locus}{segment}_BU_DD.fasta", locus=locus, segment=w.segment)

rule igblast_ref_fasta_sonarramesh:
    """Set up reference FASTA from SONAR's copy of the Ramesh sequences."""
    output: "analysis/igblast/sonarramesh.{locus}/db/{segment}.fasta"
    input: input_for_igblast_ref_fasta_sonarramesh
    shell: "cp {input} {output}"
