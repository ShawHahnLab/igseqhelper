"""
V(D)J germline allele assignment using IgDiscover'd specimen's results, for
mature antibodies belonging to the same subject.  This can help support the
manually-entered VH/DH/JH/VL/JL values in the Antibody Lineages table.

The final output here is a small CSV file per locus per specimen tallying
counts of assignments to each allele across the antibodies.
"""

rule igdiscover_igblast_summary:
    """Condense IgBLAST results for antibody seqs into a small CSV.

    There will be one or more rows each for V, D, and J (for both heavy and
    light) showing the total count of occurrences for each allele for each
    segment across the antibody sequences.  The total count for a segment can
    go higher than the count of antibody sequences if multiple candidate allele
    assignments are listed for some antibodies.
    """
    output: "analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/assignments.csv"
    input: "analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/antibodies.tsv"
    run:
        from csv import DictReader, DictWriter
        from collections import defaultdict
        cts = {
            "V": defaultdict(int),
            "D": defaultdict(int),
            "J": defaultdict(int)}
        with open(input[0]) as f_in:
            reader = DictReader(f_in, delimiter="\t")
            for row in reader:
                for allele in row["v_call"].split(","):
                    cts["V"][allele] += 1
                for allele in row["d_call"].split(","):
                    cts["D"][allele] += 1
                for allele in row["j_call"].split(","):
                    cts["J"][allele] += 1
        with open(output[0], "wt") as f_out:
            writer = DictWriter(f_out, fieldnames=["Segment", "Allele", "Count"], lineterminator="\n")
            writer.writeheader()
            for segment in cts:
                for allele in cts[segment]:
                    writer.writerow({
                        "Segment": segment,
                        "Allele": allele,
                        "Count": cts[segment][allele]})

rule igdiscover_igblast_run:
    """Run IgBLAST against this specimen's custom germline sequences.

    The query could be anything, but here we only have rules in place to
    provide an antibody fasta for antibodies belonging to this subject's
    lineages as defined in the metadata.
    """

    output: "analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/{query}.tsv"
    input:
        query="analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/{query}.fasta",
        db_v="analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/database/V.nhr",
        db_d="analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/database/D.nhr",
        db_j="analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/database/J.nhr"
    params:
        outfmt=19, # AIRR TSV format
        organism="rhesus_monkey"
    shell:
        """
            dbv={input.db_v}
            dbv=${{dbv%.nhr}}
            dbd={input.db_d}
            dbd=${{dbd%.nhr}}
            dbj={input.db_j}
            dbj=${{dbj%.nhr}}
            igblastn \
                -germline_db_V $dbv \
                -germline_db_D $dbd \
                -germline_db_J $dbj \
                -outfmt {params.outfmt} \
                -organism {params.organism} \
                -ig_seqtype Ig \
                -query {input.query} \
                -out {output}
        """

rule igdiscover_igblast_db:
    """Set up an IgBLAST database using the germline sequences inferred from this specimen."""
    output: "analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/database/{segment}.nhr"
    input: "analysis/igdiscover/{chain}.{chain_type}/{specimen}/final/database/{segment}.fasta",
    params: out="analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/database/{segment}"
    shell: "makeblastdb -dbtype nucl -parse_seqids -in {input} -out {params.out}"

rule igdiscover_igblast_query_antibodies:
    """Gather all antibody sequences for this specimen's subject and this chain."""
    output:
        query="analysis/igdiscover-igblast/{chain}.{chain_type}/{specimen}/antibodies.fasta"
    run:
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        subject = SPECIMENS[wildcards.specimen]["Subject"]
        if wildcards.chain == "heavy":
            key = "HeavySeq"
        elif wildcards.chain  == "light":
            key = "LightSeq"
        else:
            raise ValueError
        with open(output.query, "wt") as f_out:
            for mab_attrs in ANTIBODY_ISOLATES.values():
                # include any that match this subject and have a non-empty
                # sequence (sometimes we might only have one chain or the other
                # for any particular isolate)
                if mab_attrs["AntibodyLineageAttrs"]["Subject"] == subject and mab_attrs[key]:
                    SeqIO.write(SeqRecord(
                            seq=Seq(mab_attrs[key]),
                            id=mab_attrs["AntibodyIsolate"],
                            description=""),
                        f_out, "fasta-2line")
