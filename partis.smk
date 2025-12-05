"""
Run partis on all IgG+ material across specimens for one subject and chain type
(using SONAR's clustered sequences), and partition all of those sequences plus
the corresponding mature isolate sequences.

Some jargon used in the rules and scripts here:

 * Category: "isolate" or "ngs"
 * Lineage Group: automatically grouped lineage/clone identifiers coming from
   my own lineage assignments and/or partis' automatic groupings
"""

def input_for_partis_ngs_fasta(w):
    """List specimen names that should be used for partis for one subject+chain type"""
    specs = []
    for spec_attrs in SPECIMENS.values():
        if spec_attrs["Subject"] == w.subject and "IgG+" in spec_attrs["CellType"]:
            specs.append(spec_attrs["Specimen"])
    for samp_attrs in SAMPLES.values():
        spec_attrs = samp_attrs.get("SpecimenAttrs", {})
        if samp_attrs["Type"] == chain_type and \
                spec_attrs.get("Subject") == w.subject and \
                "IgG+" in spec_attrs.get("CellType"):
            specs.append(spec_attrs["Specimen"])
    targets = expand(
        "analysis/sonar/{subject}.{chain_type}/{specimen}/output/sequences/nucleotide/{specimen}_goodVJ_unique.fa",
        subject = w.subject,
        chain_type = w.chain_type,
        specimen = specs)
    return dict(zip(specs, targets))

rule partis_ngs_fasta:
    """Cross-timepoint NGS seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/ngs.fasta"
    input: unpack(input_for_partis_ngs_fasta)
    run:
        with open(output[0], "w") as f_out:
            for spec, fasta in input.items():
                for rec in SeqIO.parse(fasta, "fasta"):
                    rec.description = ""
                    rec.id = "specimen-" + spec+ "-" + rec.id
                    SeqIO.write(rec, f_out, "fasta-2line")

rule partis_isolates_fasta:
    """Cross-lineage isolate seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/isolates.fasta"
    run:
        seq_col = "LightSeq" if wildcards.chain_type in ["kappa", "lambda"] else "HeavySeq"
        with open(output[0], "w") as f_out:
            for isolate_attrs in ANTIBODY_ISOLATES.values():
                if isolate_attrs[seq_col]:
                    lineage_attrs = isolate_attrs.get("LineageAttrs", {})
                    subject = lineage_attrs.get("Subject")
                    if subject and subject == wildcards.subject:
                        seqid = "isolate-" + isolate_attrs["Isolate"]
                        seq = isolate_attrs[seq_col]
                        f_out.write(f">{seqid}\n{seq}\n")

rule partis_seqsets_fasta:
    """All seqset seqs for this subject and locus for partis"""
    output: "analysis/partis/{subject}.{chain_type}/seqsets.fasta"
    run:
        locus = {"kappa": "IGK", "lambda": "IGL"}.get(wildcards.chain_type, "IGH")
        # (noting isolate seqs so we can exclude overlaps where I've defined an
        # isolate *from* the seqset input)
        chain_col = "HeavySeq" if locus == "IGH" else "LightSeq"
        isols_by_seq = defaultdict(set)
        for isol, attrs in ANTIBODY_ISOLATES.items():
            linattrs = attrs["LineageAttrs"]
            if linattrs["Subject"] == wildcards.subject:
                isols_by_seq[attrs[chain_col]].add(isol)
                if attrs["AltName"]:
                    isols_by_seq[attrs[chain_col]].add(attrs["AltName"])
        recs = []
        for seqset, attrs in SEQSETS.items():
            path = f"analysis/seqsets/{seqset}.{locus}.fasta.gz"
            tpoint = attrs["Timepoint"]
            if attrs["Subject"] == wildcards.subject:
                with gzip.open(path, "rt", encoding="ASCII") as f_in:
                    for rec in SeqIO.parse(f_in, "fasta"):
                        seq = str(rec.seq)
                        # skip cases where the sequence matches an isolate
                        # *and* the timepoint-prefixed ID is contained within
                        # one of those isolate IDs
                        if f"wk{tpoint}-{rec.id}" in "/".join(isols_by_seq.get(seq, "")):
                            continue
                        seqid = f"seqset-{seqset}-{rec.id}"
                        recs.append((seqid, seq))
        with open(output[0], "w") as f_out:
            for seqid, seq in recs:
                f_out.write(f">{seqid}\n{seq}\n")

rule partis_combo_fasta:
    """Combined isolate+seqset+NGS seqs for partis"""
    output: "analysis/partis/{subject}.{chain_type}/combined.fasta"
    input:
        isolates="analysis/partis/{subject}.{chain_type}/isolates.fasta",
        ngs="analysis/partis/{subject}.{chain_type}/ngs.fasta",
        seqsets="analysis/partis/{subject}.{chain_type}/seqsets.fasta"
    shell: "cat {input} > {output}"

rule partis_excludes:
    """Gather sequences to exclude from partis parameter caching"""
    output:
        includes="analysis/partis/{subject}.{chain_type}/combined.includes.fasta",
        excludes="analysis/partis/{subject}.{chain_type}/combined.excludes.fasta"
    input: "analysis/igblast/custom-{subject}.IGH/partis/{subject}.{chain_type}/combined.fasta.tsv.gz"
    run:
        with \
                open(output.excludes, "w") as f_out_excl, \
                open(output.includes, "w") as f_out_incl, \
                gzip.open(input[0], "rt") as f_in:
            for row in DictReader(f_in, delimiter="\t"):
                seqid = row["sequence_id"]
                seq = row["sequence"]
                vseq = row["v_sequence_alignment"].replace("-", "")
                hndl = f_out_incl
                # I've had trouble with sequences with a big deletion in V
                # and/or a seemingly-absent junction causing partis to crash.
                # Those aren't going to be of interest for us anyway so I'll
                # just exclude them.
                if len(vseq) < 250 and (row["complete_vdj"] == "T" or not row["junction"]):
                    hndl = f_out_excl
                hndl.write(f">{seqid}\n{seq}\n")

rule partis_cache_params_fasta:
    """Make a version of the repertoire NGS FASTA for use with partis cache-parameters"""
    output: "analysis/partis/{subject}.{chain_type}/ngs.filt.fasta"
    input:
        fasta_ngs="analysis/partis/{subject}.{chain_type}/ngs.fasta",
        fasta_excludes="analysis/partis/{subject}.{chain_type}/combined.excludes.fasta"
    run:
        excludes = {rec.id: str(rec.seq) for rec in SeqIO.parse(input.fasta_excludes, "fasta")}
        with open(output[0], "w") as f_out:
            for rec in SeqIO.parse(input.fasta_ngs, "fasta"):
                if rec.id in excludes:
                    if str(rec.seq) != excludes[rec.id]:
                        raise ValueError
                else:
                    f_out.write(f">{rec.id}\n{rec.seq}\n")

def input_for_partis_germline(w):
    locus = w.locus_lower.upper()
    return {
        key: f"analysis/germline/{w.subject}.{locus}/{key}.fasta" for key in ("V", "D", "J")}

rule partis_germline:
    """Prep partis germline dir from per-subject germline files"""
    # This is just a bit of reformatting of the per-subject IgDiscover-based
    # germline files we have from elsewhere
    output:
        out_dir=directory("analysis/partis/germlines/{subject}/{locus_lower}"),
        v="analysis/partis/germlines/{subject}/{locus_lower}/{locus_lower}v.fasta",
        # D is created by the script only for heavy; the D input will be
        # ignored for light chains
        j="analysis/partis/germlines/{subject}/{locus_lower}/{locus_lower}j.fasta",
        # CSV of positions of codons for conserved AAs around CDR3
        extras="analysis/partis/germlines/{subject}/{locus_lower}/extras.csv"
    input: unpack(input_for_partis_germline)
    params:
        out_parent="analysis/partis/germlines/{subject}"
    shell: "partis_germline.py {wildcards.locus_lower} {input.V} {input.D} {input.J} {params.out_parent}"

def input_for_partis_cache_params(w):
    locus = {"kappa": "igk", "lambda": "igl"}.get(w.chain_type, "igh")
    return {
        "seqs": f"analysis/partis/{w.subject}.{w.chain_type}/ngs.filt.fasta",
        "germline": f"analysis/partis/germlines/{w.subject}/{locus}/extras.csv"}

rule partis_cache_params:
    """Run partis cache-parameters on a cross-timepoint NGS FASTA"""
    output: directory("analysis/partis/{subject}.{chain_type}/params")
    input: unpack(input_for_partis_cache_params)
    log:
        main="analysis/partis/{subject}.{chain_type}/cache_params.log.txt",
        conda="analysis/partis/{subject}.{chain_type}/cache_params.conda_build.txt"
    params:
        locus=lambda w: {"kappa": "igk", "lambda": "igl"}.get(w.chain_type, "igh"),
        species=config.get("partis_species", "macaque"),
        seed=config.get("partis_random_seed", 1),
        germline="analysis/partis/germlines/{subject}",
        leave_default_germline=config.get("partis_leave_default_germline", True),
        # I have the dependencies provided via conda but the actual partis
        # install still lives in a big ball of stuff in one directory,
        # unfortunately.  The first like of the shell commands will ensure this
        # variable is set.
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    threads: 14
    shell:
        """
            if [[ ! -n "{params.partis}" ]]; then
              echo "Need path to partis install"
              exit 1
            fi
            (
              date
              echo "PARTIS_HOME: {params.partis}"
              echo "locus: {params.locus}"
              echo "species: {params.species}"
              echo "germline: {params.germline}"
              echo "leave default germline: {params.leave_default_germline}"
              echo "random seed: {params.seed}"
            ) >> {log.main}
            germ_arg=""
            if [[ "{params.leave_default_germline}" == "True" ]]; then
              germ_arg="--leave-default-germline"
            fi
            {params.partis}/bin/partis cache-parameters --n-procs {threads} \
              --initial-germline-dir {params.germline} "$germ_arg" \
              --random-seed {params.seed} --locus {params.locus} --species {params.species} \
              --infname {input.seqs} --parameter-dir {output} 2>&1 | tee -a {log.main}
        """

rule partis_partition:
    """Run partis partition on a combined isolate+NGS FASTA with existing cached parameters"""
    output: "analysis/partis/{subject}.{chain_type}/partitions.yaml"
    input:
        params_dir="analysis/partis/{subject}.{chain_type}/params",
        fasta="analysis/partis/{subject}.{chain_type}/combined.includes.fasta"
    log:
        main="analysis/partis/{subject}.{chain_type}/partition.log.txt",
        conda="analysis/partis/{subject}.{chain_type}/partition.conda_build.txt"
    params:
        locus=lambda w: {"kappa": "igk", "lambda": "igl"}.get(w.chain_type, "igh"),
        species=config.get("partis_species", "macaque"),
        seed=config.get("partis_random_seed", 1),
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    threads: 14
    shell:
        """
            if [[ ! -n "{params.partis}" ]]; then
              echo "Need path to partis install"
              exit 1
            fi
            conda list --explicit > {log.conda}
            (
              date
              echo "PARTIS_HOME: {params.partis}"
              echo "locus: {params.locus}"
              echo "species: {params.species}"
              echo "random seed: {params.seed}"
            ) >> {log.main}
            {params.partis}/bin/partis partition --n-procs {threads} --no-naive-vsearch \
              --random-seed {params.seed} --locus {params.locus} --species {params.species} \
              --infname {input.fasta} --parameter-dir {input.params_dir} --outfname {output} 2>&1 | tee -a {log.main}
        """

rule partis_partition_airr:
    """Convert Partis' YAML output for partition data to AIRR TSV"""
    output: "analysis/partis/{subject}.{chain_type}/partitions.airr.tsv"
    input: "analysis/partis/{subject}.{chain_type}/partitions.yaml"
    log:
        main="analysis/partis/{subject}.{chain_type}/partition.log.txt",
        conda="analysis/partis/{subject}.{chain_type}/partition.conda_build.txt"
    params:
        partis=os.getenv("PARTIS_HOME", "")
    conda: "envs/partis.yaml"
    shell:
        """
            if [[ ! -n "{params.partis}" ]]; then
              echo "Need path to partis install"
              exit 1
            fi
            conda list --explicit > {log.conda}
            (
              date
              echo "PARTIS_HOME: {params.partis}"
            ) >> {log.main}
            {params.partis}/bin/parse-output.py --airr-output {input} {output} 2>&1 | tee -a {log.main}
        """

### Below is really reporting logic, but, no time for proper organization here

rule isolates_light_fasta:
    output: temp("analysis/partis/isolates_light.fasta")
    run:
        with open(output[0], "w", encoding="ASCII") as f_out:
            for isolate, attrs in ANTIBODY_ISOLATES.items():
                if attrs["LightSeq"]:
                    f_out.write(f">{isolate}\n{attrs['LightSeq']}\n")

def input_for_partis_seq_lineage_info(w):
    if w.chain_type in ("kappa", "lambda"):
        raise ValueError("light chain not supported")
    # always required the final AIRR table as input
    path = Path(expand("analysis/partis/{subject}.{chain_type}", **w)[0])
    # ...and IgBLAST's AIRR too
    igblast_path = Path(expand("analysis/igblast/custom-{subject}.IGH/partis/{subject}.{chain_type}", **w)[0])
    targets = {
        "airr": path/"partitions.airr.tsv",
        "airr_igblast": igblast_path/"combined.fasta.tsv.gz",
        "airr_isolates_light": "analysis/igblast/sonarramesh/partis/isolates_light.fasta.tsv.gz",
        "specimens": "metadata/specimens.csv",
        "seqsets": "metadata/seqsets.csv",
        "isolates": "metadata/isolates.csv"}
    # If there's a CSV provided with NGS sequence lineage info, use that also,
    # so we can assign custom lineage IDs to NGS seqs
    path_annotations = path/"ngs_lineages.csv"
    if path_annotations.exists():
        targets["ngs_annots"] = path_annotations
    return targets

rule partis_seq_lineage_info:
    """Report sequence and lineage info from partis, merging metadata for our seqs+isolates"""
    output: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    input: unpack(input_for_partis_seq_lineage_info)
    priority: 10
    run:
        cmd = "partis_seq_lineage_info.py {input.airr} {output} --metadata-isolates {input.isolates} --metadata-specimens {input.specimens} --metadata-seqsets {input.seqsets} -A {input.airr_igblast} -L {input.airr_isolates_light} --all"
        if "ngs_annots" in dict(input):
            cmd += " -n {input.ngs_annots}"
        shell(cmd)

rule partis_lineages:
    """Summarize partis info per-lineage-group, one row per group+timepoint+category"""
    output: "analysis/partis/{subject}.{chain_type}/lineage_groups.csv"
    input: "analysis/partis/{subject}.{chain_type}/seq_lineages.csv"
    priority: 10
    shell: "partis_lineages.py {input} {output}"

rule partis_lineages_summary:
    """Summarize partis info per-lineage-group further, one row per lineage group"""
    output: "analysis/partis/{subject}.{chain_type}/lineage_groups_summary.csv"
    input: "analysis/partis/{subject}.{chain_type}/lineage_groups.csv"
    priority: 10
    shell: "partis_lineages_summary.py {input} {output}"
