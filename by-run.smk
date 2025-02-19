from collections import defaultdict

rule filt:
    output:
        fqgz="analysis/filt/{run}/{sample}.fastq.gz",
        counts="analysis/filt/{run}/{sample}.filt.counts.csv"
    input:
        fqgz="analysis/merge/{run}/{sample}.fastq.gz"
    log:
        main="analysis/filt/{run}/{sample}.log.txt",
        conda="analysis/filt/{run}/{sample}.fastq.gz.conda_build.txt"
    params:
        qmin=config.get("filt_qmin", 10)
    conda: "envs/igseq.yaml"
    shell:
        """
            conda list --explicit > {log.conda}
            date >> {log.main}
            echo "Quality minimum: {params.qmin}" >> {log.main}
            fastq_qual_min.py {input.fqgz} {output.fqgz} --countsfile {output.counts} -Q {params.qmin} 2>&1 | tee -a {log.main}
        """

rule merge:
    output:
        fqgz="analysis/merge/{run}{suffix}/{sample}.fastq.gz",
        counts="analysis/merge/{run}{suffix}/{sample}.merge.counts.csv",
        log="analysis/merge/{run}{suffix}/{sample}.pear.log"
    input:
        r1="analysis/trim/{run}{suffix}/{sample}.R1.fastq.gz",
        r2="analysis/trim/{run}{suffix}/{sample}.R2.fastq.gz"
    log:
        conda="analysis/merge/{run}{suffix}/{sample}.fastq.gz.conda_build.txt"
    conda: "envs/igseq.yaml"
    threads: 4
    shell:
        """
            conda list --explicit > {log.conda}
            igseq merge -t {threads} {input.r1} {input.r2}
        """

def input_for_trim(w):
    # special case for already-prepared files from elsewhere
    if w.run == "external":
        return []
    # special case for multiple samples for one barcode pair
    # (see by-run-locus-demux.smk)
    if w.suffix == ".grouped":
        return {
            "demux": "analysis/demux/{run}{suffix}",
            "samples": "analysis/demux/{run}.samples_grouped.csv"}
    # The usual case
    for samp, attrs in SAMPLES.items():
        if samp == w.sample:
            runid = attrs["Run"]
            if runid == w.run:
                return {
                "demux": "analysis/demux/{run}{suffix}",
                "samples": ancient("metadata/samples.csv")}
            else:
                raise ValueError(f"Sample {w.sample} has run {runid}, not {w.run}")
    raise ValueError(f"Sample {w.sample} not found for run {w.run}")

rule trim:
    output:
        r1="analysis/trim/{run}{suffix}/{sample}.R1.fastq.gz",
        r2="analysis/trim/{run}{suffix}/{sample}.R2.fastq.gz",
        report1="analysis/trim/{run}{suffix}/{sample}.cutadapt1.json",
        report2="analysis/trim/{run}{suffix}/{sample}.cutadapt2.json",
        counts="analysis/trim/{run}{suffix}/{sample}.trim.counts.csv"
    input: unpack(input_for_trim)
    log:
        main="analysis/trim/{run}{suffix}/{sample}.log.txt",
        conda="analysis/trim/{run}{suffix}/{sample}.trim.conda_build.txt"
    conda: "envs/igseq.yaml"
    threads: 4
    params:
        species=config.get("species", "rhesus"),
        min_length=config.get("trim_min_length", 50),
        qual_cutoff=config.get("trim_qual_cutoff", 15),
        nextseq_qual_cutoff=config.get("trim_nextseq_qual_cutoff", 15),
        is_nextseq=lambda w: RUNS.get(w.run, {}).get("SequencerModel") == "NextSeq"
    shell:
        """
            conda list --explicit > {log.conda}
            nextseq_args=""
            if [[ "{params.is_nextseq}" == "True" ]]; then
                nextseq_args="--nextseq-trim {params.nextseq_qual_cutoff}"
            fi
            (
              date
              echo "Species: {params.species}"
              echo "Min trimmed length: {params.min_length}"
              echo "Quality cutoff: {params.qual_cutoff}"
              echo "NextSeq quality cutoff: {params.nextseq_qual_cutoff}"
              echo "Is NextSeq: {params.is_nextseq}"
            ) >> {log.main}
            (
              igseq trim -t {threads} \
                  --min-length {params.min_length} \
                  --quality-cutoff {params.qual_cutoff} \
                  --samples {input.samples} \
                  --species {params.species} \
                  {input.demux}/{wildcards.sample}.R1.fastq.gz \
                  {input.demux}/{wildcards.sample}.R2.fastq.gz \
                  $nextseq_args 2>&1
            ) | tee -a {log.main}
        """

rule phix:
    output:
        bam="analysis/phix/{run}{suffix}/phix.bam",
        counts="analysis/phix/{run}{suffix}/phix.counts.csv"
    input:
        demux="analysis/demux/{run}{suffix}"
    log:
        conda="analysis/phix/{run}{suffix}/conda_build.txt"
    conda: "envs/igseq.yaml"
    threads: 4
    shell:
        """
            conda list --explicit > {log.conda}
            igseq phix -t {threads} {input}
        """

### By-run rules

def make_run_rules():
    """Generate dynamic per-sequencing-run rules where applicable

    Demultiplexing creates varying output files on a per-run basis, so to have
    the individual files available as outputs in the workflow, we can use a
    separate rule per run.  (The alternative would be a Snakemake checkpoint,
    but why re-infer the whole DAG only after running the rule when we actually
    do know up-front what the outputs are?)"""
    by_run = defaultdict(list)
    for samp_attrs in SAMPLES.values():
        if samp_attrs["Run"]:
            by_run[samp_attrs["Run"]].append(samp_attrs)
    for runid, samp_attrs_list in by_run.items():
        samp_names = [samp["Sample"] for samp in samp_attrs_list]
        bc_pairs = [(samp["BarcodeFwd"], samp["BarcodeRev"]) for samp in samp_attrs_list]
        if len(bc_pairs) == len(set(bc_pairs)):
            # The regular case; run igseq demux with the samples CSV as-is
            # (For the alternate case, where we have fewer unique barcode pairs
            # than samples, see by-run-locus-demux.smk)
            rule:
                name: f"demux_{runid}"
                output:
                    outdir=directory(f"analysis/demux/{runid}"),
                    fqgz=expand(
                        "analysis/demux/{run}/{samp}.{rp}.fastq.gz",
                        run=runid, samp=samp_names, rp=["R1", "R2", "I1"])
                input:
                    reads=expand(
                        "analysis/reads/{run}/Undetermined_S0_L001_{rp}_001.fastq.gz",
                        run=runid, rp=["I1", "R1", "R2"]),
                    samples=ancient("metadata/samples.csv")
                log:
                    conda=f"analysis/demux/{runid}/conda_build.txt"
                conda: "envs/igseq.yaml"
                shell:
                    """
                        conda list --explicit > {log.conda}
                        igseq demux --samples {input.samples} --details {output.outdir}/details.csv.gz {input.reads}
                    """
            # These are just helpers to group outputs from other rules by run ID
            phix_targets = expand(
                "analysis/phix/{run}/phix.bam", run=runid)
            merge_targets = expand(
                "analysis/merge/{run}/{sample}.fastq.gz", run=runid, sample=samp_names)
            filt_targets = expand(
                "analysis/filt/{run}/{sample}.fastq.gz",
                run=runid, sample=samp_names)
            trim_targets = expand(
                "analysis/trim/{run}/{sample}.{rp}.fastq.gz",
                run=runid, sample=samp_names, rp=["R1", "R2"])
            fastqc_args = {
                "prefix": "analysis/fastqc",
                "run": runid,
                "sample": samp_names,
                "rp": ["I1", "R1", "R2"],
                "suffix": "_fastqc.quals.csv"}
            fastqc_targets = \
                expand("{prefix}/reads/{run}/Undetermined_S0_L001_{rp}_001{suffix}", **fastqc_args) + \
                expand("{prefix}/demux/{run}/{sample}.{rp}{suffix}", **fastqc_args)
            fastqc_args["rp"] = ["R1", "R2"] # (no I1 for the rest)
            fastqc_targets += \
                expand("{prefix}/trim/{run}/{sample}.{rp}{suffix}", **fastqc_args) + \
                expand("{prefix}/merge/{run}/{sample}{suffix}", **fastqc_args) + \
                expand("{prefix}/filt/{run}/{sample}{suffix}", **fastqc_args)
            rule:
                name: f"proc_{runid}"
                input: phix_targets + trim_targets + merge_targets + filt_targets + fastqc_targets
            rule:
                name: f"phix_{runid}"
                input: phix_targets
            rule:
                name: f"trim_{runid}"
                input: trim_targets
            rule:
                name: f"merge_{runid}"
                input: merge_targets
            rule:
                name: f"filt_{runid}"
                input: filt_targets
        rule:
            name: f"getreads_{runid}"
            input: f"analysis/reads/{runid}"

make_run_rules()

rule getreads:
    output:
        outdir=directory("analysis/reads/{run}"),
        r1="analysis/reads/{run}/Undetermined_S0_L001_R1_001.fastq.gz",
        i1="analysis/reads/{run}/Undetermined_S0_L001_I1_001.fastq.gz",
        r2="analysis/reads/{run}/Undetermined_S0_L001_R2_001.fastq.gz"
    input: ancient("/seq/runs/{run}")
    log:
        conda="analysis/reads/{run}/conda_build.txt"
    conda: "envs/igseq.yaml"
    threads: 28
    shell:
        """
            conda list --explicit > {log.conda}
            igseq getreads -t {threads} --threads-load $(({threads}<4 ? {threads} : 4)) {input}
        """

### Samples grouped by subject or specimen
# Each directory will contain a set of symbolic links pointing to the
# individual sample files

def grouped_samples_input(w, pattern="analysis/filt/{runid}/{samp}.fastq.gz"):
    targets = []
    # clumsy way of matching e.g. "igg" to "IgG+CD20+"
    def cellmatch(ct_query, ct_ref):
        for ref in ct_ref.split("+"):
            if ct_query.lower() in ref.lower():
                return True
        return False
    celltype = dict(w).get("celltype")
    if w.thing == "subject":
        for samp, attrs in SAMPLES.items():
            if attrs["SpecimenAttrs"]["Subject"] == w.name and attrs["Type"] == w.chain_type:
                if not celltype:
                    raise ValueError("cell type needed if grouping by subject")
                if cellmatch(w.celltype, attrs["SpecimenAttrs"]["CellType"]):
                    runid = attrs["Run"] or "external"
                    targets.extend(expand(pattern, runid=runid, samp=samp))
    elif w.thing == "specimen":
        for samp, attrs in SAMPLES.items():
            if attrs["Specimen"] == w.name and attrs["Type"] == w.chain_type:
                if not celltype or cellmatch(w.celltype, attrs["SpecimenAttrs"]["CellType"]):
                    runid = attrs["Run"] or "external"
                    targets.extend(expand(pattern, runid=runid, samp=samp))
    else:
        raise ValueError("unrecognized wildcards for rule grouped_samples_input")
    if not targets:
        raise ValueError("no matching inputs for rule grouped_samples_input")
    return targets

# path_link -> path_real
# (Python 3.12 adds a walk_up argument to the Path class' relative_to that I
# think would maybe take care of this.)
def relative_link_target(path_real, path_link):
    # Trim away the matching parents
    parts_real = [part for part in str(path_real).split("/") if part]
    parts_link = [part for part in str(path_link).split("/") if part]
    for idx, pair in enumerate(zip(parts_real, parts_link)):
        if pair[0] != pair[1]:
            break
    remainder_real = parts_real[idx:]
    remainder_link = parts_link[idx:]
    # However deep the remaining part of the link is, that's how many .. we
    # have to climb to point to the real path
    suffix_real = "/".join(remainder_real)
    ups = "/".join(".." for _ in range(len(remainder_link)-1))
    link_target = f"{ups}/{suffix_real}"
    return link_target

def symlink_in(path_real, link_dir):
    """Make a relative symlink from a path to the same name in a directory."""
    link_dir = Path(link_dir)
    path_link = Path(link_dir)/Path(path_real).name
    link_target = relative_link_target(path_real, path_link)
    link_dir.mkdir(parents=True, exist_ok=True)
    path_link.symlink_to(link_target)

rule grouped_samples_by_celltype:
    output: directory("analysis/samples-by-{thing}/{celltype}/{name}.{chain_type}")
    input: lambda w: grouped_samples_input(w)
    run:
        for path in input:
            symlink_in(path, output[0])

rule grouped_samples:
    output: directory("analysis/samples-by-{thing}/{name}.{chain_type}")
    input: lambda w: grouped_samples_input(w)
    run:
        for path in input:
            symlink_in(path, output[0])
