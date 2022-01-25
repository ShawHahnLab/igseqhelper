# demux any run that gets a mention in the samples table
# (not using the simpler RUNS.keys because some runs we have an entry for
# aren't used in the demux step)
RUNS_FOR_SAMPLES = set([attrs["Run"] for attrs in SAMPLES.values() if attrs["Run"]])

# do early steps on any run that's labeled IgSeq
RUNS_FOR_IGSEQ = [runid for runid in RUNS if "IgSeq" in RUNS[runid]["Protocol"]]

rule all_igblast_merge:
    input: expand("analysis/igblast/merge/{run}.done", run=RUNS_FOR_SAMPLES)

rule all_merge:
    input: expand("analysis/merge/{run}.done", run=RUNS_FOR_SAMPLES)

rule all_trim:
    input: expand("analysis/trim/{run}.done", run=RUNS_FOR_SAMPLES)

rule all_demux:
    input: expand("analysis/demux/{run}", run=RUNS_FOR_SAMPLES)

rule all_phix:
    input: expand("analysis/phix/{run}/phix.bam", run=RUNS_FOR_IGSEQ)

rule all_getreads:
    input: expand("analysis/reads/{run}", run=RUNS_FOR_IGSEQ)

### By-pair rules grouped by run

def run_igblastmerge_input(w):
    targets = []
    for sample, attrs in SAMPLES.items():
        if attrs["Run"] == w.run:
            targets.append(f"analysis/igblast/merge/{{run}}/{sample}.tsv.gz")
    return targets

rule run_igblast_merge:
    output: "analysis/igblast/merge/{run}.done"
    input: run_igblast_merge_input
    shell: "touch {output}"

def run_merge_input(w):
    targets = []
    suffixes = ["fastq.gz", "merge.counts.csv", "pear.log"]
    for sample, attrs in SAMPLES.items():
        if attrs["Run"] == w.run:
            for suffix in suffixes:
                targets.append(f"analysis/merge/{{run}}/{sample}.{suffix}")
    return targets

rule run_merge:
    output: "analysis/merge/{run}.done"
    input: run_merge_input
    shell: "touch {output}"

def run_trim_input(w):
    targets = []
    suffixes = ["R1.fastq.gz", "R2.fastq.gz", "cutadapt1.json", "cutadapt2.json", "trim.counts.csv"]
    for sample, attrs in SAMPLES.items():
        if attrs["Run"] == w.run:
            for suffix in suffixes:
                targets.append(f"analysis/trim/{{run}}/{sample}.{suffix}")
    return targets

rule run_trim:
    output: "analysis/trim/{run}.done"
    input: run_trim_input 
    shell: "touch {output}"

### By-pair rules

rule merged_igblast:
    output:
        tsv="analysis/igblast/merge/{run}/{sample}.tsv.gz"
    input:
        fqgz="analysis/merge/{run}/{sample}.fastq.gz"
    threads: 4
    shell:
        """
            igseq igblast -r rhesus/imgt -t {threads} -Q {input.fqgz} --outfmt 19 | gzip > {output}
        """

rule pair_merge:
    output:
        fqgz="analysis/merge/{run}/{sample}.fastq.gz",
        counts="analysis/merge/{run}/{sample}.merge.counts.csv",
        log="analysis/merge/{run}/{sample}.pear.log"
    input:
        r1="analysis/trim/{run}/{sample}.R1.fastq.gz",
        r2="analysis/trim/{run}/{sample}.R2.fastq.gz"
    threads: 4
    shell:
        """
            igseq merge -t {threads} {input.r1} {input.r2}
        """

def pair_trim_input(w):
    for samp, attrs in SAMPLES.items():
        if samp == w.sample:
            runid = attrs["Run"]
            if runid == w.run:
                return {
                "demux": "analysis/demux/{run}",
                "samples": ancient("metadata/samples.csv")}
            else:
                raise ValueError(f"Sample {w.sample} has run {runid}, not {w.run}")
    raise ValueError(f"Sample {w.sample} not found for run {w.run}")

rule pair_trim:
    output:
        r1="analysis/trim/{run}/{sample}.R1.fastq.gz",
        r2="analysis/trim/{run}/{sample}.R2.fastq.gz",
        report1="analysis/trim/{run}/{sample}.cutadapt1.json",
        report2="analysis/trim/{run}/{sample}.cutadapt2.json",
        counts="analysis/trim/{run}/{sample}.trim.counts.csv"
    input: unpack(pair_trim_input)
    threads: 4
    shell:
        """
            igseq trim -t {threads} \
                --samples {input.samples} \
                --species rhesus \
                {input.demux}/{wildcards.sample}.R1.fastq.gz \
                {input.demux}/{wildcards.sample}.R2.fastq.gz
        """

rule pair_phix:
    output:
        bam="analysis/phix/{run}/phix.bam",
        counts="analysis/phix/{run}/phix.counts.csv"
    input:
        demux="analysis/demux/{run}"
    threads: 4
    shell:
        """
            igseq phix -t {threads} {input}
        """

### By-run rules

rule demux:
    output: directory("analysis/demux/{run}")
    input:
        reads="analysis/reads/{run}",
        samples=ancient("metadata/samples.csv")
    shell: "igseq demux --samples {input.samples} --details {output}/details.csv.gz {input.reads}"

rule getreads:
    output: directory("analysis/reads/{run}")
    input: ancient("/seq/runs/{run}")
    threads: 28
    shell: "igseq getreads -t {threads} --threads-load $(({threads}<4 ? {threads} : 4)) {input}"

### Samples grouped by subject or specimen
# Each directory will contain a set of symbolic links pointing to the
# individual merged sample files

def grouped_samples_input(w, pattern="analysis/merge/{runid}/{samp}.fastq.gz"):
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
            if attrs["SpecimenAttrs"]["Subject"] == w.name and attrs["Type"] == w.type:
                if not celltype:
                    raise ValueError("cell type needed if grouping by subject")
                if cellmatch(w.celltype, attrs["SpecimenAttrs"]["CellType"]):
                    runid = attrs["Run"]
                    targets.extend(expand(pattern, runid=runid, samp=samp))
    elif w.thing == "specimen":
        for samp, attrs in SAMPLES.items():
            if attrs["Specimen"] == w.name and attrs["Type"] == w.type:
                if not celltype or cellmatch(w.celltype, attrs["SpecimenAttrs"]["CellType"]):
                    runid = attrs["Run"]
                    targets.extend(expand(pattern, runid=runid, samp=samp))
    else:
        raise ValueError("unrecognized wildcards for rule grouped_samples_input")
    if not targets:
        raise ValueError("no matching inputs for rule grouped_samples_input")
    return targets

def symlink_in(path_real, link_dir):
    # careful, symlink output/name TO each path!
    # ideally I'd do relative links, but pathlib's relative_to
    # apparently can't handle things like
    # Path("foo").relative_to("foo/bar/baz") -> "../.."
    link_from = (Path(link_dir)/Path(path_real).name).resolve()
    link_from.parent.mkdir(exist_ok=True, parents=True)
    link_to = Path(path_real).resolve()
    link_from.symlink_to(link_to)

rule grouped_samples_by_celltype:
    output: directory("analysis/samples-by-{thing}/{celltype}/{name}.{type}")
    input: lambda w: grouped_samples_input(w)
    run:
        for path in input:
            symlink_in(path, output[0])

rule grouped_samples:
    output: directory("analysis/samples-by-{thing}/{name}.{type}")
    input: lambda w: grouped_samples_input(w)
    run:
        for path in input:
            symlink_in(path, output[0])
