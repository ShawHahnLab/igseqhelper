"""
pRESTO rules for our protocol.

I needed to use our own custom demultiplexing first and change the defaults in
a few places to account for our weird method.  At this point individual
sequencer samples are aggregated together on a specimen and cell type basis, so
that duplicate sequences for the same phsyical specimen are correctly handled
as such.

See also:

 * https://github.com/ressy/example-presto
 * https://presto.readthedocs.io/en/stable/workflows/Greiff2014_Workflow.html
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from igseqhelper.data import amplicon_files, load_sequences

PRESTO_OPTS = {
    "assembly": {
        # The format of the sequence identifier which defines
        # shared coordinate information across paired ends.
        # (default: presto)
        "coord": "illumina",
        # Specify which read to reverse complement before
        # stitching. (default: tail)
        "rc": "tail",
        # Minimum sequence length to scan for overlap in de novo assembly.
        # (default: 8)
        # This sounds a bit low for our case, but when I manually check some of
        # these (in the 15 - 25 nt overlap range) they do look correct, so I'll
        # leave it alone for now.
        "minlen": 8,
        # Maximum sequence length to scan for overlap in de novo assembly.
        # (default: 1000)
        # The shortest expected sequence should dictate this.  With 2x309 we
        # could never even consider more than 309 overlap, and we don't expect
        # anything (heavy or light) shorter than around 300 nt anyway.
        # Also leaving this at the default for the moment.
        "maxlen": 1000
    },
    "qc": {
        "mean_qual": 30,
        "fwd_start": 0,
        "fwd_mode": "cut",
        "fwd_pf": "VPRIMER",
        "rev_start": 0,
        "rev_mode": "cut",
        "rev_pf": "CPRIMER"
    },
    "collapse": {
        "uf": "CPRIMER",
        "cf": "VPRIMER",
        "f": "DUPCOUNT",
        "num": 2
    }
}

def prep_primers_fwd(fp_csv_in, fp_fwd_out, custom=None, seqid="5PIIA"):
    """Take our primer CSV and create fwd FASTA for pRESTO."""
    # Actually we don't care which is which since we're doing this at the
    # specimen level
    if custom:
        custom = set(custom.values())
    else:
        custom = {}
    fwd = {"FwdPrimer"+str(idx+1): seq for idx, seq in enumerate(custom)}
    if not all(fwd.values()):
        fwd = {k: v for k, v in fwd.items() if v}
        sequences = load_sequences(fp_csv_in)
        extra = {seqid: sequences[seqid]["Seq"]}
        fwd.update(extra)
    with open(fp_fwd_out, "wt") as f_out:
        for primerid, seq in fwd.items():
            rec = SeqRecord(Seq(seq), id=primerid, description="")
            SeqIO.write(rec, f_out, "fasta")

TARGET_PRESTO_DATA = expand(amplicon_files("analysis/presto/data/{chain}.{chain_type}/{specimen}.{{rp}}.fastq", SAMPLES, "IgG+"), rp=["R1", "R2"])
TARGET_PRESTO_ASSEMBLY = amplicon_files("analysis/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq", SAMPLES, "IgG+")
TARGET_PRESTO_QC = amplicon_files("analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq", SAMPLES, "IgG+")
TARGET_PRESTO_ALL = amplicon_files("analysis/presto/collapse/{chain}.{chain_type}/{specimen}_atleast-2.fastq", SAMPLES, "IgG+")

rule all_presto_data:
    input: TARGET_PRESTO_DATA

rule all_presto_assembly:
    input: TARGET_PRESTO_ASSEMBLY

rule all_presto_qc:
    input: TARGET_PRESTO_QC

rule all_presto:
    input: TARGET_PRESTO_ALL

### Before pRESTO: demultiplex, trim reads. now, gather input data

rule presto_data:
    """Aggregate trimmed sequencer samples into per-specimen FASTQ files.

    pRESTO uses plaintext FASTQ so we'll leave them uncompressed, too.
    """
    output: "analysis/presto/data/{chain}.{chain_type}/{specimen}.{rp}.fastq"
    input: "analysis/reads-by-specimen/{chain}.{chain_type}/{specimen}.{rp}.fastq.gz"
    shell: "zcat {input} > {output}"

rule presto_primers:
    output: fwd="analysis/presto/qual/{chain}.{chain_type}/{specimen}.vprimers.fasta"
    input: sequences=ancient("metadata/sequences.csv")
    run:
        # Get custom forward primers, if any, for samples for this specimen
        custom = {k: v["FwdPrimer"] for k, v in SAMPLES.items() if v["Specimen"] == wildcards.specimen}
        prep_primers_fwd(input.sequences, output.fwd, custom)

### Paired-end Assembly

rule presto_assembly:
    output: "analysis/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq"
    input:
        r1="analysis/presto/data/{chain}.{chain_type}/{specimen}.R1.fastq",
        r2="analysis/presto/data/{chain}.{chain_type}/{specimen}.R2.fastq"
    log: "analysis/logs/presto/assemble.{chain}.{chain_type}.{specimen}.log"
    params:
        coord=PRESTO_OPTS["assembly"]["coord"],
        rc=PRESTO_OPTS["assembly"]["rc"],
        minlen=PRESTO_OPTS["assembly"]["minlen"],
        maxlen=PRESTO_OPTS["assembly"]["maxlen"]
    threads: 8
    shell:
        """
            AssemblePairs.py align \
                --minlen {params.minlen} --maxlen {params.maxlen} \
                -1 {input.r1} -2 {input.r2} \
                --nproc {threads} -o {output} --log {log} \
                --coord {params.coord} --rc {params.rc}
        """

### Quality Control

rule presto_qual:
    output: "analysis/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq"
    input: "analysis/presto/assemble/{chain}.{chain_type}/{specimen}_assemble-pass.fastq"
    log: "analysis/logs/presto/qual.{chain}.{chain_type}.{specimen}.log"
    params:
        mean_qual=PRESTO_OPTS["qc"]["mean_qual"]
    threads: 8
    shell:
        """
            FilterSeq.py quality \
                --nproc {threads} -s {input} -o {output} --log {log} \
                -q {params.mean_qual}
                
        """

rule presto_qual_fwd:
    output: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    input:
        data="analysis/presto/qual/{chain}.{chain_type}/{specimen}_quality-pass.fastq",
        primers="analysis/presto/qual/{chain}.{chain_type}/{specimen}.vprimers.fasta"
    log: "analysis/logs/presto/qual_fwd.{chain}.{chain_type}.{specimen}.log"
    params:
        start=PRESTO_OPTS["qc"]["fwd_start"],
        mode=PRESTO_OPTS["qc"]["fwd_mode"],
        pf=PRESTO_OPTS["qc"]["fwd_pf"]
    threads: 8
    shell:
        """
            MaskPrimers.py score \
                --nproc {threads} -s {input.data} -o {output} --log {log} \
                -p {input.primers} \
                --start {params.start} --mode {params.mode} --pf {params.pf}
        """

# Skipping this.
#
# I don't think we have any trace of the reverse primer visible in the sequence
# itself, only in the I1 read.
rule presto_qual_rev:
    output: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-REV_primers-pass.fastq"
    input:
        data="analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq",
        primers=""
    log: "analysis/logs/presto/qual_rev.{chain}.{chain_type}.{specimen}.log"
    params:
        start=PRESTO_OPTS["qc"]["rev_start"],
        mode=PRESTO_OPTS["qc"]["rev_mode"],
        pf=PRESTO_OPTS["qc"]["rev_pf"]
    threads: 8
    shell:
        """
            MaskPrimers.py score \
                    --nproc {threads} -s {input.data} -o {output} --log {log} \
                    -p {input.primers} \
                    --start {params.start} --mode {params.mode} --pf {params.pf} --revpr
        """

### Deduplication and Filtering

rule presto_collapse:
    output: "analysis/presto/collapse/{chain}.{chain_type}/{specimen}_collapse-unique.fastq"
    input: "analysis/presto/qual/{chain}.{chain_type}/{specimen}-FWD_primers-pass.fastq"
    log: "analysis/logs/presto/collapse.{chain}.{chain_type}.{specimen}.log"
    params:
        uf=PRESTO_OPTS["collapse"]["uf"],
        cf=PRESTO_OPTS["collapse"]["cf"]
    shell:
        """
            CollapseSeq.py \
                -s {input} -o {output} --log {log} \
                -n 20 --inner --uf {params.uf} \
                --cf {params.cf} --act set \
        """

rule presto_collapse_atleast:
    output: "analysis/presto/collapse/{chain}.{chain_type}/{specimen}_atleast-2.fastq"
    input: "analysis/presto/collapse/{chain}.{chain_type}/{specimen}_collapse-unique.fastq"
    params:
        f=PRESTO_OPTS["collapse"]["f"],
        num=PRESTO_OPTS["collapse"]["num"]
    shell:
        """
            SplitSeq.py group \
                -s {input} --outname {wildcards.specimen} \
                -f {params.f} --num {params.num}
        """
