"""
Processing steps from raw read data to finalized antibody results.
"""

BASEDIR = Path(str(workflow.current_basedir))
import os
import sys
sys.path.append(str(BASEDIR/"python"))
os.environ["PATH"] += f":{BASEDIR}/scripts"

wildcard_constraints:
    # project metadata
    sample="[A-Za-z0-9]+",
    specimen="[A-Za-z0-9]+",
    subject="[A-Za-z0-9]+",
    run="[-_A-Za-z0-9]+",
    # antibody concepts
    chain="(heavy|light)",
    celltype="igm|igg",
    chain_type="(alpha|delta|gamma|mu|epsilon|kappa|lambda)",
    segment="(V|D|J)",
    antibody_type="(IgA|IgD|IgG|IgM|IgE)",
    antibody_isolate=r"[-_A-Za-z0-9\.]+",
    antibody_lineage=r"[-_A-Za-z0-9\.]+",
    # other flow control
    compressed="((?!xz$).)*",
    chunk="[0-9]+",
    rp="(R1|R2|I1|I2)"

include: "snakemake-rules/metadata.smk"
include: "snakemake-rules/by-run.smk"
include: "snakemake-rules/mining-d.smk"
include: "snakemake-rules/igdiscover.smk"
include: "snakemake-rules/sonar.smk"
include: "snakemake-rules/igblast.smk"
include: "snakemake-rules/fastqc.smk"
include: "snakemake-rules/reporting.smk"
