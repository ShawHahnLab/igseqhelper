"""
Processing steps from raw read data to finalized antibody results.
"""

BASEDIR = Path(str(workflow.current_basedir))
import os
os.environ["PATH"] += f":{BASEDIR}/scripts"

wildcard_constraints:
    # project metadata
    sample="[-A-Za-z0-9]+",
    specimen="[A-Za-z0-9]+",
    subject="[A-Za-z0-9]+",
    run="[-_A-Za-z0-9]+",
    # antibody concepts
    chain="(heavy|light)",
    celltype="igm|igg",
    chain_type="(alpha|delta|gamma|mu|epsilon|kappa|lambda)",
    locus="(IGH|IGK|IGL)",
    segment="(V|D|J)",
    antibody_type="(IgA|IgD|IgG|IgM|IgE)",
    antibody_isolate=r"[-_A-Za-z0-9\.]+",
    antibody_lineage=r"[-_A-Za-z0-9\.]+",
    # other flow control
    compressed="((?!xz$).)*",
    chunk="[0-9]+",
    rp="(R1|R2|I1|I2)"

include: "workflow/metadata.smk"
include: "workflow/by-run.smk"
include: "workflow/mining-d.smk"
include: "workflow/igdiscover.smk"
include: "workflow/sonar.smk"
include: "workflow/igblast.smk"
include: "workflow/fastqc.smk"
include: "workflow/reporting.smk"
