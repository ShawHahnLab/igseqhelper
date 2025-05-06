BASEDIR = Path(str(workflow.current_basedir))
import os
os.environ["PATH"] += f":{BASEDIR}/scripts"

wildcard_constraints:
    # project metadata
    sample="[-_A-Za-z0-9]+",
    specimen="[A-Za-z0-9]+",
    subject="[A-Za-z0-9]+",
    run="[-_A-Za-z0-9]+",
    # antibody concepts
    chain="(heavy|light)",
    celltype="igm|igg",
    chain_type="(alpha|delta|gamma|mu|epsilon|kappa|lambda)",
    locus_chain="(heavy|kappa|lambda)", # ARMADiLLO's notation
    locus="(IGH|IGK|IGL)",
    segment="(V|D|J)",
    antibody_type="(IgA|IgD|IgG|IgM|IgE)",
    antibody_isolate=r"[-_A-Za-z0-9\.]+",
    antibody_lineage=r"[-_A-Za-z0-9\.]+",
    # other flow control
    compressed="((?!xz$).)*",
    num="[0-9]+", # integer
    word="[A-Za-z]+", # just letters
    suffix=r"\.?[^/.]*", # dot-whatever, which can be empty
    ext="[^/]+", # file extension (limited to one directory)
    thing="[^/]+", # any match limited to one directory
    name="[^/]+", # ditto
    rp="(R1|R2|I1|I2)"

include: "metadata.smk"
include: "by-run.smk"
include: "by-run-locus-demux.smk"
include: "mining-d.smk"
include: "igdiscover.smk"
include: "sonar.smk"
include: "sonar_recluster.smk"
include: "partis.smk"
include: "igblast.smk"
include: "fastqc.smk"
include: "armadillo.smk"
include: "reporting.smk"
include: "reporting_counts.smk"
include: "reporting_lineages.smk"
include: "reporting_sonar.smk"
include: "reporting_germline.smk"
include: "summary.smk"
