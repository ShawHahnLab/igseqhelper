"""
Generic IgBLAST rules.

To igblast a path from the analysis directory to one of the built-in rhesus
references in igseq, ask for:

analysis/igblast/{ref}/{path}.tsv.gz

...which will use analysis/{path} as input and treat {ref} as an argument to
igseq igblast -r.

To igblast a path to a custom VDJ database on disk:

analysis/igblast/custom-{ref}/{path}.tsv.gz

... will look for analysis/igblast/ref-custom-{ref} as an input and the
argument to igblast -r.

custom-{subject}.{locus} will automatically refer to the germline reference for
that subject and locus.
"""

# No input reference required unless prefixed with "custom-"
def input_for_igblast(w):
    # special handling for .fasta: use xz-compressed FASTA directly if that's
    # available on disk and the expected .fa isn't, but otherwise just use the
    # input path that's implied
    query = Path(f"analysis/{w.path}")
    if not query.exists() and query.suffix == ".fa" and Path(str(query) + ".xz").exists():
        query = str(query) + ".xz"
    targets = {"query": query}
    if w.ref.startswith("custom-"):
        targets["ref"] = f"analysis/igblast/ref-{w.ref}"
    return targets

# If there's an input named ref, use that; othewise assume it's a built-in
# rhesus reference from igseq
def ref_for_igblast(w, input):
    try:
        return input.ref
    except AttributeError:
        return f"rhesus/{w.ref}"

rule igblast:
    output: "analysis/igblast/{ref,[^/]+}/{path}.tsv.gz"
    input: unpack(input_for_igblast)
    params:
        ref=ref_for_igblast,
        outfmt=19, # AIRR TSV format
        species="rhesus"
    log:
        conda="analysis/igblast/{ref}/{path}.conda_build.txt"
    conda: "envs/igseq.yaml"
    threads: 8
    # include special handling for xz-compressed FASTA I have for some older
    # output files.  Also handle empty files (igseq's igblast wrapper, as of
    # 0.7.0rc1, will hang if given empty input.)
    shell:
        """
            conda list --explicit > {log.conda}
            if [[ {input.query} =~ \.fa\.xz$ ]]; then
              xzcat {input.query} | igseq igblast --input-format fa -Q - \
                  -t {threads} -S {params.species} -r {params.ref} \
                  -outfmt {params.outfmt} | gzip > {output}
            else
              if [[ -s {input.query} ]]; then
                igseq igblast -Q {input.query} \
                    -t {threads} -S {params.species} -r {params.ref} \
                    -outfmt {params.outfmt} | gzip > {output}
              else
                touch {output}
              fi
            fi
        """

# One automatically-available custom reference: per-subject per-locus germline
# files
rule igblast_ref_germline:
    output: directory("analysis/igblast/ref-custom-{subject}.{locus}")
    input:
        fastas=expand("analysis/germline/{{subject}}.{{locus}}/{segment}.fasta", segment=["V", "D", "J"])
    params:
        germline_dir="analysis/germline/{subject}.{locus}"
    run:
        link_target = relative_link_target(params.germline_dir, output[0])
        Path(output[0]).symlink_to(link_target)

def igblast_setup_sonar_helper_rules():
    """Define per-subject SONAR igblast rules dynamically"""
    # for each subject, take all specimens that include at least one sample for
    # the applicable locus.  note that this spans all cell types so it'll
    # include IgM+/IgD+ too.
    groups = defaultdict(list)
    for samp_attrs in SAMPLES.values():
        if samp_attrs["Skip"]:
            continue
        spec_attrs = samp_attrs.get("SpecimenAttrs")
        if not spec_attrs:
            continue
        subject = spec_attrs["Subject"]
        locus = {"kappa": "IGK", "lambda": "IGL"}.get(samp_attrs["Type"], "IGH")
        groups[(subject, locus)].append((samp_attrs["Type"], samp_attrs["Specimen"]))
    for key, pairs in groups.items():
        subject, locus = key
        chain_types = [pair[0] for pair in pairs]
        specimens = [pair[1] for pair in pairs]
        targets = expand(
            (f"analysis/igblast/custom-{subject}.{locus}/"
            f"sonar/{subject}.{{chain_type}}/{{specimen}}/"
            f"output/sequences/nucleotide/{{specimen}}_goodVJ_unique.fa.tsv.gz"),
            zip, chain_type=chain_types, specimen=specimens)
        rule:
            f"IgBLAST SONAR goodVJ unique FASTAs for one subject+locus using personal germline reference"
            name: f"igblast_sonar_{subject}_{locus}"
            input: targets

igblast_setup_sonar_helper_rules()
