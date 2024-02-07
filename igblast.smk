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
    targets = {"query": f"analysis/{w.path}"}
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
    threads: 8
    shell:
        """
            igseq igblast \
                -t {threads} \
                -S {params.species} \
                -r {params.ref} \
                -Q {input.query} \
                -outfmt {params.outfmt} | gzip > {output}
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
