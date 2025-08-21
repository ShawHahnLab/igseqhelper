"""
Rules for reporting on the individualized germline references derived from IgDiscover.
"""

rule report_germline_summary:
    """Wide-format per-subject summary table of all individualized germline details"""
    output: "analysis/reporting/germline/summary.csv"
    shell: "report_germline_summary.py analysis/reporting/germline {output}"

rule report_germline:
    """Per-subject per-segment germline details"""
    output: "analysis/reporting/germline/{subject}.{locus}/{segment}.csv"
    input:
        fasta="analysis/germline/{subject}.{locus}/{segment}.fasta",
        info="analysis/germline/{subject}.{locus}/info.csv"
    shell:
        """
            report_germline.py {wildcards.subject} {wildcards.locus} {wildcards.segment} \
                {input.fasta} {input.info} {output}
        """
