#!/usr/bin/env python

"""
Germline divergence calculation and plots.

This tool can flexibly handle sequence data input or already-prepared tabular
data, and can output tabular germline divergence data or finished plot PDFs.

Columns for tabular data in input/output:

    sequence_id
    sequence
    timepoint    x axis in plots; can parse from seq IDs if not given directly
    group        in plots, sequences with be color-coded by these labels
    locus        in plots, panels will be separated by locus
    v_call       by IgBLAST if not given
    partial      by IgBLAST if not given: "T" for sequences starting partway into V
    divergence   by IgBLAST if not given: calculated via IgBLAST if not given directly
"""

import re
import sys
import argparse
import logging
import functools
from pathlib import Path
from csv import DictReader, DictWriter
from collections import defaultdict

from pandas import DataFrame
from numpy.random import RandomState
import plotnine as p9

from igseq import igblast
from igseq import vdj
from igseq import util
from igseq.__main__ import rewrap
from igseq.record import RecordReader

LOGGER = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARNING)

DEFAULTS = {
    "plot_height": 7,
    "plot_width": 10,
    }

def nonedict():
    return defaultdict(lambda: None)

class AggrAction(argparse.Action):
    """Aggregate arguments by optional keyword"""

    def __call__(self, parser, namespace, values, option_string=None):
        existing = getattr(namespace, self.dest)
        if not existing:
            existing = nonedict()
        for value in values:
            parts = value.split("=", 1)
            if len(parts) == 2:
                existing[parts[0]] = parts[1]
            else:
                # https://stackoverflow.com/q/25314547
                func = functools.partial(lambda x: x, x=parts[0])
                existing.default_factory = func
        setattr(namespace, self.dest, existing)

class ParseInputs(argparse.Action):
    """Gather input paths with optional keyword"""

    def __call__(self, parser, namespace, values, option_string=None):
        existing = getattr(namespace, self.dest)
        if not existing:
            existing = {}
        for value in values:
            parts = value.split("=", 1)
            # use custom key if given, otherwise value is its own key
            existing[parts[0]] = parts[1] if len(parts) == 2 else parts[0]
        setattr(namespace, self.dest, existing)

def __setup_arg_parser():
    parser = argparse.ArgumentParser(
        description=rewrap(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    arg = parser.add_argument
    arg("-Q", "--query", required=True, nargs="+", action=ParseInputs, help="query files")
    arg("-o", "--output", help="output file")
    arg("--plot-width", type=int, default=DEFAULTS["plot_width"],
        help="Width of plot if output is PDF")
    arg("--plot-height", type=int, default=DEFAULTS["plot_height"],
        help="Height of plot if output is PDF")
    igblastgrp = parser.add_argument_group(
        title="IgBLAST settings")
    igblastgrp.add_argument("-r", "--reference", nargs="+",
        help="one or more FASTA/directory/builtin names pointing to V/D/J FASTA files")
    igblastgrp.add_argument("-S", "--species",
        help="species to use (human or rhesus).  Default: infer from database if possible")
    igblastgrp.add_argument("-t", "--threads", type=int, default=1,
        help="number of threads for parallel processing (default: 1)")
    querygrp = parser.add_argument_group(
        title="query attribute arguments",
        description="These arguments control how metadata is inferred for input data. "
        "Values can be prefixed with a key in order to pair them with a matching query "
        "argument, for example -Q ngs=ngs.fasta -T ngs=12 will assign timepoint 12 to "
        "all sequences in ngs.fasta.")
    def queryarg(*args, **kwargs):
        querygrp.add_argument(*args, **kwargs, nargs="+", action=AggrAction)
    queryarg("-G", "--group", default=defaultdict(lambda: "Lineage Members"),
        help="names for groups of sequences in output (e.g. \"bulk\" or \"isolate\")")
    queryarg("-P", "--pattern", default=defaultdict(lambda: "^wk([0-9]+)-"),
        help="regular expression for parsing timepoints from sequence IDs")
    queryarg("--col-timepoint", default=defaultdict(lambda: "timepoint"),
        help="column with numeric timepoints (for tabular inputs)")
    queryarg("--col-group", default=defaultdict(lambda: "group"),
        help="column with group designators (for tabular inputs)")
    queryarg("--col-seq-id", default=nonedict(),
        help="column with sequence IDs (for tabular inputs)")
    queryarg("--col-seq", default=nonedict(),
        help="column with sequences (for tabular inputs)")
    queryarg("-T", "--timepoint", default=nonedict(),
        help="default timepoint to apply to subsequent file arguments")
    return parser

def germ_div(
        paths_in, path_out, igblast_args,
        timepoint_cols=None, group_cols=None, timepoint_defs=None, timepoint_pats=None, groups=None,
        seq_id_cols=None, seq_cols=None,
        plot_width=DEFAULTS["plot_width"], plot_height=DEFAULTS["plot_height"]):
    """Calculate and/or plot germline divergence for antibody sequences."""
    def default_dict_txt(obj):
        "Just a formatting helper for logging function args here"
        out = ", ".join([f"'{key}': '{val}'" for key, val in obj.items()])
        try:
            val = obj[None]
            out = f"{out}, default: '{val}'" if out else f"default: '{val}'"
        except KeyError:
            pass
        out = "{" + out + "}"
        return out
    LOGGER.info("given ref path(s): %s", igblast_args.get("ref"))
    LOGGER.info("given query path: %s", paths_in)
    LOGGER.info("given output: %s", path_out)
    LOGGER.info("given species: %s", igblast_args.get("species"))
    LOGGER.info("given threads: %s", igblast_args.get("threads"))
    LOGGER.info("given timepoint columns: %s", default_dict_txt(timepoint_cols))
    LOGGER.info("given group columns: %s", default_dict_txt(group_cols))
    LOGGER.info("given timepoint defaults: %s", default_dict_txt(timepoint_defs))
    LOGGER.info("given timepoint patterns: %s", default_dict_txt(timepoint_pats))
    LOGGER.info("given groups: %s", default_dict_txt(groups))
    LOGGER.info("given seq ID columns: %s", default_dict_txt(seq_id_cols))
    LOGGER.info("given seq columns: %s", default_dict_txt(seq_cols))
    output = calc_germ_div(
        paths_in, igblast_args,
        timepoint_cols, group_cols, timepoint_defs, timepoint_pats, groups, seq_id_cols, seq_cols)
    germ_div_output(output, path_out, width=plot_width, height=plot_height)

def maybe_num(obj, cls=float):
    """Cast obj as float if possible, return as-is if not"""
    try:
        return cls(obj)
    except (ValueError, TypeError):
        return obj

def calc_germ_div(
        paths_in, igblast_args, timepoint_cols=None, group_cols=None,
        timepoint_defs=None, timepoint_pats=None, groups=None,
        seq_id_cols=None, seq_cols=None, partial_threshold=15):
    """Infer (with IgBLAST if needed) germline divergence for input as list of dicts"""
    ref_paths = igblast_args.get("ref")
    species = igblast_args.get("species")
    if species or ref_paths:
        if species and not ref_paths:
            # If only species is given, default to using all available reference
            # sets for that species
            ref_paths = [igblast.fuzzy_species_match(species)]
            LOGGER.info("inferred ref path: %s", ref_paths[0])
        attrs_list = vdj.parse_vdj_paths(ref_paths)
        species_det = {attrs.get("species") for attrs in attrs_list}
        species_det = {s for s in species_det if s}
        organism = igblast.detect_organism(species_det, species)
        attrs_list_grouped = vdj.group(attrs_list)
        for key, attrs_group in attrs_list_grouped.items():
            LOGGER.info("detected %s references: %d", key, len(attrs_group))
            if len(attrs_group) == 0:
                raise util.IgSeqError(f"No references for segment {key}")

    if not isinstance(paths_in, dict):
        paths_in = {path: path for path in paths_in}
    if timepoint_cols is None:
        timepoint_cols = nonedict()
    if group_cols is None:
        group_cols = nonedict()
    if timepoint_defs is None:
        timepoint_defs = nonedict()
    if timepoint_pats is None:
        timepoint_pats = nonedict()
    if groups is None:
        groups = nonedict()
    if seq_id_cols is None:
        seq_id_cols = nonedict()
    if seq_cols is None:
        seq_cols = nonedict()

    output = {}
    do_igblast = {}
    for path_key, path in paths_in.items():
        output_chunk = {}
        # Get defaults for this input
        group_col = group_cols[path_key]
        group = groups[path_key]
        tp_col = timepoint_cols[path_key]
        tp_def = timepoint_defs[path_key]
        tp_pat = timepoint_pats[path_key]
        colmap = {}
        seq_id_col = seq_id_cols[path_key]
        seq_col = seq_cols[path_key]
        if seq_id_col:
            colmap["sequence_id"] = seq_id_col
        if seq_col:
            colmap["sequence"] = seq_col
        # Set up based on existing records
        with RecordReader(path, colmap=colmap) as reader:
            for rec in reader:
                # Infer timepoint by whatever available hint
                tp_from_col = maybe_num(rec.get(tp_col))
                tp_from_pat = None
                if tp_pat:
                    match = re.search(tp_pat, rec[reader.colmap["sequence_id"]])
                    if match:
                        try:
                            tp_from_pat = match.group(1)
                        except IndexError:
                            tp_from_pat = match.group(0)
                tpoint = maybe_num(tp_from_col or tp_from_pat or tp_def)
                # Infer group
                group_from_col = rec.get(group_col)
                if group_from_col is not None:
                    group = group_from_col
                row = {
                    "sequence_id": rec[reader.colmap["sequence_id"]],
                    "timepoint": tpoint,
                    "group": group,
                    "locus": rec.get("locus", ""),
                    "v_call": rec.get("v_call", ""),
                    "partial": rec.get("partial", ""),
                    "divergence": maybe_num(rec.get("divergence", ""))}
                output_chunk[row["sequence_id"]] = row
        if any(row["divergence"] == "" for row in output_chunk.values()):
            do_igblast[path] = True
        output.update(output_chunk)

    if any(do_igblast.values()):
        with igblast.setup_db_dir(
            [str(attrs["path"]) for attrs in attrs_list]) as (db_dir, _):
            for path in paths_in.values():
                if do_igblast.get(path):
                    with igblast.run_igblast(
                            db_dir, organism, path, threads=igblast_args.get("threads", 1),
                            fmt_in=None, colmap=colmap, extra_args=["-outfmt", "19"]) as proc:
                        reader = DictReader(proc.stdout, delimiter="\t")
                        for rec in reader:
                            partial = ""
                            if int(rec["v_germline_start"]) >= partial_threshold:
                                partial = "T"
                            key = rec["sequence_id"]
                            output[key]["locus"] = rec["locus"]
                            output[key]["v_call"] = rec["v_call"]
                            output[key]["partial"] = partial
                            output[key]["divergence"] = 100 - float(rec["v_identity"])

    output = list(output.values())
    output = sorted(
        output,
        key=lambda x: (
            x["locus"],
            x["timepoint"] or 0,
            x["group"] or "",
            x["divergence"],
            x["sequence_id"]))
    output = _smooth_v_calls(output)
    return output

# for entries with multiple v calls, substitute the most abundant call across
# all entries.
def _smooth_v_calls(rows):
    v_call_tally = defaultdict(int)
    for row in rows:
        v_calls = row["v_call"].split(",")
        for v_call in v_calls:
            v_call_tally[v_call] += 1
    rows_out = []
    for row in rows:
        row = row.copy()
        v_calls = row["v_call"].split(",")
        v_call_new = sorted([(-v_call_tally[x], x) for x in v_calls])[0][1]
        row["v_call"] = v_call_new
        rows_out.append(row)
    return rows_out

def germ_div_plot(
        rows, jitter_width=None, jitter_height=None, group_colors=None,
        note_partial_seqs=False, jitter_seed=1):
    """Render a plotnine plot object for prepared germline divergence data"""
    # sort by group size to ensure members of more rare groups are shown
    # layered on top of members of more abundant groups
    # (ggplot handles the ordering of points that way, and evidently plotnine
    # does as well https://stackoverflow.com/a/16880383)
    cts = defaultdict(int)
    for row in rows:
        cts[row["group"]] += 1
    for row in rows:
        row["_groupcount"] = cts[row["group"]]
    rows.sort(key=lambda row: -row["_groupcount"])
    tp_breaks = {round(row["timepoint"]) for row in rows}
    tp_breaks.add(0)
    tp_breaks = sorted(list(tp_breaks))
    max_div = max(row["divergence"] for row in rows)

    if jitter_width is None:
        # If not specified, default the jitter width to a fraction of the x
        # axis span
        jitter_width = max(tp_breaks) / 60.0
    if jitter_height is None:
        # similar idea for jitter height
        jitter_height = max_div / 100.0

    plot_attrs = {"x": "timepoint", "y": "divergence", "color": "group"}
    if note_partial_seqs:
        plot_attrs["shape"] = "partial"
    plt = p9.ggplot(
        DataFrame(rows),
        p9.aes(**plot_attrs)) + \
        p9.geom_jitter(
            width=jitter_width,
            height=jitter_height,
            random_state=RandomState(jitter_seed)) + \
        p9.guides(
            color=p9.guide_legend(title="Category"),
            shape=False) + \
        p9.labs(
            x = "Timepoint",
            y = "V Divergence (%)") + \
        p9.coord_cartesian(
            xlim=[0, max(tp_breaks)],
            ylim=[0, max_div]) + \
        p9.scale_x_continuous(breaks=tp_breaks) + \
        p9.theme_bw()
    # Use specific default group colors if none supplied and there are exactly
    # two (otherwise colors will be assigned automatically)
    group_names = list({row["group"] for row in rows})
    if len(group_names) == 2 and not group_colors:
        group_names = sorted(group_names, key=lambda x: ("bulk" not in x.lower(), x))
        group_colors = {group_names[0]: "#000000", group_names[1]: "#ff0000"}
    if group_colors:
        plt += p9.scale_color_manual(values=group_colors)
    # If multiple loci, facet plot by locus
    loci = {row["locus"] for row in rows}
    loci = sorted(list(loci), key=lambda locus: (locus != "IGH", locus))
    if len(loci) > 1:
        LOGGER.info("Faceting plot by locus (%s)", loci)
        plt += p9.facet_grid("locus ~ .")
    return plt

def germ_div_output(rows, path_out, **kwargs):
    """Write germline divergence info to file, either as CSV or as a plot PDF"""
    ext = Path(path_out).suffix.lower()
    if ext == ".csv" or path_out == "-":
        for row in rows:
            row["divergence"] = "" if row["divergence"] is None else f"{row['divergence']:.6f}"
        try:
            handle = sys.stdout if path_out == "-" else open(path_out, "wt", encoding="UTF8")
            writer = DictWriter(handle, fieldnames = rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        finally:
            handle.close()
    elif ext == ".pdf":
        plt = germ_div_plot(rows)
        plt.save(path_out, **kwargs)

def main(arglist=None):
    """CLI for germ_div"""
    parser = __setup_arg_parser()
    args = parser.parse_args() if arglist is None else parser.parse_args(arglist)
    igblast_args = {
        "ref": args.reference,
        "species": args.species,
        "threads": args.threads}
    germ_div(
        args.query, args.output, igblast_args,
        args.col_timepoint, args.col_group, args.timepoint, args.pattern, args.group,
        args.col_seq_id, args.col_seq, args.plot_width, args.plot_height)

if __name__ == "__main__":
    main()
