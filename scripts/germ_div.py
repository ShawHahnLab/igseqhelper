#!/usr/bin/env python

"""
Germline divergence calculation and plots"""

# TODO allow toggle between boxplot and jitter versions

import re
import sys
import argparse
import logging
import functools
from pathlib import Path
from csv import DictReader, DictWriter
from collections import defaultdict

from pandas import DataFrame
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

    parser.add_argument("-r", "--reference", nargs="+",
        help="one or more FASTA/directory/builtin names pointing to V/D/J FASTA files")
    parser.add_argument("-S", "--species",
        help="species to use (human or rhesus).  Default: infer from database if possible")
    parser.add_argument("-Q", "--query", nargs="+", action=ParseInputs, help="query files")
    parser.add_argument("-o", "--output", help="output file")
    parser.add_argument(
        "-G", "--group", nargs="+", action=AggrAction,
        default=defaultdict(lambda: "Lineage Members"),
        help="names for groups of sequences in output")
    parser.add_argument(
        "-P", "--pattern", nargs="+", action=AggrAction,
        default=defaultdict(lambda: "^wk([0-9]+)-"),
        help="regular expression for parsing timepoints from sequence IDs")
    parser.add_argument(
        "--col-timepoint", nargs="+", action=AggrAction,
        default=defaultdict(lambda: "timepoint"),
        help="column with numeric timepoints (for tabular inputs)")
    parser.add_argument(
        "--col-group", nargs="+", action=AggrAction,
        default=defaultdict(lambda: "group"),
        help="column with group designators (for tabular inputs)")
    parser.add_argument(
        "--col-seq-id", nargs="+", action=AggrAction, default=nonedict(),
        help="column with sequence IDs (for tabular inputs)")
    parser.add_argument(
        "--col-seq", nargs="+", action=AggrAction, default=nonedict(),
        help="column with sequences (for tabular inputs)")
    parser.add_argument(
        "-T", "--timepoint", nargs="+", action=AggrAction, default=nonedict(),
        help="default timepoint to apply to subsequent file arguments")
    parser.add_argument(
        "--plot-width", default=DEFAULTS["plot_width"], help="Width of plot if output is PDF")
    parser.add_argument(
        "--plot-height", default=DEFAULTS["plot_height"], help="Height of plot if output is PDF")
    return parser

def germ_div(
        ref_paths, paths_in, path_out, species=None,
        timepoint_cols=None, group_cols=None, timepoint_defs=None, timepoint_pats=None, groups=None,
        seq_id_cols=None, seq_cols=None,
        plot_width=DEFAULTS["plot_width"], plot_height=DEFAULTS["plot_height"]):
    def default_dict_txt(obj):
        out = ", ".join([f"'{key}': '{val}'" for key, val in obj.items()])
        try:
            val = obj[None]
            if out:
                out = f"{out}, default: '{val}'"
            else:
                out = f"default: '{val}'"
        except KeyError:
            pass
        out = "{" + out + "}"
        return out
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given query path: %s", paths_in)
    LOGGER.info("given output: %s", path_out)
    LOGGER.info("given species: %s", species)
    LOGGER.info("given timepoint columns: %s", default_dict_txt(timepoint_cols))
    LOGGER.info("given group columns: %s", default_dict_txt(group_cols))
    LOGGER.info("given timepoint defaults: %s", default_dict_txt(timepoint_defs))
    LOGGER.info("given timepoint patterns: %s", default_dict_txt(timepoint_pats))
    LOGGER.info("given groups: %s", default_dict_txt(groups))
    LOGGER.info("given seq ID columns: %s", default_dict_txt(seq_id_cols))
    LOGGER.info("given seq columns: %s", default_dict_txt(seq_cols))
    output = calc_germ_div(
        ref_paths, paths_in, species,
        timepoint_cols, group_cols, timepoint_defs, timepoint_pats, groups, seq_id_cols, seq_cols)
    germ_div_output(output, path_out, width=plot_width, height=plot_height)

def maybe_num(obj, cls=float):
    try:
        return cls(obj)
    except (ValueError, TypeError):
        return obj

def calc_germ_div(ref_paths, paths_in, species=None, timepoint_cols=None, group_cols=None, timepoint_defs=None, timepoint_pats=None, groups=None, seq_id_cols=None, seq_cols=None, partial_threshold=15):
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
                if tp_from_col is not None:
                    tp = tp_from_col
                elif tp_from_pat is not None:
                    tp = tp_from_pat
                else:
                    tp = tp_def
                if tp is not None:
                    tp = float(tp)
                # Infer group
                group_from_col = rec.get(group_col)
                if group_from_col is not None:
                    group = group_from_col
                row = {
                    "sequence_id": rec[reader.colmap["sequence_id"]],
                    "timepoint": tp,
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
                    with igblast.run_igblast(db_dir, organism, path, threads=1, fmt_in=None, colmap=colmap, extra_args=["-outfmt", "19"]) as proc:
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
        key=lambda x: [x["locus"], x["timepoint"] or 0, x["group"] or "", x["divergence"], x["sequence_id"]])
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

def germ_div_plot(rows, jitter_width=None, jitter_height=None):
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

    plt = p9.ggplot(
        DataFrame(rows),
        p9.aes(x = "timepoint", y = "divergence", color = "group", shape = "partial")) + \
        p9.geom_jitter(width=jitter_width, height=jitter_height) + \
        p9.labs(
            x = "Timepoint",
            y = "V Divergence (%)") + \
        p9.coord_cartesian(
            xlim=[0, max(tp_breaks)],
            ylim=[0, max_div]) + \
        p9.scale_x_continuous(breaks=tp_breaks) + \
        p9.theme_bw()
    loci = {row["locus"] for row in rows}
    if len(loci) > 1:
        LOGGER.info("Faceting plot by locus (%s)", loci)
        plt += p9.facet_grid("locus ~ .")
    return plt

def germ_div_output(rows, path_out, **kwargs):
    ext = Path(path_out).suffix.lower()
    if ext == ".csv":
        for row in rows:
            if row["divergence"] is None:
                row["divergence"] = ""
            else:
                row["divergence"] = f"{row['divergence']:.6f}"
        with open(path_out, "wt") as f_out:
            writer = DictWriter(f_out, fieldnames = rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    elif ext == ".pdf":
        plt = germ_div_plot(rows)
        plt.save(path_out, **kwargs)

def main(arglist=None):
    parser = __setup_arg_parser()
    if arglist is None:
        args, args_extra = parser.parse_known_args()
    else:
        args, args_extra = parser.parse_known_args(arglist)
    if not vars(args):
        parser.print_help()
        sys.exit(0)
    germ_div(
        args.reference, args.query, args.output, args.species,
        args.col_timepoint, args.col_group, args.timepoint, args.pattern, args.group,
        args.col_seq_id, args.col_seq, args.plot_width, args.plot_height)

if __name__ == "__main__":
    main()
