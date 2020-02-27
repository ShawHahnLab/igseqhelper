#!/usr/bin/env python

"""
Immunoglobulin Sequence Data Analysis

Various Python helpers for our analysis here, particularly to use in Snakemake
rules.  See also the igseq R package.
"""

from pathlib import Path
import logging

def __find_r_pkg():
    LOGGER.debug("finding R package")
    igseq_inst = Path(__file__).parent.parent.parent
    if igseq_inst.name == "inst":
        igseq = igseq_inst.parent
    else:
        igseq_inst = igseq_inst.parent
        igseq = igseq_inst
    igseq_inst = igseq_inst.resolve()
    igseq = igseq.resolve()
    LOGGER.debug("finding R package: inst: %s", str(igseq_inst))
    LOGGER.debug("finding R package: pkg: %s", str(igseq))
    return igseq_inst, igseq

LOGGER = logging.getLogger(__name__)

R_PKG_INST, R_PKG_PATH = __find_r_pkg()
