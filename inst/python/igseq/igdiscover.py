"""
Helpers for running IgDiscover.
"""

# Mapping of chain heavy/light names to heavy/light short specifiers
CHAINS = {"heavy": "H", "light": "L"}
# Mapping of chain types to heavy/light short specifiers per locus.  Note the
# potential confusion between L here and L above though.
CHAIN_TYPES = {
    "alpha": "H",
    "delta": "H",
    "gamma": "H",
    "mu": "H",
    "epsilon": "H",
    "kappa": "K",
    "lambda": "L"}
# Mapping of chain heavy/light names to list of segments included
SEGMENTS = {"heavy": ["V", "D", "J"], "light": ["V", "J"]}
