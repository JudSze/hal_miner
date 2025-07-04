import re

# DOI: 10.1039/D0CS01551B
GENERAL_FDH_MOTIFS = {
    "G.G..G": list(range(0, 100)),      # flavin-binding
    "W.W.I.": list(range(206, 212)),    # prevents monooxygenation activity
    "F.*P.*S.G": list(range(280, 305))  # lines the tunnel between the cofactor and the site of halogenaiton
}


# doi/full/10.1021/acscatal.2c01184
C_TERMINAL_MOTIF = {
    "(RNA{2}|RNG{2}|Y{2}[GA]{2})": list(range(200, 290)),
}

# We can differentiate NHFeHals from NHFe-dependent hydroxylases
# and dioxygenases by looking at the catalytic triad.
# HXD/E and HXG/A are quite small motifs that can occur more
# in the alignment. Searching these motifs as regex will probably return several hits, so
# worth looking at the alignment.
CATALYTIC_TRIAD = {
    "H.[GA]": list(range(150, 350)),
    "H.[DE]": list(range(150, 350))
}

# DIMETAL-CARBOXYLATE
# doi/10.1021/acs.biochem.4c00720
DIMETAL_MOTIFS = {
    "E..QE..H": list(range(120, 180)),
    "H..DE..H": list(range(300, 350))
}

# VCPO
SELECTIVITY_RESIDUES = {
    "KS": [191, 295]
}

# VBPO
# doi.org/10.1016/j.bioorg.2012.05.003
VBPO_ACTIVE_SITE_CONTAINING = {
    "HP.Y.GHA": list(range(500, 545)),
    "KW..H...RPEA": list(range(410, 455))
}

VBPO_INTERMOLECULAR_BRIDGE = {
    "C...D"
    "C....D...C"
}

# these might be taxonomy-specific
VBPO_INTRAMOLECULAR_BRIDGE = {
    "C.P.P"
    "I.*C.*LT.EGE.NK"
    "C.G..TG...C"
}
