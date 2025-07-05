import re

PROFILES = {
    "selective_chloroperoxidases": "/home/szenei/hal_miner/halogenase_detection/phmm_db/VHPO_VCPO_selective.hmm",
    "non-selective_chloroperoxidases": "/home/szenei/hal_miner/halogenase_detection/phmm_db/VHPO_VCPO_nonselective.hmm",
    "bromoperoxidases": "/home/szenei/hal_miner/halogenase_detection/phmm_db/VHPO_VBPO.hmm",
    "iodoperoxidases": "halogenase_detection/phmm_db/VHPO_VIPO.hmm"
}

# DOI: 10.1039/D0CS01551B
FDHs = {
    "flavin-binding": {"signature": "G.G..G",
                       "region": list(range(0, 100))
                       },
    "no_monoox": {"signature": "W.W.I.",
                  "region": list(range(206, 212))
                  },
    "tunnel_line": {"signature": "F.*P.*S.G",
                    "region": list(range(280, 305))
                    }
}


# doi/full/10.1021/acscatal.2c01184
SAM = {
    "c_terminal_motig": {"signature": "(RNA{2}|RNG{2}|Y{2}[GA]{2})",
                         "region": list(range(200, 290))
                         },
}

# We can differentiate NHFeHals from NHFe-dependent hydroxylases
# and dioxygenases by looking at the catalytic triad.
# HXD/E and HXG/A are quite small motifs that can occur more
# in the alignment. Searching these motifs as regex will probably return several hits, so
# worth looking at the alignment.
NONHEME_IRON = {
    "halogenase": {"signature": "H.[GA]",
                   "motif": list(range(150, 350))
                   },
    "non_halogenase": {"signature": "H.[DE]",
                       "motif": list(range(150, 350))
                       }
}

# DIMETAL-CARBOXYLATE
# doi/10.1021/acs.biochem.4c00720
DIMETAL = {
    "first_motif": {"signature": "E..QE..H",
                    "region": list(range(120, 180))
                    },
    "second_motif": {"signature": "H..DE..H",
                     "region": list(range(300, 350))
                     }
}

# VCPO
VCPO = {
    "selectivity_residues": {"signature": "KS",
                           "region": [191, 295]
                           }
}

# VBPO
# doi.org/10.1016/j.bioorg.2012.05.003
VBPO = {
    "first_active_site": {"signature": "K...H...RPEA",
                          "region": list(range(300, 350))
                          },
    "second_active_site": {"signature": "HP.Y..GHA",
                           "region": list(range(355, 425))
                           },
    "intermolecular_bridges": {"first_motif": "C...D",
                               "second_motif": "C....D...C",
                               "region": list(range(100, 425))
                               },
    "intramolecular_bridges":  {"first_signature": "C.P.P",
                               "second_signature": "I.*C.*LT.EGE.NK",
                               "third_signature": "C.G..TG...C",
                               "region": list(range(100, 425))
                               },
}

VIPO = {
    "catalytic_residues": {"signature": "YRFHRH",
                           "region": [261, 329, 351, 358, 408, 414]
                           }
}