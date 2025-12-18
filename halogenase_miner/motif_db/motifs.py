# do predictive regexes, strict or loose mode
import pkg_resources

class Profiles:
    fdh_all_conventional = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/FDH_all_conventional.hmm')
    fdh_unconventional = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/FDH_unconventional.hmm')
    fdh_cycline_orsellinic = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/FDH_cycline_orsellinic.hmm')
    fdh_pyrrole = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/FDH_pyrrole.hmm')
    fdh_trp_5 = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/FDH_trp_5.hmm')
    fdh_trp_6_7 = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/FDH_trp_6_7.hmm')
    fdh_tyrosine = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/FDH_tyrosine-like_hpg.hmm')
    nhfe_indole_alkaloid = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/NHFe_indole_alkaloid.hmm')
    nhfe_nucleotide = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/NHFe_nucleotide.hmm')
    nhfe_small_amino_acid = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/NHFe_small_amino_acids.hmm')
    nhfe_variant_B = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/NHFe_variant_B.hmm')
    sam_chlorinases = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/SAM_chlorinase.hmm')
    sam_fluorinases = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/SAM_fluorinase.hmm')
    vhpo_general = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/VHPO_general.hmm')
    vhpo_selective_chloroperoxidases = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/VHPO_VCPO_selective.hmm')
    vhpo_non_selective_chloroperoxidases = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/VHPO_VCPO_nonselective.hmm')
    vhpo_bromoperoxidases = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/VHPO_VBPO.hmm')
    vhpo_iodoperoxidases = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/VHPO_VIPO.hmm')
    dimetal_carboxylate = pkg_resources.resource_filename('halogenase_miner', 'phmm_db/dimetal-carboxylate.hmm')
    copper_dependent = pkg_resources.resource_filename('halogenase_miner', "phmm_db/copper-dependent.hmm")

# DOI: 10.1039/D0CS01551B
FDHs = {
    "flavin_binding": {"signature": "G.G..G",
                       "region": list(range(0, 100))
                       },
    "no_monoox": {"signature": "W.W.I",
                  "region": list(range(206, 212))
                  },
    "tunnel_line": {"signature": "F.*P.*S.*G",
                    "region": list(range(280, 350))
                    },
    "unvonv_flavin_binding": {"signature": "G.G..A",
                              "region": list(range(0, 100))},
    "unconv_broad_motif": {"signature": "F.*P.*S.*G",
                           "region": list(range(10, 500))}
}


# doi/full/10.1021/acscatal.2c01184
SAM = {
    "c_terminal_motif": {"signature": "(RNA{2}|RNG{2}|Y{2}[GA]{2})",
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
                   "region": list(range(100, 180))
                   },
    "non_halogenase": {"signature": "H.[DE]",
                       "region": list(range(100, 150))
                       },
    "nucleotide": {"signature": "H.[GA]",
                   "region": list(range(15, 50))},
    "unique_variant_b": {"signature": "H.S.P",
                         "region": list(range(200, 250))}
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
# active site and vanadium binding site should be noted
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
    "intramolecular_bridges":  {"first_motif": "C.P.P", # look for structures to figure if they are shifted among the nonmatches
                               "second_motif": "I.*C.*LT.EGE.NK", # hmm visualization tool for retrieving the pHMM positions, include statistics
                               "third_motif": "C.G..TG...C",
                               "region": list(range(100, 425))
                               },
}

VIPO = {
    "catalytic_residues": {"signature": "YRFHRH",
                           "region": [261, 329, 351, 358, 408, 414]
                           }
}

# Copper-dependent halogenase
# doi.org/10.1038/s41586-024-08362-4
# two HxxHC motifs, putative metal binding sites
# regioselective iodination, thiocyanation, selenocyanation
COPPER = { "first_motif": {"signature": "H..HC",
                           "region": list(range(65, 80))
                           },
            "second_motif": {"signature": "H..HC",
                             "region": list(range(95, 110))}
}