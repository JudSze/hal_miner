import re
from halogenase_miner.categorization.signature_search import (
    get_catalytic_residues,
    search_motif
)
from halogenase_miner.categorization.family_categorization import EnzymeFamily
from halogenase_miner.motif_db.motifs import VBPO

res = EnzymeFamily("test/fdh_unconventional.fasta")

# Get overview based on the active sites and other catalytic residues
res.vanadium_dependent.vhpo_categorize_selective_chlorinases()
res.vanadium_dependent.vhpo_categorize_iodinases()
res.vanadium_dependent.vhpo_categorize_brominases()
res.vanadium_dependent.vhpo_categorize_selective_chlorinases()
res.vanadium_dependent.vhpo_intermol_categorize_brominases()

# Chech for sequences containing the Cs forming intramolecular bridges in VBPOs
intramolecular_matches = []
for protein in res.vanadium_dependent.brominase_hits:
    if re.search(VBPO["intramolecular_bridges"]["first_motif"]):
        intramolecular_matches = (res.vanadium_dependent.brominase_hits, VBPO, "intramolecular_bridges")
intramolecular_matches

l = dict.fromkeys(VBPO["intramolecular_bridges"])
l.pop('region')
l

VBPO["intramolecular_bridges"].keys()

from halogenase_miner.motif_db import specific