import re
from halogenase_miner.categorization.signature_search import (
    get_catalytic_residues,
    search_motif
)
from halogenase_miner.categorization.family_categorization import EnzymeFamily
from halogenase_miner.motif_db.motifs import VBPO, COPPER

res = EnzymeFamily("/home/szenei/hal_miner/copper-dependent_hits.fasta")
test_coppers=res.copper_dependent.copper_binding_motifs(mode="strict")
test_coppers

first_motif = search_motif(res.copper_dependent.hits, COPPER, "first_motif")
second_motif = search_motif(res.copper_dependent.hits, COPPER, "second_motif")

both_motifs = list(set(first_motif) & set(second_motif))
len(both_motifs)


with open("copper-dependent_both_motifs.txt", "w") as f:
    for accession in both_motifs:
        f.write(f"{accession}\n")

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

VBPO["intramolecular_bridges"].keys()
