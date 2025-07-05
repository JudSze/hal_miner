import re

from motif_db.motifs import (
    VBPO,
    VCPO,
    VIPO
)
from signature_search import (
    align_to_phmm,
    get_catalytic_residues
)

class VanadiumDependent:
    def __init__(self):
        self.fasta = str
        self.brominase_hits = align_to_phmm()
        self.iodinase_hits = align_to_phmm()
        self.selective_chloroperoxidase_hits = align_to_phmm()
        self.non_selective_chloroperoxidase_hits = align_to_phmm()

    def categorize_brominases(self):
        brominases = []
        vbpo_cat_1 = get_catalytic_residues(self.brominase_hits, VBPO["first_active_site"]["region"])
        vbpo_cat_2 = get_catalytic_residues(self.brominase_hits, VBPO["second_active_site"]["region"])
        for signature_1, signature_2 in vbpo_cat_1, vbpo_cat_2:
            if (re.search(VBPO["first_active_site"]["signature"], signature_1)
                and re.search(["second_active_site"]["signature"], signature_2)):
                brominases.append
        return brominases

    def intermol_categorize_brominases(self):
        brominanses_intermol = []
        molecular_bridges = get_catalytic_residues(self.brominase_hits, VBPO["intermolecular-bridges"]["region"])
        for hit in self.brominase_hits:
            for protein in molecular_bridges:
                if (re.search(VBPO["intermolecular-bridges"]["first_motif"], molecular_bridges[protein])
                    or re.search(VBPO["intermolecular-bridges"]["second_motif"], molecular_bridges[protein])):
                    brominanses_intermol.append(protein)

        return molecular_bridges

    def categorize_iodinases(self):
        iodinases = []
        catalytic_residues = get_catalytic_residues(self.iodinase_hits, VIPO["catalytic_residues"]["region"])
        for protein in catalytic_residues:
            if re.search(VIPO["catalytic_residues"]["signature"], catalytic_residues[protein]):
                iodinases.append(protein.decode("utf-8"))

        return iodinases

    def categorize_selective_chlorinases(self):
        selective_chlorinases = []
        catalytic_residues = get_catalytic_residues(self.selective_chloroperoxidase_hits, VCPO["selectivity_residues"]["region"])
        for protein in catalytic_residues:
            if re.search(VCPO["selectivity_residues"]["signature"], catalytic_residues[protein]):
                selective_chlorinases.append(protein.decode("utf-8"))

        return selective_chlorinases