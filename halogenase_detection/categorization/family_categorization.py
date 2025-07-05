import re
from collections import Counter

from halogenase_detection.motif_db.motifs import (
    Profiles,
    VBPO,
    VCPO,
    VIPO,
    SAM,
    DIMETAL,
    FDHs
)

from signature_search import (
    align_to_phmm,
    get_catalytic_residues,
    search_motif
)

class FlavinDependent:
    def __init__(self):
        self.fasta = str
        self.general_hits = align_to_phmm(Profiles.fdh_all_conventional, self.fasta)
        self.unconventional_hits = align_to_phmm(Profiles.fdh_unconventional, self.fasta)

    def motif_based_categorization(self):
        motif_matches = []
        all_matches = []

        for motif in FDHs:
            motif_matches.append(search_motif(self.general_hits, FDHs, motif))

        for motif_match in motif_matches:
            all_matches.append(motif_match.keys())

        complete_matches = Counter(all_matches)
        return [protein for protein, count in complete_matches.items() if count == 3]

class DimetalCarboxylate:
    def __init__(self):
        self.fasta = str
        self.dimetal_hits = align_to_phmm(Profiles.dimetal_carboxylate, self.fasta)
        self.dimetal_motif_matches = search_motif(self.dimetal_hits, DIMETAL, "second_motif")

class SAMDependent:
    def __init__(self):
        self.fasta = str
        self.chlorinase_hits = align_to_phmm(Profiles.sam_chlorinases, self.fasta)
        self.fluorinase_hits = align_to_phmm(Profiles.sam_fluorinases, self.fasta)
        self.fluorinase_motif_matches = search_motif(self.fluorinase_hits, SAM, "c_terminal_motif")

class VanadiumDependent:
    def __init__(self):
        self.fasta = str
        self.brominase_hits = align_to_phmm(Profiles.bromoperoxidases, self.fasta)
        self.iodinase_hits = align_to_phmm(Profiles.iodoperoxidases, self.fasta)
        self.selective_chloroperoxidase_hits = align_to_phmm(Profiles.selective_chloroperoxidases, self.fasta)
        self.non_selective_chloroperoxidase_hits = align_to_phmm(Profiles.non_selective_chloroperoxidases, self.fasta)

    def categorize_brominases(self):
        brominases = []
        vbpo_cat_1 = get_catalytic_residues(self.brominase_hits, VBPO["first_active_site"]["region"])
        vbpo_cat_2 = get_catalytic_residues(self.brominase_hits, VBPO["second_active_site"]["region"])
        for signature_1, signature_2 in vbpo_cat_1, vbpo_cat_2:
            if (re.search(VBPO["first_active_site"]["signature"], vbpo_cat_1[signature_1])
                and re.search(["second_active_site"]["signature"], vbpo_cat_2[signature_2])):
                brominases.append
        return brominases

    def intermol_categorize_brominases(self):
        brominanses_intermol = []
        molecular_bridges = get_catalytic_residues(self.brominase_hits, VBPO["intermolecular-bridges"]["region"])
        for protein, signature in molecular_bridges.items():
            if (re.search(VBPO["intermolecular-bridges"]["first_motif"], signature)
                or re.search(VBPO["intermolecular-bridges"]["second_motif"], signature)):
                brominanses_intermol.append(protein)

        return molecular_bridges

    def categorize_iodinases(self):
        iodinases = []
        catalytic_residues = get_catalytic_residues(self.iodinase_hits, VIPO["catalytic_residues"]["region"])
        for protein, signature in catalytic_residues.items():
            if re.search(VIPO["catalytic_residues"]["signature"], signature):
                iodinases.append(protein.decode("utf-8"))

        return iodinases

    def categorize_selective_chlorinases(self):
        selective_chlorinases = []
        catalytic_residues = get_catalytic_residues(self.selective_chloroperoxidase_hits, VCPO["selectivity_residues"]["region"])
        for protein, signature in catalytic_residues.items():
            if re.search(VCPO["selectivity_residues"]["signature"], signature):
                selective_chlorinases.append(protein.decode("utf-8"))

        return selective_chlorinases
