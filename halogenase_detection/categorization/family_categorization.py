import re
from collections import Counter

from halogenase_detection.motif_db import motifs
from halogenase_detection.motif_db.motifs import (
    Profiles,
    VBPO,
    VCPO,
    VIPO,
    SAM,
    DIMETAL,
    FDHs,
    NONHEME_IRON
)

from halogenase_detection.categorization.signature_search import (
    align_to_phmm,
    get_catalytic_residues,
    search_motif
)

class FlavinDependent:
    def __init__(self, fasta):
        self.fasta = fasta
        self.general_hits = align_to_phmm(motifs.Profiles.fdh_all_conventional, self.fasta)
        self.unconventional_hits = align_to_phmm(motifs.Profiles.fdh_unconventional, self.fasta)

    def fdh_motif_based_categorization(self):
        motif_matches = []
        all_matches = []

        for motif in FDHs:
            motif_matches.append(search_motif(self.general_hits, FDHs, motif))

        for motif_match in motif_matches:
            all_matches.append(motif_match.keys())

        complete_matches = Counter(all_matches)
        return [protein for protein, count in complete_matches.items() if count == 3]

class DimetalCarboxylate:
    def __init__(self, fasta):
        self.fasta = fasta
        self.dimetal_hits = align_to_phmm(motifs.Profiles.dimetal_carboxylate, self.fasta)
        self.dimetal_motif_matches = search_motif(self.dimetal_hits, DIMETAL, "second_motif")

class SAMDependent:
    def __init__(self, fasta):
        self.fasta = fasta
        self.chlorinase_hits = align_to_phmm(motifs.Profiles.sam_chlorinases, self.fasta)
        self.fluorinase_hits = align_to_phmm(motifs.Profiles.sam_fluorinases, self.fasta)
        self.fluorinase_motif_matches = search_motif(self.fluorinase_hits, SAM, "c_terminal_motif")

class VanadiumDependent:
    def __init__(self, fasta):
        self.fasta = fasta
        self.general_hits = align_to_phmm(motifs.Profiles.vhpo_general, self.fasta)
        self.brominase_hits = align_to_phmm(motifs.Profiles.vhpo_bromoperoxidases, self.fasta)
        self.iodinase_hits = align_to_phmm(motifs.Profiles.vhpo_iodoperoxidases, self.fasta)
        self.selective_chloroperoxidase_hits = align_to_phmm(motifs.Profiles.vhpo_selective_chloroperoxidases, self.fasta)
        self.non_selective_chloroperoxidase_hits = align_to_phmm(motifs.Profiles.vhpo_non_selective_chloroperoxidases, self.fasta)

    def vhpo_categorize_brominases(self):
        brominases = []
        vbpo_cat_1 = get_catalytic_residues(self.brominase_hits, VBPO["first_active_site"]["region"])
        vbpo_cat_2 = get_catalytic_residues(self.brominase_hits, VBPO["second_active_site"]["region"])
        proteins = list({*vbpo_cat_1.keys(), *vbpo_cat_2.keys()})

        for protein in proteins:
            try:
                if (re.search(VBPO["first_active_site"]["signature"], vbpo_cat_1[protein])
                    and re.search(VBPO["second_active_site"]["signature"], vbpo_cat_2[protein])):
                    brominases.append(protein)
            except KeyError:
                print("protein not in matches")
        return brominases

    def vhpo_intermol_categorize_brominases(self):
        brominanses_intermol = []
        molecular_bridges = get_catalytic_residues(self.brominase_hits, VBPO["intermolecular_bridges"]["region"])
        for protein, signature in molecular_bridges.items():
            if (re.search(VBPO["intermolecular_bridges"]["first_motif"], signature)
                or re.search(VBPO["intermolecular_bridges"]["second_motif"], signature)):
                brominanses_intermol.append(protein)

        return brominanses_intermol

    def vhpo_categorize_iodinases(self):
        iodinases = []
        catalytic_residues = get_catalytic_residues(self.iodinase_hits, VIPO["catalytic_residues"]["region"])
        for protein, signature in catalytic_residues.items():
            if re.search(VIPO["catalytic_residues"]["signature"], signature):
                iodinases.append(protein.decode("utf-8"))

        return iodinases

    def vhpo_categorize_selective_chlorinases(self):
        selective_chlorinases = []
        catalytic_residues = get_catalytic_residues(self.selective_chloroperoxidase_hits, VCPO["selectivity_residues"]["region"])
        for protein, signature in catalytic_residues.items():
            if re.search(VCPO["selectivity_residues"]["signature"], signature):
                selective_chlorinases.append(protein.decode("utf-8"))

        return selective_chlorinases

class NonHemeIronDependent:
    def __init__(self, fasta):
        self.fasta = fasta
        self.variant_a_hits = align_to_phmm(motifs.Profiles.nhfe_variant_A, self.fasta)
        self.variant_b_hits = align_to_phmm(motifs.Profiles.nhfe_variant_B, self.fasta)
        self.indole_alkaloid_hits = align_to_phmm(motifs.Profiles.nhfe_indole_alkaloid, self.fasta)
        self.nucleotide_hits = align_to_phmm(motifs.Profiles.nhfe_nucleotide, self.fasta)

    def nonheme_halogenase_catalytic_triad(self, hits):
        return search_motif(hits, NONHEME_IRON, "halogenase")

    def nonheme_non_halogenase_catalytic_triad(self, hits):
        return search_motif(hits, NONHEME_IRON, "non_halogenase")

class EnzymeFamily(FlavinDependent, DimetalCarboxylate, SAMDependent, VanadiumDependent, NonHemeIronDependent):
    def __init__(self, fasta):
        self.fasta = fasta
        self.non_hemes = NonHemeIronDependent(fasta)
        self.vanadium_dependent = VanadiumDependent(fasta)
        self.dimetal_carboxylate = DimetalCarboxylate(fasta)
        self.flavin_dependent = FlavinDependent(fasta)
        self.sam_dependent = SAMDependent(fasta)

