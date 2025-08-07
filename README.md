
# Halogenase detection
The workflow detects and categorizes halogenases across five enzyme families:
* Flavin-dependent
* S-Adenosyl-L-Methionine-dependent
* Nonheme-iron alphaketoglutarate dependent
* Dimetal-carboxylate
* Vanadium-dependent


# General categorization
For the flavin-dependent family, the general workflow relies on the common WxWxI and the Fx.Px.Sx.G motifs for the conventional members,
and for GxGxxA and Fx.Px.Sx.G motifs for the unconventional members e.g. AetF, JamD, VatD and PhmJ.

For the SAM-dependent family, it searches for the C-terminal motif (RNAA, YYGG and RNGG and YYAA).

For the NHFe-dependent enzymes it tries to detect the catalytic triad (HXG/A for halogenases and HXD/E for dioxygenases/hydroxylases).

For the Dimetal-carboxylate family, it relies on E..QE..H.

For the Vanadium-dependent enzymes, the workflow branches into three separate paths for chloroperoxidases, bromoperoxidases and iodoperoxidases.
For the chloroperoxidases it looks for K and S to identify selective enzymes. For the bromoperoxidases it looks for two motifs that contain the active site.
For iodoperoxidases, it searches for several catalytic residues.

# Mining features.
/In the halogenase_detection/motif_db/motifs.py file you will find all the motifs you can search for in your fasta file.
Use search_motif function from halogenase_detection/categorization/signature_search.py to specify which family and motif you want to look at.
It is recommended to run the general workflow by defining an EnzymeFamily("file.fasta") instance.
If you want to check the hits and motifs for each halogenase family in that EnzymeFamily instance, you have to define the family first.
E.g. your_instance.flavin_dependent or your_instance.vanadium_dependent and then choose the hits or the results of the categorization.

Moreover, in the halogenase_detection/motif_db/specific_enzymes.json file, you find individual enzymes and positions of confirmed or putative catalytic residues within
a certain pHMM. If you are unsure which enzyme you want to compare your target-sequences to, please navigate to the halogenase database.
You will find the name of the enzymes with more information, including the role of the residues andinformation about the halide-acceptance.

If you are using this workflow, please cite the source.

_Credits to Simon Shaw for providing a scalable solution and setting software development standards for the workflow in antiSMASH. The package is built on [PyHMMER (Larralde et al., 2022), a Python library binding to HMMER (Eddy, 2011)](https://pyhmmer.readthedocs.io/en/stable/index.html)._

__Citation__