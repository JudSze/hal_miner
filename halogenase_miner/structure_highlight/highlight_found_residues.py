import pandas as pd
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
import numpy as np

ppdb = PandasPdb().read_pdb('test/test.pdb')

atoms = ppdb.df['ATOM'].copy()

print("Original columns:", atoms.columns.tolist())

residue_numbers = sorted(atoms['residue_number'].unique())
colors = plt.cm.viridis(np.linspace(0, 1, len(residue_numbers)))

color_map = {}
for i, res_num in enumerate(residue_numbers):
    color_map[res_num] = i / (len(residue_numbers) - 1) * 90  # Scale to 0-100

atoms['b_factor'] = atoms['residue_number'].map(color_map)

ppdb.df['ATOM'] = atoms
ppdb.df['ATOM']

ppdb.to_pdb(path='modified_test.pdb',
           records=['ATOM', 'HETATM', 'OTHERS'],
           append_newline=True)

# Loading from AlphaFold
from biopandas.mmcif import PandasMmcif

ppdb = PandasMmcif().fetch_mmcif(uniprot_id="Q5VSL9", source="alphafold2-v4")
ppdb.df['ATOM']