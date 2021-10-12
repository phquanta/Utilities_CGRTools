from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import rdChemReactions
import subprocess
import xml.etree.ElementTree as ET


s="[N:1].[H][O:2][H].[O:3]=[O:4]>>[O:2].[H][N:1]([H])[H].[H][O:4][O:3][H]"
s1 = AllChem.ReactionFromSmarts(s)
s2=AllChem.ReactionToSmarts(s1)
#s = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]')

#m = AllChem.ReactionFromSmarts(s)
print(s2)
#sma = Chem.MolToSmiles(m)

#print(sma)