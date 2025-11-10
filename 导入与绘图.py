
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
# %%


m = Chem.MolFromSmiles('C(C[C@@H](C(=O)O)N)CN=C(N)N')

img = Draw.MolToImage(m)

# %%

template = Chem.MolFromSmiles('c1nccc2n1ccc2')
AllChem.Compute2DCoords(template)

ms = [Chem.MolFromSmiles(smi) for smi in ('OCCc1ccn2cnccc12','C1CC1Oc1cc2ccncn2c1','CNC(=O)c1nccc2cccn12')]
for m in ms:
    _ = AllChem.GenerateDepictionMatching2DStructure(m,template)

