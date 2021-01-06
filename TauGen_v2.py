# imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

# function to generate tautomers
def TauGen(mol):
	enumerator = rdMolStandardize.TautomerEnumerator()
    return enumerator.Enumerate(mol)

def CanonicalTauGen(mol):
	enumerator = rdMolStandardize.TautomerEnumerator()
	return enumerator.Canonicalize(mol)

def TauGenFromSmiles(smi):
    mol = Chem.MolFromSmiles(smi)
    return TauGen(mol)

def CanonicalTauGenFromSmiles(smi):
    mol = Chem.MolFromSmiles(smi)
    return CanonicalTauGen(mol)

# function for 3D tautomers generation
def TauGen3D(mol):
    res = []
    for mol in TauGen(mol):
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)
        AllChem.MMFFOptimizeMolecule(mol2)
        res.append(mol2)
        
    return res

# test
mol = Chem.MolFromSmiles('NC1=NNC(C1)=O')
my = TauGen3D(mol)
for i in my:
    new = Chem.RemoveHs(i)
    print(Chem.MolToSmiles(new))