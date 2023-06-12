from openbabel import openbabel as ob
from openbabel import pybel

def getOBCanonicalSMILES(smiles:str):
    """This function takes an input as a SMILES string and
    returns a Canonical SMILES.
    Args (str): SMILES string.
    Returns (str): Canonical SMILES string.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
        
    # Create an Open Babel molecule object
    mol = ob.OBMol()
    
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.ReadString(mol, smiles)
    
    canSMILES = conv.WriteString(mol)
    canSMILES = canSMILES.strip()  # Remove leading/trailing whitespace
    return canSMILES   


def getOBInChI(smiles:str,InChIKey:bool=False):
    """This function takes an input as a SMILES string and
    returns a InChI
    Args (str): SMILES string.
    Returns (str): InChI string.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")

    # Create an Open Babel molecule object
    mol = ob.OBMol()

    # Create OBConversion
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "inchi")
    conv.ReadString(mol, smiles)
    
    inchi = conv.WriteString(mol)
    inchi = inchi.strip()  # Remove leading/trailing whitespace
    if InChIKey:
        conv.SetOptions("K", conv.OUTOPTIONS)
        inchikey_ = conv.WriteString(mol).rstrip()
        return inchikey_
    return inchi

def getOBMol(smiles:str, threeD:bool=False):
    """This function takes an input as a SMILES string and
    returns a 2D/3D mol block.
    Args (str): SMILES string.
    Returns (str): Mol block (2D/3D).
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")

    if threeD:
        mol = pybel.readstring("smi", smiles)
        mol.addh()
        mol.make3D()
        mol.removeh()
        return mol.write("mol")
        
    # Create an Open Babel molecule object
    mol = ob.OBMol()
    
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "mol")
    conv.ReadString(mol, smiles)
    
    # Generate 2D coordinates
    obBuilder = ob.OBBuilder()
    obBuilder.Build(mol)    
    
    mol_block = conv.WriteString(mol)
    mol_block = mol_block.strip()  # Remove leading/trailing whitespace
    return mol_block



