from openbabel import openbabel as ob
from openbabel import pybel


def getOBCanonicalSMILES(smiles: str) -> str:
    """
    Convert a SMILES string to Canonical SMILES.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        str: Canonical SMILES string.
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


def getOBInChI(smiles: str, InChIKey: bool = False) -> str:
    """
    Convert a SMILES string to InChI.

    Args:
        smiles (str): Input SMILES string.
        InChIKey (bool, optional): Whether to return InChIKey. Defaults to False.

    Returns:
        str: InChI string or InChIKey string if InChIKey is True.
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


def getOBMol(smiles: str, threeD: bool = False, depict: bool = False) -> str:
    """
    Convert a SMILES string to a 2D/3D mol block.

    Args:
        smiles (str): Input SMILES string.
        threeD (bool, optional): Generate 3D structure. Defaults to False.
        depict (bool, optional): Generate 3D structure for depiction. Defaults to False.

    Returns:
        str: Mol block (2D/3D).
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")

    if threeD:
        mol = pybel.readstring("smi", smiles)
        mol.addh()
        mol.make3D()
        gen3d = ob.OBOp.FindType("gen3D")
        gen3d.Do(mol.OBMol, "--best")
        if depict:
            return mol.write("mol")
        else:
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
