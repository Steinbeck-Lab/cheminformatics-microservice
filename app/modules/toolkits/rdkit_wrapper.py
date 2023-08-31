from typing import List, Union
from chembl_structure_pipeline import standardizer
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops
from app.modules.toolkits.cdk_wrapper import getCDKSDGMol
from hosegen import HoseGenerator


def checkSMILES(smiles: str) -> Chem.Mol:
    """
    Check whether the input SMILES string is valid. If not, attempt to standardize
    the molecule using the ChEMBL standardization pipeline.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        Chem.Mol or None: Valid molecule object or None if an error occurs.
            If an error occurs during SMILES parsing, an error message is returned.
    """

    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)
    try:
        if mol:
            return mol
        elif standardizer.standardize_molblock(mol):
            return standardizer.standardize_molblock(mol)
    except Exception:
        return "Error reading SMILES string, check again."


def checkRo5Violations(mol) -> int:
    """
    Check the molecule for violations of Lipinski's Rule of Five.

    Args:
        mol (Chem.Mol): RDKit molecule object.

    Returns:
        int: Number of Lipinski Rule violations.
    """
    num_of_violations = 0
    if Descriptors.MolLogP(mol) > 5:
        num_of_violations += 1
    if Descriptors.MolWt(mol) > 500:
        num_of_violations += 1
    if Lipinski.NumHAcceptors(mol) > 10:
        num_of_violations += 1
    if Lipinski.NumHDonors(mol) > 5:
        num_of_violations += 1
    return num_of_violations


def getRDKitDescriptors(smiles: str) -> Union[tuple, str]:
    """
    Calculate a selected set of molecular descriptors for the input SMILES string.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        dict: Dictionary of calculated molecular descriptors.
            If an error occurs during SMILES parsing, an error message is returned.
    """
    mol = checkSMILES(smiles)
    if mol:
        AtomC = rdMolDescriptors.CalcNumAtoms(mol)
        HeavyAtomsC = rdMolDescriptors.CalcNumHeavyAtoms(mol)
        MolWt = "%.2f" % Descriptors.MolWt(mol)
        ExactMolWt = "%.5f" % Descriptors.ExactMolWt(mol)
        ALogP = "%.2f" % QED.properties(mol).ALOGP
        NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        PSA = "%.2f" % rdMolDescriptors.CalcTPSA(mol)
        HBA = Descriptors.NumHAcceptors(mol)
        HBD = Descriptors.NumHDonors(mol)
        Lipinski_HBA = Lipinski.NumHAcceptors(mol)
        Lipinski_HBD = Lipinski.NumHDonors(mol)
        Ro5Violations = checkRo5Violations(mol)
        AromaticRings = rdMolDescriptors.CalcNumAromaticRings(mol)
        QEDWeighted = "%.2f" % QED.qed(mol)
        FormalCharge = rdmolops.GetFormalCharge(mol)
        fsp3 = "%.3f" % rdMolDescriptors.CalcFractionCSP3(mol)
        NumRings = rdMolDescriptors.CalcNumRings(mol)
        VABCVolume = None
        return (
            AtomC,
            HeavyAtomsC,
            float(MolWt),
            float(ExactMolWt),
            float(ALogP),
            NumRotatableBonds,
            float(PSA),
            HBA,
            HBD,
            Lipinski_HBA,
            Lipinski_HBD,
            Ro5Violations,
            AromaticRings,
            float(QEDWeighted),
            FormalCharge,
            float(fsp3),
            NumRings,
            str(VABCVolume),
        )
    else:
        return "Error reading SMILES string, check again."


def get3Dconformers(smiles, depict=True) -> Chem.Mol:
    """
    Convert a SMILES string to an RDKit Mol object with 3D coordinates.

    Args:
        smiles (str): SMILES string.
        depict (bool, optional): If True, returns the molecule's 3D structure in MolBlock format.
            If False, returns the 3D molecule without hydrogen atoms.

    Returns:
        str or rdkit.Chem.rdchem.Mol: If `depict` is True, returns the 3D structure in MolBlock format.
            Otherwise, returns an RDKit Mol object.

    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles.__contains__("R"):
        mol_str = getCDKSDGMol(smiles)
        mol = Chem.MolFromMolBlock(mol_str)
    else:
        mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, maxAttempts=5000, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception as e:
            print(e)
            AllChem.EmbedMolecule(mol, maxAttempts=5000, useRandomCoords=True)
        if depict:
            return Chem.MolToMolBlock(mol)
        else:
            mol = Chem.RemoveHs(mol)
            return Chem.MolToMolBlock(mol)
    else:
        return "Error reading SMILES string, check again."


def getTanimotoSimilarityRDKit(smiles1, smiles2) -> Union[float, str]:
    """
    Calculate the Tanimoto similarity index between two SMILES strings using Morgan Fingerprints.

    Args:
        smiles1 (str): First SMILES string.
        smiles2 (str): Second SMILES string.

    Returns:
        float or str: Tanimoto similarity index if SMILES strings are valid.
            If an error occurs during SMILES parsing, an error message is returned.

    """
    # Create two example molecules
    mol1 = checkSMILES(smiles1)
    mol2 = checkSMILES(smiles2)
    if mol1 and mol2:
        # Generate Morgan fingerprints for each molecule
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)

        # Calculate the Tanimoto similarity between the fingerprints
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return similarity
    else:
        return "Check SMILES strings for Errors"


async def getRDKitHOSECodes(smiles: str, noOfSpheres: int) -> List[str]:
    """
    Calculate and retrieve RDKit HOSE codes for a given SMILES string.
    This function takes a SMILES string as input and returns the calculated HOSE codes.

    Args:
        smiles (str): The SMILES string of the molecule.
        no_of_spheres (int): Number of spheres for which to generate HOSE codes.

    Returns:
        List[str]: List of HOSE codes generated for each atom.

    Raises:
        ValueError: If the input SMILES string is empty or contains whitespace.
    """

    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)
    gen = HoseGenerator()
    hosecodes = []
    for i in range(0, len(mol.GetAtoms()) - 1):
        hosecode = gen.get_Hose_codes(mol, i, noOfSpheres)
        hosecodes.append(hosecode)
    return hosecodes


def is_valid_molecule(input_text) -> Union[str, bool]:
    """
    Check whether the input text represents a valid molecule in SMILES or Molblock format.

    Args:
        input_text (str): SMILES string or Molblock.

    Returns:
        str: "smiles" if the input is a valid SMILES, "mol" if the input is a valid Molblock, otherwise False.
    """
    try:
        molecule = Chem.MolFromSmiles(input_text)
        if molecule:
            return "smiles"
        else:
            molecule = Chem.MolFromMolBlock(input_text)
            if molecule:
                return "mol"
            else:
                return False
    except Exception:
        return False


def has_stereochemistry(smiles: str) -> bool:
    """
    Check if the given SMILES string contains stereochemistry information.

    Args:
        smiles (str): SMILES string.

    Returns:
        bool: True if the SMILES contains stereochemistry information, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return False

    for atom in mol.GetAtoms():
        chiral_tag = atom.GetChiralTag()
        if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
            return True

    return False


def get2Dmol(smiles: str) -> str:
    """
    Generate a 2D Mol block representation from a given SMILES string.

    Args:
        smiles (str): The input SMILES string.

    Returns:
        str: 2D Mol block representation.
        If an error occurs during SMILES parsing, an error message is returned.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        AllChem.Compute2DCoords(mol)
        molfile = Chem.MolToMolBlock(mol)
        return molfile
    else:
        return "Error reading SMILES string, check again."


def getRDKitCXSMILES(smiles: str) -> str:
    """
    Generate CXSMILES representation with coordinates from a given SMILES string.

    Args:
        smiles (str): The input SMILES string.

    Returns:
        str: CXSMILES representation with coordinates.
        If an error occurs during SMILES parsing, an error message is returned.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        AllChem.Compute2DCoords(mol)
        return Chem.MolToCXSmiles(mol)
    else:
        return "Error reading SMILES string, check again."


def getProperties(sdf_file) -> dict:
    """
    Extracts properties from a single molecule contained in an SDF file.

    This function uses the RDKit library to read an SDF (Structure-Data File) and extract properties
    from the first molecule in the file. It checks if the supplied SDF file contains a valid molecule
    and retrieves its properties as a dictionary.

    Args:
        sdf_file (str): The path to the SDF file containing the molecule.

    Returns:
        Dict or None: A dictionary containing the properties of the molecule. If the SDF file contains
        a valid molecule, the dictionary will have property names as keys and property values as values.
        If no valid molecule is found, or if there are no properties associated with the molecule, None
        is returned.

    Raises:
        ValueError: If the SDF file is not found or cannot be read.
    """
    # Create an SDMolSupplier to read the SDF file
    suppl = Chem.SDMolSupplier()
    suppl.SetData(sdf_file.encode("utf-8"))

    # Check if the SDF file contains a valid molecule
    if len(suppl) == 1 and suppl[0]:
        # Extract properties as a dictionary
        properties = suppl[0].GetPropsAsDict()
        return properties
    else:
        return {"Error": "No properties found"}
