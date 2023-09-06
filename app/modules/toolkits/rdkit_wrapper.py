from typing import List, Union
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops
from hosegen import HoseGenerator


def checkRo5Violations(molecule: any) -> int:
    """
    Check the molecule for violations of Lipinski's Rule of Five.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        int: Number of Lipinski Rule violations.
    """
    num_of_violations = 0
    if Descriptors.MolLogP(molecule) > 5:
        num_of_violations += 1
    if Descriptors.MolWt(molecule) > 500:
        num_of_violations += 1
    if Lipinski.NumHAcceptors(molecule) > 10:
        num_of_violations += 1
    if Lipinski.NumHDonors(molecule) > 5:
        num_of_violations += 1
    return num_of_violations


def getRDKitDescriptors(molecule: any) -> Union[tuple, str]:
    """
    Calculate a selected set of molecular descriptors for the input SMILES string.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        dict: Dictionary of calculated molecular descriptors.
            If an error occurs during SMILES parsing, an error message is returned.
    """
    if molecule:
        AtomC = rdMolDescriptors.CalcNumAtoms(molecule)
        HeavyAtomsC = rdMolDescriptors.CalcNumHeavyAtoms(molecule)
        MolWt = "%.2f" % Descriptors.MolWt(molecule)
        ExactMolWt = "%.5f" % Descriptors.ExactMolWt(molecule)
        ALogP = "%.2f" % QED.properties(molecule).ALOGP
        NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)
        PSA = "%.2f" % rdMolDescriptors.CalcTPSA(molecule)
        HBA = Descriptors.NumHAcceptors(molecule)
        HBD = Descriptors.NumHDonors(molecule)
        Lipinski_HBA = Lipinski.NumHAcceptors(molecule)
        Lipinski_HBD = Lipinski.NumHDonors(molecule)
        Ro5Violations = checkRo5Violations(molecule)
        AromaticRings = rdMolDescriptors.CalcNumAromaticRings(molecule)
        QEDWeighted = "%.2f" % QED.qed(molecule)
        FormalCharge = rdmolops.GetFormalCharge(molecule)
        fsp3 = "%.3f" % rdMolDescriptors.CalcFractionCSP3(molecule)
        NumRings = rdMolDescriptors.CalcNumRings(molecule)
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


def get3Dconformers(molecule: any, depict=True) -> Chem.Mol:
    """
    Convert a SMILES string to an RDKit Mol object with 3D coordinates.

    Args:
        molecule (Chem.Mol): RDKit molecule object.
        depict (bool, optional): If True, returns the molecule's 3D structure in MolBlock format.
            If False, returns the 3D molecule without hydrogen atoms.

    Returns:
        str or rdkit.Chem.rdchem.Mol: If `depict` is True, returns the 3D structure in MolBlock format.
            Otherwise, returns an RDKit Mol object.

    """
    if molecule:
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule, maxAttempts=5000, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(molecule)
        except Exception:
            AllChem.EmbedMolecule(molecule, maxAttempts=5000, useRandomCoords=True)
        if depict:
            return Chem.MolToMolBlock(molecule)
        else:
            molecule = Chem.RemoveHs(molecule)
            return Chem.MolToMolBlock(molecule)


def getTanimotoSimilarityRDKit(mol1, mol2) -> Union[float, str]:
    """
    Calculate the Tanimoto similarity index between two SMILES strings using Morgan Fingerprints.

    Args:
        mol1 (Chem.Mol): First molecule object.
        mol2 (Chem.Mol): Second molecule object.

    Returns:
        float or str: Tanimoto similarity index if molecules are valid.
    """
    if mol1 and mol2:
        # Generate Morgan fingerprints for each molecule
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)

        # Calculate the Tanimoto similarity between the fingerprints
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return similarity
    else:
        return "Check SMILES strings for Errors"


async def getRDKitHOSECodes(molecule: any, noOfSpheres: int) -> List[str]:
    """
    Calculate and retrieve RDKit HOSE codes for a given SMILES string.
    This function takes a SMILES string as input and returns the calculated HOSE codes.

    Args:
        molecule (Chem.Mol): RDKit molecule object.
        no_of_spheres (int): Number of spheres for which to generate HOSE codes.

    Returns:
        List[str]: List of HOSE codes generated for each atom.

    Raises:
        ValueError: If the input SMILES string is empty or contains whitespace.
    """

    gen = HoseGenerator()
    hosecodes = []
    for i in range(0, len(molecule.GetAtoms()) - 1):
        hosecode = gen.get_Hose_codes(molecule, i, noOfSpheres)
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


def has_stereochemistry(molecule: any) -> bool:
    """
    Check if the given SMILES string contains stereochemistry information.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        bool: True if the SMILES contains stereochemistry information, False otherwise.
    """
    if molecule is None:
        return False

    for atom in molecule.GetAtoms():
        chiral_tag = atom.GetChiralTag()
        if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
            return True

    return False


def get2Dmol(molecule: any) -> str:
    """
    Generate a 2D Mol block representation from a given SMILES string.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        str: 2D Mol block representation.
        If an error occurs during SMILES parsing, an error message is returned.
    """

    if molecule:
        AllChem.Compute2DCoords(molecule)
        molfile = Chem.MolToMolBlock(molecule)
        return molfile


def getRDKitCXSMILES(molecule: any) -> str:
    """
    Generate CXSMILES representation with coordinates from a given SMILES string.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        str: CXSMILES representation with coordinates.
        If an error occurs during SMILES parsing, an error message is returned.
    """

    if molecule:
        AllChem.Compute2DCoords(molecule)
        return Chem.MolToCXSmiles(molecule)


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
