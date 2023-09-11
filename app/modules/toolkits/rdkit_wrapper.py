from typing import List, Union
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem,
    Descriptors,
    QED,
    Lipinski,
    rdMolDescriptors,
    rdmolops,
    rdFingerprintGenerator,
    MACCSkeys,
)
from hosegen import HoseGenerator


def check_RO5_violations(molecule: any) -> int:
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


def get_rdkit_descriptors(molecule: any) -> Union[tuple, str]:
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
        Ro5Violations = check_RO5_violations(molecule)
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


def get_3d_conformers(molecule: any, depict=True) -> Chem.Mol:
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


def get_tanimoto_similarity_rdkit(
    mol1, mol2, fingerprinter="ECFP", radius=2, nBits=2048
) -> Union[float, str]:
    """
    Calculate the Tanimoto similarity index between two molecular structures represented as RDKit Mol objects.

    This function computes the Tanimoto similarity index, a measure of structural similarity, between two chemical compounds
    using various fingerprinting methods available in RDKit.

    Args:
        mol1 (Chem.Mol): The RDKit Mol object representing the first molecule.
        mol2 (Chem.Mol): The RDKit Mol object representing the second molecule.
        fingerprinter (str, optional): The type of fingerprint to use. Defaults to "ECFP".
        radius (int, optional): The radius parameter for ECFP fingerprints. Ignored for other fingerprint types.
        nBits (int, optional): The number of bits for fingerprint vectors. Ignored for MACCS keys.

    Returns:
        Union[float, str]: The Tanimoto similarity index between the two molecules if they are valid.
            If molecules are not valid, returns a string indicating an error.

    Note:
        - Supported fingerprinter options: "ECFP", "RDKit", "Atompairs", "MACCS".
        - ECFP fingerprints are based on atom environments up to a specified radius.
        - RDKit and Atom Pair fingerprints are based on different molecular descriptors.
        - MACCS keys are a fixed-length binary fingerprint.
    """
    if mol1 and mol2:
        if fingerprinter == "ECFP":
            # Generate Morgan fingerprints for each molecule
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits)
        elif fingerprinter == "RDKit":
            # Generate RDKit fingerprints for each molecule
            rdkgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=nBits)
            fp1 = rdkgen.GetFingerprint(mol1)
            fp2 = rdkgen.GetFingerprint(mol2)
        elif fingerprinter == "Atompairs":
            # Generate Atompairs fingerprints for each molecule
            apgen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=nBits)
            fp1 = apgen.GetFingerprint(mol1)
            fp2 = apgen.GetFingerprint(mol2)
        elif fingerprinter == "MACCS":
            # Generate MACCSkeys for each molecule
            fp1 = MACCSkeys.GenMACCSKeys(mol1)
            fp2 = MACCSkeys.GenMACCSKeys(mol2)
        else:
            return "Unsupported fingerprinter!"

        # Calculate the Tanimoto similarity between the fingerprints
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return similarity
    else:
        return "Check SMILES strings for Errors"


async def get_rdkit_HOSE_codes(molecule: any, noOfSpheres: int) -> List[str]:
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


def get_2d_mol(molecule: any) -> str:
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


def get_rdkit_CXSMILES(molecule: any) -> str:
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


def get_properties(sdf_file) -> dict:
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
