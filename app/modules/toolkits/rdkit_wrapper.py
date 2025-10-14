from __future__ import annotations

from typing import List
from typing import Tuple
from typing import Union

from hosegen import HoseGenerator
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import MACCSkeys
from rdkit.Chem import QED
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
from rdkit.Chem.FilterCatalog import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Contrib.efgs import efgs
from rdkit.Contrib.IFG import ifg
from rdkit.Contrib.SA_Score import sascorer
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator
from mapchiral.mapchiral import encode, jaccard_similarity


def check_RO5_violations(molecule: any) -> int:
    """Check the molecule for violations of Lipinski's Rule of Five.

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


def check_RO5_violations_detailed(molecule: any) -> dict:
    """Check the molecule for violations of Lipinski's Rule of Five with detailed information.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        dict: Dictionary containing violation details with keys:
            - violations: int (number of violations)
            - details: list of violation descriptions
            - properties: dict of actual property values
            - passes: bool (True if no violations)
    """
    violations = []
    properties = {}

    # Calculate properties
    mw = Descriptors.MolWt(molecule)
    logp = Descriptors.MolLogP(molecule)
    hba = Lipinski.NumHAcceptors(molecule)
    hbd = Lipinski.NumHDonors(molecule)

    properties = {
        "molecular_weight": round(mw, 2),
        "logp": round(logp, 2),
        "hb_acceptors": hba,
        "hb_donors": hbd,
    }

    # Check violations
    if logp > 5:
        violations.append(f"LogP = {logp:.2f} (> 5)")
    if mw > 500:
        violations.append(f"MW = {mw:.1f} Da (> 500)")
    if hba > 10:
        violations.append(f"HBA = {hba} (> 10)")
    if hbd > 5:
        violations.append(f"HBD = {hbd} (> 5)")

    return {
        "violations": len(violations),
        "details": violations,
        "properties": properties,
        "passes": len(violations) == 0,
    }


def get_MolVolume(molecule: any) -> float:
    """
    Calculate the volume of a molecule.

    This function calculates the volume of a molecule using RDKit's molecular modeling functionalities.
    It adds hydrogens to the molecule, embeds it into 3D space, and computes the molecular volume.

    Args:
        molecule (any): The molecule for which the volume needs to be calculated.

    Returns:
        float: The volume of the molecule.
    """
    molecule = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(molecule, useRandomCoords=True)
    volume = AllChem.ComputeMolVolume(molecule, gridSpacing=0.2)
    return volume


def get_rdkit_descriptors(molecule: any) -> Union[tuple, str]:
    """Calculate a selected set of molecular descriptors for the input SMILES.

    string.

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
        VABCVolume = "%.2f" % get_MolVolume(molecule)
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
            float(VABCVolume),
        )


def get_3d_conformers(molecule: any, depict=True) -> Chem.Mol:
    """Convert a SMILES string to an RDKit Mol object with 3D coordinates.

    Args:
        molecule (Chem.Mol): RDKit molecule object.
        depict (bool, optional): If True, returns the molecule's 3D structure in MolBlock format. If False, returns the 3D molecule without hydrogen atoms.

    Returns:
        str or rdkit.Chem.rdchem.Mol: If `depict` is True, returns the 3D structure in MolBlock format. Otherwise, returns an RDKit Mol object.
    """
    if molecule:
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule, maxAttempts=5000, useRandomCoords=True)
        try:
            AllChem.MMFFOptimizeMolecule(molecule)
        except Exception:
            AllChem.EmbedMolecule(
                molecule,
                maxAttempts=5000,
                useRandomCoords=True,
            )
        if depict:
            return Chem.MolToMolBlock(molecule)
        else:
            molecule = Chem.RemoveHs(molecule)
            return Chem.MolToMolBlock(molecule)


def get_tanimoto_similarity_rdkit(
    mol1,
    mol2,
    fingerprinter="ECFP",
    radius=2,
    nBits=2048,
) -> Union[float, str]:
    """Calculate the Tanimoto similarity index between two molecular.

    structures.

    represented as RDKit Mol objects.

    This function computes the Tanimoto similarity index, a measure of structural similarity, between two chemical compounds
    using various fingerprinting methods available in RDKit.

    Args:
        mol1 (Chem.Mol): The RDKit Mol object representing the first molecule.
        mol2 (Chem.Mol): The RDKit Mol object representing the second molecule.
        fingerprinter (str, optional): The type of fingerprint to use. Options are "ECFP", "RDKit", "AtomPairs", "MACCS". Defaults to "ECFP".
        radius (int, optional): The radius parameter for ECFP fingerprints (e.g. radius 2 for generating ECFP4 fingerprints, default value).
        Ignored for all other fingerprinter options than "ECFP".

    Returns:
        Union[float, str]: The Tanimoto similarity index between the two molecules if they are valid. If molecules are not valid, returns a string indicating an error.

    Note:
        - Supported fingerprinter options: "ECFP", "RDKit", "Atompairs", "MACCS".
        - ECFP fingerprints are based on atom environments up to a specified radius.
        - RDKit and Atom Pair fingerprints are based on different molecular descriptors.
        - MACCS keys are a fixed-length binary fingerprint.
        - MAPC (MinHashed Atom-Pair Fingerprint Chiral): https://github.com/reymond-group/mapchiral
    """
    if mol1 and mol2:
        if fingerprinter == "ECFP":
            # Generate Morgan fingerprints for each molecule
            morgan_fps = rdFingerprintGenerator.GetMorganGenerator(
                radius, fpSize=nBits, includeChirality=True
            )
            fp1 = morgan_fps.GetFingerprint(mol1)
            fp2 = morgan_fps.GetFingerprint(mol2)
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
        elif fingerprinter == "MAPC":
            # Generate MAPC for each molecule
            fp1 = encode(mol1, max_radius=radius, n_permutations=nBits, mapping=False)
            fp2 = encode(mol2, max_radius=radius, n_permutations=nBits, mapping=False)
            similarity = jaccard_similarity(fp1, fp2)
            return similarity
        else:
            return "Unsupported fingerprinter!"

        # Calculate the Tanimoto similarity between the fingerprints
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return similarity
    else:
        return "Check SMILES strings for Errors"


async def get_rdkit_HOSE_codes(molecule: any, noOfSpheres: int) -> List[str]:
    """Calculate and retrieve RDKit HOSE codes for a given SMILES string.

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
    """Check whether the input text represents a valid molecule in SMILES or.

    Molblock format.

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


def has_stereo_defined(molecule: Chem.Mol) -> bool:
    """
    Checks if a molecular structure represented by an RDKit molecule object has any chiral centers defined.

    Args:
        molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.

    Returns:
        bool: True if the molecule has at least one chiral center defined, False otherwise.
    """
    if molecule is None:
        return False

    for atom in molecule.GetAtoms():
        chiral_tag = atom.GetChiralTag()
        if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
            return True

    return False


def has_potential_stereochemistry(molecule: Chem.Mol) -> bool:
    """
    Checks if a molecular structure represented by an RDKit molecule object has any stereochemistry information.

    Args:
        molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.

    Returns:
        bool: True if the molecule has stereochemistry information, False otherwise.

    This function uses the RDKit's FindPotentialStereo function to identify potential stereochemistry information
    in the molecule. If any stereochemistry information is found, the function returns True, otherwise False.
    """

    if molecule is None:
        return False

    stereo_info = Chem.FindPotentialStereo(molecule)
    if len(list(stereo_info)) > 0:
        return True
    else:
        return False


def get_2d_mol(molecule: any) -> str:
    """Generate a 2D Mol block representation from a given SMILES string.

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
    """Generate CXSMILES representation with coordinates from a given SMILES.

    string.

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
    """Extracts properties from a single molecule contained in an SDF file.

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


def get_sas_score(molecule: any) -> float:
    """Calculate the Synthetic Accessibility Score (SAS) for a given molecule.

    The Synthetic Accessibility Score is a measure of how easy or difficult it is to synthesize a given molecule.
    A higher score indicates a molecule that is more challenging to synthesize, while a lower score suggests a molecule
    that is easier to synthesize.

    Parameters:
        molecule (Chem.Mol): An RDKit molecule object representing the chemical structure.

    Returns:
        float: The Synthetic Accessibility Score rounded to two decimal places.

    Note:
        - The SAS is calculated using the sascorer.calculateScore() function from the RDKit Contrib library.
        - The SAS score can be used as a factor in drug design and compound optimization, with lower scores often
          indicating more drug-like and synthesizable molecules.

    See Also:
        - RDKit Contrib: https://rdkit.org/docs_contribs/index.html

    References:
        - Ertl, P., & Schuffenhauer, A. (2009). Estimation of synthetic accessibility score of drug-like molecules based
          on molecular complexity and fragment contributions. Journal of Cheminformatics, 1(1), 8.
          DOI: 10.1186/1758-2946-1-8
        - RDKit Documentation: https://www.rdkit.org/docs/index.html
    """
    if molecule:
        sas_score = sascorer.calculateScore(molecule)
        return round(sas_score, 2)


def get_PAINS(molecule: any) -> Union[bool, Tuple[str, str]]:
    """Check if a molecule contains a PAINS (Pan Assay INterference compoundS)substructure.

    Parameters:
    molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
    Union[bool, Tuple[str, str]]: The function returns a tuple with the PAINS family and its description if a PAINS substructure is detected in the molecule. Otherwise, it returns False.

    This function uses the RDKit library to check if the given molecule contains
    any PAINS substructure. PAINS are known substructures that may interfere
    with various biological assays.
    """
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    entry = catalog.GetFirstMatch(molecule)
    if entry:
        family = entry.GetProp("Scope")
        description = entry.GetDescription().capitalize()
        return family, description
    else:
        return False


def get_PAINS_detailed(molecule: any) -> dict:
    """Check if a molecule contains a PAINS substructure with detailed information.

    Parameters:
    molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
    dict: Dictionary containing PAINS analysis with keys:
        - contains_pains: bool (True if PAINS found - this is BAD for drug-likeness)
        - family: str or None (PAINS family if found)
        - description: str or None (PAINS description if found)
        - passes: bool (True if NO PAINS found - this is GOOD for drug-likeness)
        - details: str (human-readable explanation)
    """
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    entry = catalog.GetFirstMatch(molecule)
    if entry:
        family = entry.GetProp("Scope")
        description = entry.GetDescription().capitalize()
        return {
            "contains_pains": True,
            "family": family,
            "description": description,
            "passes": False,  # Finding PAINS is BAD, so passes = False
            "details": f"PAINS match found: {family} - {description}",
        }
    else:
        return {
            "contains_pains": False,
            "family": None,
            "description": None,
            "passes": True,  # No PAINS found is GOOD, so passes = True
            "details": "No PAINS substructures detected",
        }


def get_GhoseFilter(molecule: any) -> bool:
    """Determine if a molecule satisfies Ghose's filter criteria.

    Ghose's filter is a set of criteria for drug-like molecules.
    This function checks if a given molecule meets the criteria defined by Ghose.

    Parameters:
    molecule (any):  A molecule represented as an RDKit Mol object.

    Returns:
    bool: True if the molecule meets Ghose's criteria, False otherwise.

    Ghose's criteria:
    - Molecular Weight (MW) should be between 160 and 480.
    - LogP (Partition Coefficient) should be between 0.4 and 5.6.
    - Number of Atoms (NoAtoms) should be between 20 and 70.
    - Molar Refractivity (MolarRefractivity) should be between 40 and 130.
    """
    MW = Descriptors.ExactMolWt(molecule)
    logP = Descriptors.MolLogP(molecule)
    NoAtoms = rdMolDescriptors.CalcNumAtoms(molecule)
    MolarRefractivity = Chem.Crippen.MolMR(molecule)

    # Check if the molecule satisfies Ghose's criteria.
    if (
        (160 <= MW <= 480)
        and (0.4 <= logP <= 5.6)
        and (20 <= NoAtoms <= 70)
        and (40 <= MolarRefractivity <= 130)
    ):
        return True
    else:
        return False


def get_GhoseFilter_detailed(molecule: any) -> dict:
    """Determine if a molecule satisfies Ghose's filter criteria with detailed information.

    Parameters:
    molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
    dict: Dictionary containing Ghose filter analysis with keys:
        - passes: bool (True if passes Ghose criteria)
        - violations: list of violation descriptions
        - properties: dict of actual property values
        - details: str (human-readable explanation)
    """
    MW = Descriptors.ExactMolWt(molecule)
    logP = Descriptors.MolLogP(molecule)
    NoAtoms = rdMolDescriptors.CalcNumAtoms(molecule)
    MolarRefractivity = Chem.Crippen.MolMR(molecule)

    violations = []
    if not (160 <= MW <= 480):
        violations.append(f"MW = {MW:.1f} (not in 160-480)")
    if not (0.4 <= logP <= 5.6):
        violations.append(f"LogP = {logP:.2f} (not in 0.4-5.6)")
    if not (20 <= NoAtoms <= 70):
        violations.append(f"Atoms = {NoAtoms} (not in 20-70)")
    if not (40 <= MolarRefractivity <= 130):
        violations.append(f"MR = {MolarRefractivity:.1f} (not in 40-130)")

    return {
        "passes": len(violations) == 0,
        "violations": violations,
        "properties": {
            "molecular_weight": round(MW, 1),
            "logp": round(logP, 2),
            "atom_count": NoAtoms,
            "molar_refractivity": round(MolarRefractivity, 1),
        },
        "details": "No violations" if len(violations) == 0 else "; ".join(violations),
    }


def get_VeberFilter(molecule: any) -> bool:
    """Apply the Veber filter to evaluate the drug-likeness of a molecule.

    The Veber filter assesses drug-likeness based on two criteria: the number of
    rotatable bonds and the polar surface area (TPSA). A molecule is considered
    drug-like if it has 10 or fewer rotatable bonds and a TPSA of 140 or less.

    Parameters:
        molecule (any):  A molecule represented as an RDKit Mol object.

    Returns:
        bool: True if the molecule passes the Veber filter criteria, indicating
          drug-likeness; False otherwise.

    Note:
    The function relies on RDKit functions to calculate the number of rotatable
    bonds and TPSA, and it returns a boolean value to indicate whether the input
    molecule passes the Veber filter criteria.

    Reference:
    Veber, D. F., Johnson, S. R., Cheng, H. Y., Smith, B. R., Ward, K. W., & Kopple,
    K. D. (2002). Molecular properties that influence the oral bioavailability of
    drug candidates. Journal of Medicinal Chemistry, 45(12), 2615-2623.
    DOI: 10.1021/jm020017n
    """
    NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)
    tpsa = Descriptors.TPSA(molecule)
    if NumRotatableBonds <= 10 and tpsa <= 140:
        return True
    else:
        return False


def get_VeberFilter_detailed(molecule: any) -> dict:
    """Apply the Veber filter with detailed information about violations.

    Parameters:
        molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
        dict: Dictionary containing Veber filter analysis with keys:
            - passes: bool (True if passes Veber criteria)
            - violations: list of violation descriptions
            - properties: dict of actual property values
            - details: str (human-readable explanation)
    """
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)
    tpsa = Descriptors.TPSA(molecule)

    violations = []
    if rotatable_bonds > 10:
        violations.append(f"Rotatable bonds = {rotatable_bonds} (> 10)")
    if tpsa > 140:
        violations.append(f"TPSA = {tpsa:.1f} (> 140)")

    return {
        "passes": len(violations) == 0,
        "violations": violations,
        "properties": {"rotatable_bonds": rotatable_bonds, "tpsa": round(tpsa, 1)},
        "details": "No violations" if len(violations) == 0 else "; ".join(violations),
    }


def get_REOSFilter(molecule: any) -> bool:
    """Determine if a molecule passes the REOS (Rapid Elimination Of Swill).

    filter.

    The REOS filter is a set of criteria that a molecule must meet to be considered
    a viable drug-like compound. This function takes a molecule as input and checks
    its properties against the following criteria:

    - Molecular Weight (MW): Must be in the range [200, 500].
    - LogP (Partition Coefficient): Must be in the range [-5, 5].
    - Hydrogen Bond Donors (HBD): Must be in the range [0, 5].
    - Hydrogen Bond Acceptors (HBA): Must be in the range [0, 10].
    - Formal Charge: Must be in the range [-2, 2].
    - Number of Rotatable Bonds: Must be in the range [0, 8].
    - Number of Heavy Atoms (non-hydrogen atoms): Must be in the range [15, 50].

    Parameters:
        molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
        bool: True if the molecule passes the REOS filter, False otherwise.
    """
    MW = Descriptors.ExactMolWt(molecule)
    logP = Descriptors.MolLogP(molecule)
    HBD = Descriptors.NumHDonors(molecule)
    HBA = Descriptors.NumHAcceptors(molecule)
    FormalCharge = rdmolops.GetFormalCharge(molecule)
    NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)
    HeavyAtomsC = rdMolDescriptors.CalcNumHeavyAtoms(molecule)

    if (
        200 <= MW <= 500
        and -5 <= logP <= 5
        and 0 <= HBD <= 5
        and 0 <= HBA <= 10
        and -2 <= FormalCharge <= 2
        and 0 <= NumRotatableBonds <= 8
        and 15 <= HeavyAtomsC <= 50
    ):
        return True
    else:
        return False


def get_REOSFilter_detailed(molecule: any) -> dict:
    """Determine if a molecule passes the REOS filter with detailed information.

    Parameters:
        molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
        dict: Dictionary containing REOS filter analysis with keys:
            - passes: bool (True if passes REOS criteria)
            - violations: list of violation descriptions
            - properties: dict of actual property values
            - details: str (human-readable explanation)
    """
    MW = Descriptors.ExactMolWt(molecule)
    logP = Descriptors.MolLogP(molecule)
    HBD = Descriptors.NumHDonors(molecule)
    HBA = Descriptors.NumHAcceptors(molecule)
    FormalCharge = rdmolops.GetFormalCharge(molecule)
    NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)
    HeavyAtomsC = rdMolDescriptors.CalcNumHeavyAtoms(molecule)

    violations = []
    if not (200 <= MW <= 500):
        violations.append(f"MW = {MW:.1f} (not in 200-500)")
    if not (-5 <= logP <= 5):
        violations.append(f"LogP = {logP:.2f} (not in -5 to 5)")
    if not (0 <= HBD <= 5):
        violations.append(f"HBD = {HBD} (not in 0-5)")
    if not (0 <= HBA <= 10):
        violations.append(f"HBA = {HBA} (not in 0-10)")
    if not (-2 <= FormalCharge <= 2):
        violations.append(f"Charge = {FormalCharge} (not in -2 to 2)")
    if not (0 <= NumRotatableBonds <= 8):
        violations.append(f"RotBonds = {NumRotatableBonds} (not in 0-8)")
    if not (15 <= HeavyAtomsC <= 50):
        violations.append(f"HeavyAtoms = {HeavyAtomsC} (not in 15-50)")

    return {
        "passes": len(violations) == 0,
        "violations": violations,
        "properties": {
            "molecular_weight": round(MW, 1),
            "logp": round(logP, 2),
            "hb_donors": HBD,
            "hb_acceptors": HBA,
            "formal_charge": FormalCharge,
            "rotatable_bonds": NumRotatableBonds,
            "heavy_atoms": HeavyAtomsC,
        },
        "details": "No violations" if len(violations) == 0 else "; ".join(violations),
    }


def get_RuleofThree(molecule: any) -> bool:
    """Check if a molecule meets the Rule of Three criteria.

    The Rule of Three is a guideline for drug-likeness in chemical compounds.
    It suggests that a molecule is more likely to be a good drug candidate if it
    meets the following criteria:
    1. Molecular Weight (MW) <= 300
    2. LogP (partition coefficient) <= 3
    3. Number of Hydrogen Bond Donors (HBD) <= 3
    4. Number of Hydrogen Bond Acceptors (HBA) <= 3
    5. Number of Rotatable Bonds <= 3

    Parameters:
        molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
        bool: True if the molecule meets the Rule of Three criteria, False otherwise.
    """
    MW = Descriptors.ExactMolWt(molecule)
    logP = Descriptors.MolLogP(molecule)
    HBD = Descriptors.NumHDonors(molecule)
    HBA = Descriptors.NumHAcceptors(molecule)
    NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)

    if MW <= 300 and logP <= 3 and HBD <= 3 and HBA <= 3 and NumRotatableBonds <= 3:
        return True
    else:
        return False


def get_RuleofThree_detailed(molecule: any) -> dict:
    """Check if a molecule meets the Rule of Three criteria with detailed information.

    Parameters:
        molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
        dict: Dictionary containing Rule of Three analysis with keys:
            - passes: bool (True if passes Rule of Three criteria)
            - violations: list of violation descriptions
            - properties: dict of actual property values
            - details: str (human-readable explanation)
    """
    MW = Descriptors.ExactMolWt(molecule)
    logP = Descriptors.MolLogP(molecule)
    HBD = Descriptors.NumHDonors(molecule)
    HBA = Descriptors.NumHAcceptors(molecule)
    NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)

    violations = []
    if MW > 300:
        violations.append(f"MW = {MW:.1f} (> 300)")
    if logP > 3:
        violations.append(f"LogP = {logP:.2f} (> 3)")
    if HBD > 3:
        violations.append(f"HBD = {HBD} (> 3)")
    if HBA > 3:
        violations.append(f"HBA = {HBA} (> 3)")
    if NumRotatableBonds > 3:
        violations.append(f"RotBonds = {NumRotatableBonds} (> 3)")

    return {
        "passes": len(violations) == 0,
        "violations": violations,
        "properties": {
            "molecular_weight": round(MW, 1),
            "logp": round(logP, 2),
            "hb_donors": HBD,
            "hb_acceptors": HBA,
            "rotatable_bonds": NumRotatableBonds,
        },
        "details": "No violations" if len(violations) == 0 else "; ".join(violations),
    }


def get_ertl_functional_groups_ifg(molecule: any) -> list:
    """
    This function takes an organic molecule as input and uses the algorithm
    proposed by Peter Ertl to identify functional groups within the molecule.
    The implementation by Richard Hall and Guillaume Godin (IFG) is used here
    (https://github.com/rdkit/rdkit/tree/master/Contrib/IFG).

    Parameters:
    - molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
    - list: A list of dictionaries (one dict for each identified functional group in the molecule) with structured data about each detected functional group (FG).
      For each dictionary, the "atomIds" key value contains the atom IDs of the FG atoms marked according to the Ertl algorithm.
      "atoms" contains a SMILES representation of the FG (again, only the marked atoms) and "type"
      contains a SMILES representation of the FG including also the unmarked, environmental carbon atoms connected to the marked atoms.
      "description" contains the full string representation of the IFG object for display purposes.
      Example: "[{'atomIds': [1, 3, 21], 'atoms': 'OCO', 'type': 'CC1(C)OCCO1', 'description': "IFG(atomIds=(1, 3, 21), atoms='OCO', type='CC1(C)OCCO1')"}, {...}]"
      The IFG implementation does not provide a generalization of the FG environments.
      If no functional groups are found, the function returns a list with a single element:
      "[{'None': 'No fragments found'}]"

    References:
    - Ertl, P. An algorithm to identify functional groups in organic molecules. J Cheminform 9, 36 (2017). https://doi.org/10.1186/s13321-017-0225-z
    """
    if molecule:
        fragments = ifg.identify_functional_groups(molecule)
        if fragments:
            # Convert IFG objects to structured dictionaries for better frontend handling
            structured_groups = []
            for fragment in fragments:
                try:
                    # Extract information from IFG object
                    group_data = {
                        "atomIds": (
                            list(fragment.atomIds)
                            if hasattr(fragment, "atomIds")
                            else []
                        ),
                        "atoms": (
                            str(fragment.atoms) if hasattr(fragment, "atoms") else ""
                        ),
                        "type": str(fragment.type) if hasattr(fragment, "type") else "",
                        "description": str(
                            fragment
                        ),  # Full string representation for display
                    }
                    structured_groups.append(group_data)
                except Exception:
                    # Fallback to string representation if structured extraction fails
                    structured_groups.append(
                        {
                            "atomIds": [],
                            "atoms": "",
                            "type": "",
                            "description": str(fragment),
                        }
                    )
            return structured_groups
        else:
            return [{"None": "No fragments found"}]
    else:
        return [{"Error": "Check input SMILES"}]


# TODO: finish this implementation and add it to the router as another option; add tests also!
def get_ertl_functional_groups_efgs(molecule: any) -> list:
    """
    This function takes an organic molecule as input and uses the algorithm
    proposed by Peter Ertl to identify functional groups within the molecule.
    The implementation by Gonzalo Colmenarejo (EFGs) is used here
    (https://github.com/rdkit/rdkit/tree/master/Contrib/efgs).

    Parameters:
    - molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
    - list: A list of dictionaries (one dict for each identified functional group in the molecule) with structured data about each detected functional group (FG).
      For each dictionary, the "fgs" key value contains the atom IDs of the FG atoms marked according to the Ertl algorithm.
      "psmis" contains a (canonical) "pseudo" (see Ertl's paper) SMILES representation of the "generalized" (see Ertl's paper) FG and "fg_mols"
      contains an RDKit Mol object of the generalized FG.
      "description" contains the full string representation of the EFGs object for display purposes.
      Example: "[{'fgs': [1, 3, 21], 'psmis': '[R][O]C[O][R]', 'fg_mols': <rdkit.Chem.rdchem.Mol object at 0x00000213936BA030>, 'description': "EFGS([[1, 3, 21], '[R][O]C[O][R]', <rdkit.Chem.rdchem.Mol object at 0x00000213936BA030>])"}, {...}]"
      If no functional groups are found, the function returns a list with a single element:
      "[{'None': 'No fragments found'}]"

    References:
    - Ertl, P. An algorithm to identify functional groups in organic molecules. J Cheminform 9, 36 (2017). https://doi.org/10.1186/s13321-017-0225-z
    - Colmenarejo, G. EFGs: A Complete and Accurate Implementation of Ertl's Functional Group Detection Algorithm in RDKit. J Chem Inf Model 65, 3 (2025). https://doi.org/10.1021/acs.jcim.4c02268
    """
    if molecule:
        # we ignore the img_text output here
        img_text, fgs, psmis, fg_mols = efgs.get_dec_fgs(molecule)
        if fgs:
            # Convert EFGS objects to structured dictionaries for better frontend handling
            structured_groups = []
            for i in range(len(fgs)):
                try:
                    # Extract information from EFGS object
                    group_data = {
                        "fgs": (
                            list(fgs[i])
                            if fgs[i]
                            else []
                        ),
                        "psmis": (
                            str(psmis[i]) if psmis[i] else ""
                        ),
                        "fg_mols": fg_mols[i] if fg_mols[i] else None,
                        "description": "EFGS(" + str([
                            list(fgs[i]) if fgs[i] else [],
                            str(psmis[i]) if psmis[i] else "",
                            fg_mols[i] if fg_mols[i] else None,
                        ]) + ")",  # Full string representation for display
                    }
                    structured_groups.append(group_data)
                except Exception as e:
                    # Fallback to string representation if structured extraction fails
                    structured_groups.append(
                        {
                            "fgs": [],
                            "psmis": "",
                            "fg_mols": None,
                            "description": f"{str(e)}",
                        }
                    )
            return structured_groups
        else:
            return [{"None": "No fragments found"}]
    else:
        return [{"Error": "Check input SMILES"}]


def get_ertl_functional_groups_efgs_depiction(molecule: any) -> str:
    """
    This function takes an organic molecule as input and uses the algorithm
    proposed by Peter Ertl to identify functional groups within the molecule.
    The implementation by Gonzalo Colmenarejo (EFGs) is used here
    (https://github.com/rdkit/rdkit/tree/master/Contrib/efgs).
    This function returns a depiction of the molecule where the identified functional groups are highlighted.

    Parameters:
    - molecule (any): A molecule represented as an RDKit Mol object.

    Returns:
    - str: The PNG string of the molecule depiction with highlighted functional groups.
    # TODO: check what happens if no FGs are found!

    References:
    - Ertl, P. An algorithm to identify functional groups in organic molecules. J Cheminform 9, 36 (2017). https://doi.org/10.1186/s13321-017-0225-z
    - Colmenarejo, G. EFGs: A Complete and Accurate Implementation of Ertl's Functional Group Detection Algorithm in RDKit. J Chem Inf Model 65, 3 (2025). https://doi.org/10.1021/acs.jcim.4c02268
    """
    if molecule:
        # we ignore all the other outputs here
        img_text, fgs, psmis, fg_mols = efgs.get_dec_fgs(molecule)
        if img_text:
            return img_text
        else:
            return "None: No depiction generated"
    else:
        return "Error: Check input SMILES"


def get_standardized_tautomer(
    molecule: any,
    isomeric: bool = True,
) -> str:
    """Generate the standardized tautomer SMILES for a given molecule.

    Args:
        molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.
        isomeric (bool, optional): Flag to generate isomeric SMILES. Defaults to True.

    Returns:
        str: The standardized tautomer SMILES, or an error message.
    """

    if molecule:
        [a.SetAtomMapNum(0) for i, a in enumerate(molecule.GetAtoms())]
        initial_smiles = Chem.MolToSmiles(
            molecule, isomericSmiles=isomeric, kekuleSmiles=True
        )
        canonical_mol = Chem.MolFromSmiles(Chem.CanonSmiles(initial_smiles))

        if canonical_mol:
            te = TautomerEnumerator()
            standardized_mol = te.Canonicalize(canonical_mol)
            new_smiles = Chem.MolToSmiles(
                standardized_mol, isomericSmiles=isomeric, kekuleSmiles=True
            )

            return new_smiles
    else:
        return "Error Check input SMILES"


def has_cis_trans_stereochemistry(molecule: any) -> bool:
    """
    Detect whether a molecule has cis/trans (E/Z) stereochemistry assigned.

    Parameters:
    -----------
    molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.

    Returns:
    --------
    bool
        True if cis/trans stereochemistry is assigned, False otherwise
    """
    if molecule is None:
        return False

    # Check each bond for stereochemistry
    for bond in molecule.GetBonds():
        # Check if bond is a double bond
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Check if stereochemistry is assigned
            stereo = bond.GetStereo()
            if stereo in [
                Chem.BondStereo.STEREOE,
                Chem.BondStereo.STEREOZ,
                Chem.BondStereo.STEREOTRANS,
                Chem.BondStereo.STEREOCIS,
            ]:
                return True

    return False
