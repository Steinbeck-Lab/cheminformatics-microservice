from chembl_structure_pipeline import standardizer
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops
from app.modules.toolkits.cdk_wrapper import getCDKSDGMol
from hosegen import HoseGenerator


def checkSMILES(smiles: str):
    """This functions checks whether or not the SMILES
    is valid. If not, it attempts to standardize the molecule
    using the ChEMBL standardization pipeline.

    Args (str):
        SMILES string as input.
    Returns (str):
        Molblock/ Error message.

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
        print("Error reading SMILES string, check again.")


def checkRo5Violations(mol):
    """Takes a molecule and checks whether the molecule violates
    Lipinski's Rule of Five.

    Args :
        molecules rdkit.Chem.rdmol.Mol: rdkit_mol Objects
    Returns (int):
        A number of violations of Lipinski Rules.

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


def getRDKitDescriptors(smiles: str):
    """Take an input SMILES and generate a selected set of molecular
    descriptors as a dictionary.

    Args (str):
        SMILES string
    Returns (dict):
        a dictionary of calculated descriptors

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
        )
    else:
        return "Error reading SMILES string, check again."


def get3Dconformers(smiles, depict=True):
    """Convert SMILES to Mol with 3D coordinates.

    Args (str):
        SMILES string.
    Returns (rdkil.mol):
        A mol object with 3D coordinates. Optimized with MMFF94 forcefield.

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


def getTanimotoSimilarityRDKit(smiles1, smiles2):
    """
    Take two SMILES strings and calculate
    Tanimoto similarity index using Morgan
    Fingerprints.

    Args (str,str):
        SMILES strings.
    Returns (float):
        Tanimoto similarity.
    """
    # create two example molecules
    mol1 = checkSMILES(smiles1)
    mol2 = checkSMILES(smiles2)
    if mol1 and mol2:
        # generate Morgan fingerprints for each molecule
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)

        # calculate the Tanimoto similarity between the fingerprints
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return similarity
    else:
        return "Check SMILES strings for Errors"


async def getRDKitHOSECodes(smiles: str, noOfSpheres: int):
    """
    This function takes a SMILES string as input and
    returns the calculated HOSEcodes

    Args (smiles: str, noOfSpheres: int):
        SMILES string and No of Spheres as int.
    Returns:
        HOSECodes

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


def is_valid_molecule(input_text):
    """
    This functions checks whether the input text
    is a molblock or SMILES.

    Args (str):
        SMILES string or molblock.
    Returns (str):
        SMILES/Mol flag.
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


def has_stereochemistry(smiles: str):
    """
    This function checks whether the input has stereochemistry or not.

    Args (str) :
        SMILES string.
    Returns (bool):
        True or false.

    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return False

    for atom in mol.GetAtoms():
        chiral_tag = atom.GetChiralTag()
        if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
            return True

    return False


def get2Dmol(smiles: str):
    """This function takes an input as a SMILES string and
    returns a 2D mol block.

    Args (str):
        SMILES string.
    Returns (str):
        2D Mol block.

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


def getRDKitCXSMILES(smiles: str):
    """This function takes an input as a SMILES string and
    returns a CXSMILES with coordinates.

    Args (str):
        SMILES string.
    Returns (str):
        CXSMILES with coordinates.

    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        AllChem.Compute2DCoords(mol)
        return Chem.MolToCXSmiles(mol)
    else:
        return "Error reading SMILES string, check again."
