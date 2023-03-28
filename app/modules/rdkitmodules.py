from chembl_structure_pipeline import standardizer
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops


def checkSMILES(smiles: str):
    """This functions checks whether or not the SMILES
    is valid. If not, it attempts to standardize the molecule
    using the ChEMBL standardization pipeline.
    """
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
        Args : molecules rdkit.Chem.rdmol.Mol: rdkit_mol Objects
        Returns (int): A number of violations of Lipinski Rules.
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
    Args (str): SMILES string
    Returns (dict): a dictionary of calculated descriptors
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
    """Convert SMILES to Mol with 3D coordinates
    Args (str): SMILES string.
    Returns (rdkil.mol): A mol object with 3D coordinates.
    optimized with MMFF94 forcefield.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        AllChem.Compute2DCoords(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        if depict:
            return Chem.MolToMolBlock(mol)
        else:
            mol = Chem.RemoveHs(mol)
            return Chem.MolToMolBlock(mol)
    else:
        return "Error reading SMILES string, check again."
