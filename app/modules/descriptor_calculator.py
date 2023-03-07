from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops


def checkRo5Violations(mol):
    """Takes a molecules and checks whether the molecule violates
    the Lipinski's Rule of Five.
    Args : molecules rdkit.Chem.rdmol.Mol: rdkit_mol Objects
    Returns (int): A number of violations on Lipinski Rules.
    """
    num_of_violations = 0
    if Descriptors.MolLogP(mol) > 5:
        rule_break += 1
    if Descriptors.MolWt(mol) > 500:
        rule_break += 1
    if Lipinski.NumHAcceptors(mol) > 10:
        rule_break += 1
    if Lipinski.NumHDonors(mol) > 5:
        rule_break += 1
    return num_of_violations


def GetBasicDescriptors(smiles):
    """Take an input SMILES and generates a selected set of molecular
    descriptors as a dictionary
    Args (str): SMILES string
    Returns (dict): a dictionary of calculated descriptors
    """
    mol = Chem.MolFromSmiles(smiles)
    AtomC = rdMolDescriptors.CalcNumAtoms(mol)
    HeavyAtomsC = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    MolWt = "%.2f" % Descriptors.MolWt(mol)
    ExactMolWt = "%.2f" % Descriptors.ExactMolWt(mol)
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
    FormalCharge = "%.2f" % rdmolops.GetFormalCharge(mol)

    AllDescriptors = {
        "Atom count": AtomC,
        "Heavy atom count": HeavyAtomsC,
        "Molecular weight": MolWt,
        "Exact molecular weight": ExactMolWt,
        "ALogP": ALogP,
        "Rotatable bond count": NumRotatableBonds,
        "Topological polar surface area": PSA,
        "Hydrogen bond acceptors": HBA,
        "Hydrogen bond donors": HBD,
        "Hydrogen bond acceptors(Lipinski)": Lipinski_HBA,
        "Hydrogen bond donors(Lipinski)": Lipinski_HBD,
        "Lipinski's rule of five violations": Ro5Violations,
        "Aromatic rings count": AromaticRings,
        "QED drug likeliness": QEDWeighted,
        "Formal Charge": FormalCharge,
    }

    return AllDescriptors
