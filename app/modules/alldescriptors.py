from rdkit.Chem import Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops
from app.modules.toolkits.rdkitmodules import (
    checkRo5Violations,
    checkSMILES,
    getTanimotoSimilarityRDKit,
)
from app.modules.toolkits.cdk import (
    getCDKSDG,
    JClass,
    cdk_base,
    getAromaticRingCount,
    getTanimotoSimilarityCDK,
)


def getAllRDKitDescriptors(smiles: str):
    """Take an input SMILES and generate a selected set of molecular
    descriptors generated using RDKit as a list.
    Args (str): SMILES string
    Returns (list): a list of calculated descriptors
    """
    mol = checkSMILES(smiles)
    if mol:
        AtomC = rdMolDescriptors.CalcNumAtoms(mol)
        BondC = mol.GetNumBonds()
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
        fsp3 = "%.2f" % rdMolDescriptors.CalcFractionCSP3(mol)
        NumRings = rdMolDescriptors.CalcNumRings(mol)
        return (
            AtomC,
            BondC,
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


def getAllCDKDescriptors(smiles: str):
    """Take an input SMILES and generate a selected set of molecular
    descriptors generated using CDK as a list.
    Args (str): SMILES string
    Returns (list): a list of calculated descriptors
    """
    Mol = getCDKSDG(smiles)
    if Mol:
        AtomCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.AtomCountDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        BondCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.BondCountDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        HeavyAtomsC = Mol.getAtomCount()
        WeightDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.WeightDescriptor")()
            .calculate(Mol)
            .getValue()
            .toString()
        )
        TotalExactMass = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator"
        ).getTotalExactMass(Mol)
        ALogP = (
            JClass(cdk_base + ".qsar.descriptors.molecular.ALOGPDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        NumRotatableBonds = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.RotatableBondsCountDescriptor"
            )()
            .calculate(Mol)
            .getValue()
        )
        TPSADescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.TPSADescriptor")()
            .calculate(Mol)
            .getValue()
            .toString()
        )
        HBondAcceptorCountDescriptor = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.HBondAcceptorCountDescriptor"
            )()
            .calculate(Mol)
            .getValue()
        )
        HBondDonorCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.HBondDonorCountDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        RuleOfFiveDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.RuleOfFiveDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        AromaticRings = getAromaticRingCount(Mol)
        QEDWeighted = None
        FormalCharge = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator"
        ).getTotalFormalCharge(Mol)
        FractionalCSP3Descriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.FractionalCSP3Descriptor")()
            .calculate(Mol)
            .getValue()
            .toString()
        )
        NumRings = JClass(cdk_base + ".graph.Cycles").mcb(Mol).numberOfCycles()

        return (
            int(str(AtomCountDescriptor)),
            int(str(BondCountDescriptor)),
            HeavyAtomsC,
            float("{:.2f}".format(float(str(WeightDescriptor)))),
            float("{:.5f}".format(float(str(TotalExactMass)))),
            float("{:.2f}".format(float(str(ALogP).split(",")[0]))),
            int(str(NumRotatableBonds)),
            float("{:.2f}".format(float(str(TPSADescriptor)))),
            int(str(HBondAcceptorCountDescriptor)),
            int(str(HBondDonorCountDescriptor)),
            int(str(HBondAcceptorCountDescriptor)),
            int(str(HBondDonorCountDescriptor)),
            int(str(RuleOfFiveDescriptor)),
            AromaticRings,
            str(QEDWeighted),
            FormalCharge,
            float("{:.2f}".format(float(str(FractionalCSP3Descriptor)))),
            NumRings,
        )


def getCDKRDKitcombinedDescriptors(smiles: str):
    """Take an input SMILES and generates a selected set of molecular
    descriptors using CDK and RDKit and returns as a dictionary.
    Args (str): SMILES string.
    Returns (dist): a dictionary of calculated descriptors.
    """
    RDKitDescriptors = getAllRDKitDescriptors(smiles)
    CDKDescriptors = getAllCDKDescriptors(smiles)
    AllDescriptors = (
        "Atom count",
        "Bond count",
        "Heavy atom count",
        "Molecular weight",
        "Exact molecular weight",
        "Calculated LogP",
        "Rotatable bond count",
        "Topological polar surface area",
        "Hydrogen bond acceptors",
        "Hydrogen bond donors",
        "Hydrogen bond acceptors(Lipinski)",
        "Hydrogen bond donors(Lipinski)",
        "Lipinski's rule of five violations",
        "Aromatic rings count",
        "QED drug likeliness",
        "Formal Charge",
        "FractionCSP3",
        "Number of Minimal Rings",
    )

    if len(AllDescriptors) == len(RDKitDescriptors) == len(CDKDescriptors):
        combinedDict = {
            AllDescriptors[i]: (RDKitDescriptors[i], CDKDescriptors[i])
            for i in range(len(AllDescriptors))
        }
        return combinedDict
    else:
        return "Error dictionary lenth invalid"


def getTable(tanimoto_values: list):
    """
    This function converts a list into
    a html table.
    Args (list): Tanimoto similarity values.
    Returns (html): HTML table.
    """
    table_html = "<table>"
    # Add header row with column indexes
    table_html += "<tr><th></th>"
    for j in range(len(tanimoto_values[0])):
        table_html += f"<th>{j}</th>"
    table_html += "</tr>"
    # Add data rows with row indexes
    for i, row in enumerate(tanimoto_values):
        table_html += "<tr>"
        table_html += f"<td>{i}</td>"
        for cell in row:
            table_html += f"<td>{cell}</td>"
        table_html += "</tr>"
    table_html += "</table>"
    return table_html


def getTanimotoSimilarity(smileslist: str, toolkit: str = "cdk"):
    """
    Take a list of SMILES strings split them
    and calculate Tanimoto similarity index
    between every SMILES string pair and return a matrix.
    Args (str): SMILES strings.
    Returns (list): Tanimoto similarity matrix.
    """
    # Parse comma-separated list of SMILES into a list of SMILES strings
    smiles_list = smileslist.split(",")

    # Create empty matrix
    matrix = [[0] * len(smiles_list) for i in range(len(smiles_list))]

    # Multiply each value with each other and store in matrix
    for i in range(len(smiles_list)):
        for j in range(len(smiles_list)):
            if toolkit == "rdkit":
                matrix[i][j] = getTanimotoSimilarityRDKit(
                    smiles_list[i], smiles_list[j]
                )
            else:
                matrix[i][j] = getTanimotoSimilarityCDK(smiles_list[i], smiles_list[j])
    return getTable(matrix)
