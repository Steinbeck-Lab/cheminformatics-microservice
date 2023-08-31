from rdkit.Chem import Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops
from typing import Dict, List, Union, Tuple
from app.modules.toolkits.rdkit_wrapper import (
    checkRo5Violations,
    checkSMILES,
    getTanimotoSimilarityRDKit,
)
from app.modules.toolkits.cdk_wrapper import (
    getCDKSDG,
    JClass,
    cdk_base,
    getAromaticRingCount,
    getTanimotoSimilarityCDK,
    getVanderWaalsVolume,
)


def getAllRDKitDescriptors(smiles: str) -> Union[tuple, str]:
    """
    Calculate a selected set of molecular descriptors using RDKit.
    This function takes an input SMILES string and calculates various molecular descriptors
    using RDKit.

    Args:
        smiles (str): The input SMILES string representing the molecule.

    Returns:
        tuple: A tuple containing calculated molecular descriptors. If an error occurs during processing, an error message is returned.
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
        VABCVolume = None
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
            VABCVolume,
        )
    else:
        return "Error reading SMILES string, check again."


def getAllCDKDescriptors(smiles: str) -> Union[tuple, str]:
    """
    Calculate a set of molecular descriptors using the CDK.
    This function takes a SMILES string as input and calculates various molecular descriptors
    using the CDK. The calculated descriptors are returned as a tuple.

    Args:
        smiles (str): The input SMILES string representing the molecular structure.

    Returns:
        tuple: A tuple containing calculated molecular descriptors. If an error occurs during processing, an error message is returned.
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
        VABCVolume = getVanderWaalsVolume(Mol)

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
            float(str(VABCVolume)),
        )
    else:
        return "Error reading SMILES string, check again."


def getCDKRDKitcombinedDescriptors(
    smiles: str,
) -> Union[dict, str]:
    """
    Calculate a selected set of molecular descriptors using CDK and RDKit for a given SMILES string.

    Args:
        smiles (str): A SMILES string representing a chemical compound.

    Returns:
        Union[Dict[str, Tuple[float, float]], str]:
            - If successful, a dictionary containing calculated descriptors. Each descriptor is a key
              mapped to a tuple of its values calculated by RDKit and CDK.
            - If unsuccessful due to descriptor calculation errors, returns an error message as a string.
    """
    # Calculate RDKit and CDK descriptors
    RDKitDescriptors = getAllRDKitDescriptors(smiles)
    CDKDescriptors = getAllCDKDescriptors(smiles)

    # List of descriptors to calculate
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
        "Hydrogen bond acceptors (Lipinski)",
        "Hydrogen bond donors (Lipinski)",
        "Lipinski's rule of five violations",
        "Aromatic rings count",
        "QED drug likeliness",
        "Formal Charge",
        "FractionCSP3",
        "Number of Minimal Rings",
        "Van der Waals Volume",
    )

    if len(AllDescriptors) == len(RDKitDescriptors) == len(CDKDescriptors):
        combinedDict = {
            AllDescriptors[i]: (RDKitDescriptors[i], CDKDescriptors[i])
            for i in range(len(AllDescriptors))
        }
        return combinedDict
    else:
        return "Error: Descriptor calculation failed or descriptor count mismatch."


def get_table(tanimoto_values: list) -> str:
    """
    Convert a list of Tanimoto similarity values into an HTML table.

    Args:
        tanimoto_values (list): A list of lists containing Tanimoto similarity values.

    Returns:
        str: HTML representation of the table.
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


def getTanimotoSimilarity(smileslist: str, toolkit: str = "cdk") -> list:
    """
    Calculate the Tanimoto similarity index between pairs of SMILES strings.

    This function takes a list of SMILES strings, splits them, and calculates
    the Tanimoto similarity index between every pair of SMILES strings.

    Args:
        smileslist (str): A comma-separated list of SMILES strings.
        toolkit (str, optional): The toolkit to use for calculating similarity.
            Can be "cdk" (Chemistry Development Kit) or "rdkit" (RDKit).
            Defaults to "cdk".

    Returns:
        list: A matrix containing Tanimoto similarity scores.
            Rows and columns correspond to SMILES strings in the input list.

    Raises:
        ValueError: If an unsupported toolkit is provided.
    """

    # Parse comma-separated list of SMILES into a list of SMILES strings
    smiles_list = smileslist.split(",")

    # Create an empty matrix to store similarity scores
    matrix = [[0.0] * len(smiles_list) for _ in range(len(smiles_list))]

    # Populate the matrix with similarity scores
    for i in range(len(smiles_list)):
        for j in range(len(smiles_list)):
            if toolkit == "rdkit":
                matrix[i][j] = getTanimotoSimilarityRDKit(
                    smiles_list[i], smiles_list[j]
                )
            elif toolkit == "cdk":
                matrix[i][j] = getTanimotoSimilarityCDK(smiles_list[i], smiles_list[j])
            else:
                raise ValueError("Unsupported toolkit:", toolkit)

    return matrix
