from __future__ import annotations

from typing import Union

from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

from app.modules.toolkits.cdk_wrapper import cdk_base
from app.modules.toolkits.cdk_wrapper import get_aromatic_ring_count
from app.modules.toolkits.cdk_wrapper import get_CDK_SDG
from app.modules.toolkits.cdk_wrapper import get_tanimoto_similarity_CDK
from app.modules.toolkits.cdk_wrapper import get_vander_waals_volume
from app.modules.toolkits.cdk_wrapper import JClass
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import check_RO5_violations
from app.modules.toolkits.rdkit_wrapper import get_tanimoto_similarity_rdkit


def get_all_rdkit_descriptors(molecule: any) -> Union[tuple, str]:
    """Calculate a selected set of molecular descriptors using RDKit.

    This function takes an input SMILES string and calculates various molecular descriptors
    using RDKit.

    Args:
        molecule (Chem.mol): RDKit molecule object.

    Returns:
        tuple: A tuple containing calculated molecular descriptors. If an error occurs during processing, an error message is returned.
    """

    if molecule:
        AtomC = rdMolDescriptors.CalcNumAtoms(molecule)
        BondC = molecule.GetNumBonds()
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
        fsp3 = "%.2f" % rdMolDescriptors.CalcFractionCSP3(molecule)
        NumRings = rdMolDescriptors.CalcNumRings(molecule)
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


def get_all_cdk_descriptors(molecule: any) -> Union[tuple, str]:
    """Calculate a set of molecular descriptors using the CDK.

    This function takes a SMILES string as input and calculates various molecular descriptors
    using the CDK. The calculated descriptors are returned as a tuple.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        tuple: A tuple containing calculated molecular descriptors. If an error occurs during processing, an error message is returned.
    """
    SDGMol = get_CDK_SDG(molecule)
    if SDGMol:
        AtomCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.AtomCountDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        BondCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.BondCountDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        HeavyAtomsC = SDGMol.getAtomCount()
        WeightDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.WeightDescriptor")()
            .calculate(SDGMol)
            .getValue()
            .toString()
        )
        TotalExactMass = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator",
        ).getTotalExactMass(SDGMol)
        ALogP = (
            JClass(cdk_base + ".qsar.descriptors.molecular.ALOGPDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        NumRotatableBonds = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.RotatableBondsCountDescriptor",
            )()
            .calculate(SDGMol)
            .getValue()
        )
        TPSADescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.TPSADescriptor")()
            .calculate(SDGMol)
            .getValue()
            .toString()
        )
        HBondAcceptorCountDescriptor = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
            )()
            .calculate(SDGMol)
            .getValue()
        )
        HBondDonorCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.HBondDonorCountDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        RuleOfFiveDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.RuleOfFiveDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        AromaticRings = get_aromatic_ring_count(SDGMol)
        QEDWeighted = None
        FormalCharge = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator",
        ).getTotalFormalCharge(SDGMol)
        FractionalCSP3Descriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.FractionalCSP3Descriptor")()
            .calculate(SDGMol)
            .getValue()
            .toString()
        )
        NumRings = (
            JClass(
                cdk_base + ".graph.Cycles",
            )
            .mcb(SDGMol)
            .numberOfCycles()
        )
        VABCVolume = get_vander_waals_volume(SDGMol)

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


def get_cdk_rdkit_combined_descriptors(
    smiles: str,
) -> Union[dict, str]:
    """Calculate a selected set of molecular descriptors using CDK and RDKit
    for a given SMILES string.

    Args:
        smiles (str): A SMILES string representing a chemical compound.

    Returns:
        Union[Dict[str, Tuple[float, float]], str]:
            - If successful, a dictionary containing calculated descriptors. Each descriptor is a key
              mapped to a tuple of its values calculated by RDKit and CDK.
            - If unsuccessful due to descriptor calculation errors, returns an error message as a string.
    """
    # Calculate RDKit and CDK descriptors
    rdkitMol = parse_input(smiles, "rdkit", False)
    rdkit_descriptors = get_all_rdkit_descriptors(rdkitMol)
    cdkMol = parse_input(smiles, "cdk", False)
    cdk_descriptors = get_all_cdk_descriptors(cdkMol)

    # List of descriptors to calculate
    all_descriptors = (
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

    if len(all_descriptors) == len(rdkit_descriptors) == len(cdk_descriptors):
        combined_dict = {
            descriptor: (rdkit_desc, cdk_desc)
            for descriptor, rdkit_desc, cdk_desc in zip(
                all_descriptors,
                rdkit_descriptors,
                cdk_descriptors,
            )
        }
        return combined_dict
    else:
        return "Error: Dictionary length is invalid"


def get_table(tanimoto_values: list) -> str:
    """Convert a list of Tanimoto similarity values into an HTML table.

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


def get_tanimoto_similarity(smileslist: str, toolkit: str = "cdk") -> list:
    """Calculate the Tanimoto similarity index between pairs of SMILES strings.

    This function takes a list of SMILES strings, splits them, and calculates
    the Tanimoto similarity index between every pair of SMILES strings.

    Args:
        smileslist (str): A comma-separated list of SMILES strings.
        toolkit (str, optional): The toolkit to use for calculating similarity.
            Can be "cdk" (Chemistry Development Kit) or "rdkit" (RDKit).Defaults to "cdk".

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
                mol1 = parse_input(smiles_list[i], "rdkit", False)
                mol2 = parse_input(smiles_list[j], "rdkit", False)
                matrix[i][j] = get_tanimoto_similarity_rdkit(mol1, mol2)
            elif toolkit == "cdk":
                mol1 = parse_input(smiles_list[i], "cdk", False)
                mol2 = parse_input(smiles_list[j], "cdk", False)
                matrix[i][j] = get_tanimoto_similarity_CDK(mol1, mol2)
            else:
                raise ValueError("Unsupported toolkit:", toolkit)

    return get_table(matrix)
