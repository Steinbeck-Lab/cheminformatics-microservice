from __future__ import annotations

from typing import Dict
from typing import Union

from app.modules.all_descriptors import get_cdk_rdkit_combined_descriptors
from app.modules.npscorer import get_np_score
from app.modules.toolkits.cdk_wrapper import get_CDK_descriptors
from app.modules.toolkits.cdk_wrapper import get_CDK_MolecularFormula
from app.modules.toolkits.cdk_wrapper import get_murcko_framework
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import get_rdkit_descriptors
from app.modules.tools.sugar_removal import get_sugar_info


def get_descriptors(smiles: str, toolkit: str) -> Union[tuple, str]:
    """Calculate descriptors using RDKit or CDK toolkit for the given SMILES.

    Args:
        smiles (str): SMILES input.
        toolkit (str): Toolkit choice ("rdkit" or "cdk").

    Returns:
        dict or str: Dictionary of descriptors and their values if successful,
                     or an error message if the toolkit choice is invalid or SMILES is invalid.
    """

    mol = parse_input(smiles, toolkit, False)
    if mol:
        if toolkit == "rdkit":
            Descriptors = get_rdkit_descriptors(mol)
            return Descriptors
        elif toolkit == "cdk":
            Descriptors = get_CDK_descriptors(mol)
            return Descriptors


def get_COCONUT_descriptors(smiles: str, toolkit: str) -> Union[Dict[str, float], str]:
    """Calculate COCONUT descriptors using RDKit or CDK toolkit for the given.

    SMILES.

    Args:
        smiles (str): SMILES input.
        toolkit (str): Toolkit choice ("rdkit" or "cdk").

    Returns:
        dict or str: Dictionary of COCONUT descriptors and their values if successful,
                     or an error message if the toolkit choice is invalid or SMILES is invalid.
    """
    if toolkit == "all":
        AllDescriptors = get_cdk_rdkit_combined_descriptors(smiles)
        return AllDescriptors
    else:
        Descriptors = get_descriptors(smiles, toolkit)

        rdkitMolecule = parse_input(smiles, "rdkit", False)
        nplikeliness = float(get_np_score(rdkitMolecule))

        cdkMolecule = parse_input(smiles, "cdk", False)
        hasLinearSugar, hasCircularSugars = get_sugar_info(cdkMolecule)
        framework = get_murcko_framework(cdkMolecule)
        molFormula = get_CDK_MolecularFormula(cdkMolecule)

        CombinedDescriptors = list(Descriptors)
        CombinedDescriptors.extend(
            [hasLinearSugar, hasCircularSugars, framework, nplikeliness, molFormula],
        )

        DescriptorList = (
            "atom_count",
            "heavy_atom_count",
            "molecular_weight",
            "exact_molecular_weight",
            "alogp",
            "rotatable_bond_count",
            "topological_polar_surface_area",
            "hydrogen_bond_acceptors",
            "hydrogen_bond_donors",
            "hydrogen_bond_acceptors_lipinski",
            "hydrogen_bond_donors_lipinski",
            "lipinski_rule_of_five_violations",
            "aromatic_rings_count",
            "qed_drug_likeliness",
            "formal_charge",
            "fractioncsp3",
            "number_of_minimal_rings",
            "van_der_waals_volume",
        )

    combinedDescriptors = dict(zip(DescriptorList, Descriptors))
    combinedDescriptors.update(
        {
            "linear_sugars": hasLinearSugar,
            "circular_sugars": hasCircularSugars,
            "murcko_framework": framework,
            "nplikeness": nplikeliness,
            "molecular_formula": str(molFormula),
        },
    )

    return (
        combinedDescriptors
        if len(DescriptorList) == len(Descriptors)
        else "Error Calculating Descriptors"
    )
