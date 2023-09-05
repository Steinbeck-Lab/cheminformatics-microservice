from app.modules.toolkits.rdkit_wrapper import getRDKitDescriptors
from typing import Dict, Union
from app.modules.toolkits.cdk_wrapper import (
    getMurkoFramework,
    getCDKDescriptors,
)
from app.modules.all_descriptors import getCDKRDKitcombinedDescriptors
from app.modules.npscorer import getNPScore
from app.modules.tools.sugar_removal import getSugarInfo
from app.modules.toolkits.helpers import parseInput


def getDescriptors(smiles: str, toolkit: str) -> Union[tuple, str]:
    """
    Calculate descriptors using RDKit or CDK toolkit for the given SMILES.

    Args:
        smiles (str): SMILES input.
        toolkit (str): Toolkit choice ("rdkit" or "cdk").

    Returns:
        dict or str: Dictionary of descriptors and their values if successful,
                     or an error message if the toolkit choice is invalid or SMILES is invalid.
    """

    mol = parseInput(smiles, toolkit, False)
    if mol:
        if toolkit == "rdkit":
            Descriptors = getRDKitDescriptors(mol)
            return Descriptors
        elif toolkit == "cdk":
            Descriptors = getCDKDescriptors(mol)
            return Descriptors


def getCOCONUTDescriptors(smiles: str, toolkit: str) -> Union[Dict[str, float], str]:
    """
    Calculate COCONUT descriptors using RDKit or CDK toolkit for the given SMILES.

    Args:
        smiles (str): SMILES input.
        toolkit (str): Toolkit choice ("rdkit" or "cdk").

    Returns:
        dict or str: Dictionary of COCONUT descriptors and their values if successful,
                     or an error message if the toolkit choice is invalid or SMILES is invalid.
    """
    if toolkit == "all":
        AllDescriptors = getCDKRDKitcombinedDescriptors(smiles)
        return AllDescriptors
    else:
        Descriptors = getDescriptors(smiles, toolkit)

        rdkitMolecule = parseInput(smiles, "rdkit", False)
        nplikeliness = float(getNPScore(rdkitMolecule))

        cdkMolecule = parseInput(smiles, "cdk", False)
        hasLinearSugar, hasCircularSugars = getSugarInfo(cdkMolecule)
        framework = getMurkoFramework(cdkMolecule)

        CombinedDescriptors = list(Descriptors)
        CombinedDescriptors.extend(
            [hasLinearSugar, hasCircularSugars, framework, nplikeliness]
        )

        DescriptorList = (
            "atom_count",
            "heavy_atom_count",
            "molecular_weight",
            "exactmolecular_weight",
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
            "van_der_walls_volume",
        )

    combinedDescriptors = dict(zip(DescriptorList, Descriptors))
    combinedDescriptors.update(
        {
            "linear_sugars": hasLinearSugar,
            "circular_sugars": hasCircularSugars,
            "murko_framework": framework,
            "nplikeness": nplikeliness,
        }
    )

    return (
        combinedDescriptors
        if len(DescriptorList) == len(Descriptors)
        else "Error Calculating Descriptors"
    )
