from app.modules.toolkits.rdkit_wrapper import getRDKitDescriptors, checkSMILES
from typing import Dict, Union
from app.modules.toolkits.cdk_wrapper import (
    getMurkoFramework,
    getCDKDescriptors,
)
from app.modules.alldescriptors import getCDKRDKitcombinedDescriptors
from app.modules.npscorer import getNPScore
from app.modules.tools.sugarremoval import getSugarInfo


def getDescriptors(smiles: str, toolkit: str) -> Union[Dict[str, float], str]:
    """
    Calculate descriptors using RDKit or CDK toolkit for the given SMILES.

    Args:
        smiles (str): SMILES input.
        toolkit (str): Toolkit choice ("rdkit" or "cdk").

    Returns:
        dict or str: Dictionary of descriptors and their values if successful,
                     or an error message if the toolkit choice is invalid or SMILES is invalid.
    """

    mol = checkSMILES(smiles)
    if mol:
        if toolkit == "rdkit":
            Descriptors = getRDKitDescriptors(smiles)
            return Descriptors
        elif toolkit == "cdk":
            Descriptors = getCDKDescriptors(smiles)
            return Descriptors
        else:
            return "Error calculating Descriptors"
    else:
        return "Error reading SMILES check again."


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

        hasLinearSugar, hasCircularSugars = getSugarInfo(smiles)
        framework = getMurkoFramework(smiles)
        nplikeliness = float(getNPScore(smiles))
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
            "linear_sugars",
            "circular_sugars",
            "murko_framework",
            "nplikeness",
        )

        if len(DescriptorList) == len(CombinedDescriptors):
            combinedDict = {
                DescriptorList[i]: CombinedDescriptors[i]
                for i in range(len(DescriptorList))
            }
            return combinedDict
        else:
            return "Error Calculating Descriptors"
