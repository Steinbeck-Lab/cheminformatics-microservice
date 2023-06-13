from app.modules.toolkits.rdkitmodules import getRDKitDescriptors, checkSMILES
from app.modules.toolkits.cdkmodules import (
    getSugarInfo,
    getMurkoFramework,
    getCDKDescriptors,
)
from app.modules.alldescriptors import getCDKRDKitcombinedDescriptors
from app.modules.npscorer import getNPScore


def getDescriptors(smiles: str, toolkit: str):
    """This function takes a user input as
    SMILES string and decides whether to get
    Descriptor values from RDKit or CDK.
    Args (str): SMILES input.
    Returns (list): Decriptor list as list.
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


def getCOCONUTDescriptors(smiles: str, toolkit: str):
    """This function takes a user input as
    SMILES string and returns descriptors
    those are available in COCONUT. Uses
    RDKit and CDK at the backend.
    Args (str): SMILES input.
    Returns (dict): Decriptor list as dictionary.
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
            "linear_sugars",
            "circular_sugars",
            "murko_framework",
            "nplikeliness",
        )

        if len(DescriptorList) == len(CombinedDescriptors):
            combinedDict = {
                DescriptorList[i]: CombinedDescriptors[i]
                for i in range(len(DescriptorList))
            }
            return combinedDict
        else:
            return "Error Calculating Descriptors"
