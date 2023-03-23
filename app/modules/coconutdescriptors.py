from app.modules.rdkitmodules import getDescriptors, checkSMILES
from app.modules.cdkmodules import getSugarInfo, getMurkoFramework
from app.modules.npscorer import getNPScore


def getCOCONUTDescriptors(smiles: str):
    """This function takes a user input as
    SMILES string and returns descriptors
    those are available in COCONUT. Uses
    RDKit and CDK at the backend.
    Args (str): SMILES input.
    Returns (dict): Decriptor list as dictionary.
    """
    mol = checkSMILES(smiles)
    if mol:
        (
            AtomC,
            HeavyAtomsC,
            MolWt,
            ExactMolWt,
            ALogP,
            NumRotatableBonds,
            PSA,
            HBA,
            HBD,
            Lipinski_HBA,
            Lipinski_HBD,
            Ro5Violations,
            AromaticRings,
            QEDWeighted,
            FormalCharge,
            fsp3,
            NumRings,
        ) = getDescriptors(smiles)
        hasLinearSugar, hasCircularSugars = getSugarInfo(smiles)
        framework = getMurkoFramework(smiles)
        nplikeliness = getNPScore(smiles)

        AllDescriptors = {
            "atom_count": AtomC,
            "heavy_atom_count": HeavyAtomsC,
            "molecular_weight": MolWt,
            "exact molecular_weight": ExactMolWt,
            "alogp": ALogP,
            "rotatable_bond_count": NumRotatableBonds,
            "topological_polar_surface_area": PSA,
            "hydrogen_bond_acceptors": HBA,
            "hydrogen_bond_donors": HBD,
            "hydrogen_bond_acceptors_lipinski": Lipinski_HBA,
            "hydrogen_bond_donors_lipinski": Lipinski_HBD,
            "lipinski_rule_of_five_violations": Ro5Violations,
            "aromatic_rings_count": AromaticRings,
            "qed_drug_likeliness": QEDWeighted,
            "formal_charge": FormalCharge,
            "fractioncsp3": fsp3,
            "number_of_minimal_rings": NumRings,
            "linear_sugars": hasLinearSugar,
            "circular_sugars": hasCircularSugars,
            "murko_framework": framework,
            "nplikeliness": nplikeliness,
        }

        return AllDescriptors
    else:
        return "Error reading SMILES check again."
