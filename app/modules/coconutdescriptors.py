from rdkitmodules import getDescriptors
from cdkmodules import getSugarInfo, getMurkoFramework


def getCOCONUTDescriptors(smiles:str):
	"""This function takes a user input as
	SMILES string and returns descriptors
	those are available in COCONUT. Uses
	RDKit and CDK at the backend.
	Args (str): SMILES input.
	Returns (dict): Decriptor list as dictionary.
	"""

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
        "FractionCSP3": fsp3,
        "Number of Minimal Rings": NumRings,
        "Linear Sugars": hasLinearSugar,
        "Circular Sugars": hasCircularSugars,
        "Murko Framework": framework,
    }
    
    return AllDescriptors