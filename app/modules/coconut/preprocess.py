from rdkit import Chem
import app.modules.toolkits.rdkit_wrapper as rdkitmodules
import app.modules.toolkits.cdk_wrapper as cdk
from app.modules.coconut.descriptors import getCOCONUTDescriptors


def getMolBlock(input_text: str) -> str:
    """
    Generate a Molblock from input text using CDK.

    Args:
        input_text (str): Input text (Mol/SMILES).

    Returns:
        str: Molblock representation.

    Raises:
        ValueError: If input_text is not a valid Mol or SMILES.

    """
    check = rdkitmodules.is_valid_molecule(input_text)

    if check == "smiles":
        mol_block = cdk.getCDKSDGMol(input_text, V3000=False).replace("$$$$\n", "")
        return mol_block
    elif check == "mol":
        return input_text
    else:
        return "Error!, Check the input text."


def getMolculeHash(smiles: str) -> dict:
    """
    Return various molecule hashes for the provided SMILES.

    Args:
        smiles (str): Standardized SMILES string.

    Returns:
        dict: Dictionary containing Formula, Isomeric SMILES, and Canonical SMILES.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        Formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        Isomeric_SMILES = Chem.MolToSmiles(mol, kekuleSmiles=True)
        Canonical_SMILES = Chem.MolToSmiles(
            mol, kekuleSmiles=True, isomericSmiles=False
        )
        return {
            "Formula": Formula,
            "Isomeric_SMILES": Isomeric_SMILES,
            "Canonical_SMILES": Canonical_SMILES,
        }


def getRepresentations(smiles: str) -> dict:
    """
    Return COCONUT representations for the provided SMILES.

    Args:
        smiles (str): SMILES string.

    Returns:
        dict: Dictionary containing InChI, InChi Key, and Murko framework.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        InChI = Chem.inchi.MolToInchi(mol)
        InChI_Key = Chem.inchi.MolToInchiKey(mol)
        Murko = cdk.getMurkoFramework(smiles)
        return {"InChI": InChI, "InChI_Key": InChI_Key, "Murko": Murko}


def getCOCONUTpreprocessing(input_text: str) -> dict:
    """
    Preprocess user input text suitable for the COCONUT database submission data.

    Args:
        input_text (str): Input text (Mol/str).

    Returns:
        dict: COCONUT preprocessed data.

    """
    original_mol = getMolBlock(input_text)
    standarised_mol_block = rdkitmodules.standardizer.standardize_molblock(original_mol)
    standardised_SMILES = Chem.MolToSmiles(
        Chem.MolFromMolBlock(standarised_mol_block), kekuleSmiles=True
    )
    molecule_hash = getMolculeHash(standardised_SMILES)
    parent_canonical_smiles = molecule_hash["Canonical_SMILES"]
    parent_2D_molblock = cdk.getCDKSDGMol(parent_canonical_smiles, V3000=False).replace(
        "$$$$\n", ""
    )
    parent_3D_molblock = rdkitmodules.get3Dconformers(parent_canonical_smiles)
    parent_2D_molblock_v3 = cdk.getCDKSDGMol(
        parent_canonical_smiles, V3000=True
    ).replace("$$$$\n", "")
    parent_representations = getRepresentations(parent_canonical_smiles)
    parent_descriptors = getCOCONUTDescriptors(parent_canonical_smiles, "rdkit")

    if rdkitmodules.has_stereochemistry(standardised_SMILES):
        variant_isomeric_smiles = molecule_hash["Isomeric_SMILES"]
        variant_2D_molblock = cdk.getCDKSDGMol(
            variant_isomeric_smiles, V3000=False
        ).replace("$$$$\n", "")
        variant_2D_molblock_v3 = cdk.getCDKSDGMol(
            variant_isomeric_smiles, V3000=True
        ).replace("$$$$\n", "")
        variant_3D_molblock = rdkitmodules.get3Dconformers(variant_isomeric_smiles)
        variant_representations = getRepresentations(variant_isomeric_smiles)
        variant_descriptors = getCOCONUTDescriptors(variant_isomeric_smiles, "rdkit")

        return {
            "original_mol": original_mol,
            "standardised_mol": standarised_mol_block,
            "standardised_SMILES": standardised_SMILES,
            "molecule_hash": molecule_hash,
            "parent": {
                "2D_mol": parent_2D_molblock,
                "3D_mol": parent_3D_molblock,
                "v3000": parent_2D_molblock_v3,
                "representations": parent_representations,
                "descriptors": parent_descriptors,
            },
            "stereochemical_variants": True,
            "variants": {
                "2D_mol": variant_2D_molblock,
                "3D_mol": variant_3D_molblock,
                "v3000": variant_2D_molblock_v3,
                "representations": variant_representations,
                "descriptors": variant_descriptors,
            },
        }

    return {
        "original_mol": original_mol,
        "standardised_mol": standarised_mol_block,
        "standardised_SMILES": standardised_SMILES,
        "molecule_hash": molecule_hash,
        "parent": {
            "2D_mol": parent_2D_molblock,
            "3D_mol": parent_3D_molblock,
            "v3000": parent_2D_molblock_v3,
            "representations": parent_representations,
            "descriptors": parent_descriptors,
        },
        "stereochemical_variants": False,
    }
