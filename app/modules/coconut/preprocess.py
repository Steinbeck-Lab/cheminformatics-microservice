from rdkit import Chem
import app.modules.toolkits.rdkit_wrapper as rdkitmodules
from chembl_structure_pipeline import standardizer
import app.modules.toolkits.cdk_wrapper as cdk
from app.modules.coconut.descriptors import getCOCONUTDescriptors
from app.modules.toolkits.helpers import parseInput


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
        molecule = parseInput(input_text, "cdk", False)
        mol_block = cdk.getCDKSDGMol(molecule, V3000=False).replace("$$$$\n", "")
        return mol_block
    elif check == "mol":
        return input_text
    else:
        return "Error!, Check the input text."


def getMoleculeHash(molecule: any) -> dict:
    """
    Return various molecule hashes for the provided SMILES.

    Args:
        smiles (str): Standardized SMILES string.

    Returns:
        dict: Dictionary containing Formula, Isomeric SMILES, and Canonical SMILES.

    """
    if molecule:
        Formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)
        Isomeric_SMILES = Chem.MolToSmiles(molecule, kekuleSmiles=True)
        Canonical_SMILES = Chem.MolToSmiles(
            molecule, kekuleSmiles=True, isomericSmiles=False
        )
        return {
            "Formula": Formula,
            "Isomeric_SMILES": Isomeric_SMILES,
            "Canonical_SMILES": Canonical_SMILES,
        }
    else:
        return {"Error": "Check input SMILES"}


def getRepresentations(molecule: any) -> dict:
    """
    Return COCONUT representations for the provided SMILES.

    Args:
        smiles (str): SMILES string.

    Returns:
        dict: Dictionary containing InChI, InChi Key, and Murko framework.

    """
    if molecule:
        InChI = Chem.inchi.MolToInchi(molecule)
        InChI_Key = Chem.inchi.MolToInchiKey(molecule)
        cdkMolecule = parseInput(Chem.MolToSmiles(molecule), "cdk", False)
        Murko = cdk.getMurkoFramework(cdkMolecule)
        return {"InChI": InChI, "InChI_Key": InChI_Key, "Murko": Murko}
    else:
        return {"Error": "Check input SMILES"}


def getCOCONUTpreprocessing(input_text: str) -> dict:
    """
    Preprocess user input text suitable for the COCONUT database submission data.

    Args:
        input_text (str): Input text (Mol/str).

    Returns:
        dict: COCONUT preprocessed data.

    """
    original_mol = getMolBlock(input_text)
    standarised_mol_block = standardizer.standardize_molblock(original_mol)
    standardised_SMILES = Chem.MolToSmiles(
        Chem.MolFromMolBlock(standarised_mol_block), kekuleSmiles=True
    )

    rdkitMol = parseInput(standardised_SMILES, "rdkit", False)
    molecule_hash = getMoleculeHash(rdkitMol)

    parent_canonical_smiles = molecule_hash["Canonical_SMILES"]
    cdkParentMol = parseInput(parent_canonical_smiles, "cdk", False)
    parent_2D_molblock = cdk.getCDKSDGMol(cdkParentMol, V3000=False).replace(
        "$$$$\n", ""
    )
    parent_2D_molblock_v3 = cdk.getCDKSDGMol(cdkParentMol, V3000=True).replace(
        "$$$$\n", ""
    )
    rdkitParentMol = parseInput(parent_canonical_smiles, "rdkit", False)
    parent_3D_molblock = rdkitmodules.get3Dconformers(rdkitParentMol)

    parent_representations = getRepresentations(rdkitParentMol)
    parent_descriptors = getCOCONUTDescriptors(parent_canonical_smiles, "rdkit")

    if rdkitmodules.has_stereochemistry(rdkitMol):
        variant_isomeric_smiles = molecule_hash["Isomeric_SMILES"]
        cdkVariantMol = parseInput(variant_isomeric_smiles, "cdk", False)
        variant_2D_molblock = cdk.getCDKSDGMol(cdkVariantMol, V3000=False).replace(
            "$$$$\n", ""
        )
        variant_2D_molblock_v3 = cdk.getCDKSDGMol(cdkVariantMol, V3000=True).replace(
            "$$$$\n", ""
        )
        rdkitVariantMol = parseInput(standardised_SMILES, "rdkit", False)
        variant_3D_molblock = rdkitmodules.get3Dconformers(rdkitVariantMol)
        variant_representations = getRepresentations(rdkitVariantMol)
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
