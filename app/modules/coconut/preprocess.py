from __future__ import annotations

from chembl_structure_pipeline import standardizer
from rdkit import Chem

import app.modules.toolkits.cdk_wrapper as cdk
import app.modules.toolkits.rdkit_wrapper as rdkitmodules
from app.modules.coconut.descriptors import get_COCONUT_descriptors
from app.modules.toolkits.helpers import parse_input


def get_mol_block(input_text: str) -> str:
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
        molecule = parse_input(input_text, "cdk", False)
        mol_block = cdk.get_CDK_SDG_mol(
            molecule,
            V3000=False,
        ).replace("$$$$\n", "")
        return mol_block
    elif check == "mol":
        return input_text
    else:
        return "Error!, Check the input text."


def get_molecule_hash(molecule: any) -> dict:
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
            molecule,
            kekuleSmiles=True,
            isomericSmiles=False,
        )
        return {
            "Formula": Formula,
            "Isomeric_SMILES": Isomeric_SMILES,
            "Canonical_SMILES": Canonical_SMILES,
        }
    else:
        return {"Error": "Check input SMILES"}


def get_representations(molecule: any) -> dict:
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
        cdkMolecule = parse_input(Chem.MolToSmiles(molecule), "cdk", False)
        Murko = cdk.get_murko_framework(cdkMolecule)
        return {"InChI": InChI, "InChI_Key": InChI_Key, "Murko": Murko}
    else:
        return {"Error": "Check input SMILES"}


def get_COCONUT_preprocessing(input_text: str) -> dict:
    """
    Preprocess user input text suitable for the COCONUT database submission data.

    Args:
        input_text (str): Input text (Mol/str).

    Returns:
        dict: COCONUT preprocessed data.
    """
    original_mol = get_mol_block(input_text)
    standarised_mol_block = standardizer.standardize_molblock(original_mol)
    standardised_SMILES = Chem.MolToSmiles(
        Chem.MolFromMolBlock(standarised_mol_block),
        kekuleSmiles=True,
    )

    rdkitMol = parse_input(standardised_SMILES, "rdkit", False)
    molecule_hash = get_molecule_hash(rdkitMol)

    parent_canonical_smiles = molecule_hash["Canonical_SMILES"]
    cdkParentMol = parse_input(parent_canonical_smiles, "cdk", False)
    parent_2D_molblock = cdk.get_CDK_SDG_mol(cdkParentMol, V3000=False).replace(
        "$$$$\n",
        "",
    )
    parent_2D_molblock_v3 = cdk.get_CDK_SDG_mol(cdkParentMol, V3000=True).replace(
        "$$$$\n",
        "",
    )
    rdkitParentMol = parse_input(parent_canonical_smiles, "rdkit", False)
    parent_3D_molblock = rdkitmodules.get_3d_conformers(rdkitParentMol)

    parent_representations = get_representations(rdkitParentMol)
    parent_descriptors = get_COCONUT_descriptors(
        parent_canonical_smiles,
        "rdkit",
    )

    if rdkitmodules.has_stereochemistry(rdkitMol):
        variant_isomeric_smiles = molecule_hash["Isomeric_SMILES"]
        cdkVariantMol = parse_input(variant_isomeric_smiles, "cdk", False)
        variant_2D_molblock = cdk.get_CDK_SDG_mol(cdkVariantMol, V3000=False).replace(
            "$$$$\n",
            "",
        )
        variant_2D_molblock_v3 = cdk.get_CDK_SDG_mol(cdkVariantMol, V3000=True).replace(
            "$$$$\n",
            "",
        )
        rdkitVariantMol = parse_input(standardised_SMILES, "rdkit", False)
        variant_3D_molblock = rdkitmodules.get_3d_conformers(rdkitVariantMol)
        variant_representations = get_representations(rdkitVariantMol)
        variant_descriptors = get_COCONUT_descriptors(
            variant_isomeric_smiles,
            "rdkit",
        )

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
