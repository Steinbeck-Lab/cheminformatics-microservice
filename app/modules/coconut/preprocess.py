from __future__ import annotations

from chembl_structure_pipeline import checker
from chembl_structure_pipeline import standardizer
from rdkit import Chem

import app.modules.toolkits.cdk_wrapper as cdk
import app.modules.toolkits.rdkit_wrapper as rdkitmodules
from app.modules.coconut.descriptors import get_COCONUT_descriptors
from app.modules.toolkits.helpers import InvalidInputException
from app.modules.toolkits.helpers import parse_input


def get_mol_block(input_text: str) -> str:
    """Generate a Molblock from input text using CDK.

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
    """Return various molecule hashes for the provided SMILES.

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
    """Return COCONUT representations for the provided SMILES.

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
        return {
            "standard_inchi": InChI,
            "standard_inchikey": InChI_Key,
            "murko_framework": Murko,
        }
    else:
        return {"Error": "Check input SMILES"}


def get_COCONUT_preprocessing(input_text: str) -> dict:
    """Preprocess user input text suitable for the COCONUT database submission.

    Args:
        input_text (str): Input text (Mol/str).

    Returns:
        dict: COCONUT preprocessed data.
    """
    try:
        original_mol = parse_input(input_text, "rdkit", False)
        original_mol_block = get_mol_block(input_text)
        original_mol_hash = get_molecule_hash(original_mol)
        original_representations = get_representations(original_mol)
        original_descriptors = get_COCONUT_descriptors(input_text, "rdkit")
        standarised_mol_block = standardizer.standardize_molblock(original_mol_block)

        standardized_SMILES = Chem.MolToSmiles(
            Chem.MolFromMolBlock(standarised_mol_block),
            kekuleSmiles=True,
        )

        standardized_mol = parse_input(standardized_SMILES, "rdkit", False)
        standardized_representations = get_representations(standardized_mol)
        standardized_descriptors = get_COCONUT_descriptors(standardized_SMILES, "rdkit")

        parent_canonical_smiles = original_mol_hash["Canonical_SMILES"]
        rdkitParentMol = parse_input(parent_canonical_smiles, "rdkit", False)
        parent_3D_molblock = rdkitmodules.get_3d_conformers(rdkitParentMol)

        parent_representations = get_representations(rdkitParentMol)
        parent_descriptors = get_COCONUT_descriptors(parent_canonical_smiles, "rdkit")

        return {
            "original": {
                "representations": {
                    "2D_MOL": original_mol_block,
                    "3D_MOL": rdkitmodules.get_3d_conformers(original_mol),
                    "cannonical_smiles": original_mol_hash["Isomeric_SMILES"],
                    **original_representations,
                },
                "has_stereo": rdkitmodules.has_stereochemistry(original_mol),
                "descriptors": original_descriptors,
                "errors": checker.check_molblock(original_mol_block),
            },
            "standardized": {
                "representations": {
                    "2D_MOL": original_mol_block,
                    "3D_MOL": rdkitmodules.get_3d_conformers(standardized_mol),
                    "cannonical_smiles": standardized_SMILES,
                    **standardized_representations,
                },
                "has_stereo": rdkitmodules.has_stereochemistry(standardized_mol),
                "descriptors": standardized_descriptors,
                "errors": checker.check_molblock(standarised_mol_block),
            },
            "parent": {
                "representations": {
                    "3D_MOL": parent_3D_molblock,
                    "cannonical_smiles": parent_canonical_smiles,
                    **parent_representations,
                },
                "has_stereo": rdkitmodules.has_stereochemistry(rdkitParentMol),
                "descriptors": parent_descriptors,
            },
        }
    except InvalidInputException as e:
        return {"Error": f"Invalid input SMILES {e}"}
