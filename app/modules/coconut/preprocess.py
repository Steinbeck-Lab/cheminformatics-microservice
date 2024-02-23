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

    try:
        molecule = parse_input(input_text, "cdk", False)
        mol_block = cdk.get_CDK_SDG_mol(
            molecule,
            V3000=False,
        ).replace("$$$$\n", "")
        return mol_block
    except InvalidInputException:
        raise InvalidInputException(f"Invalid input SMILES: {input_text}")


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


def get_COCONUT_preprocessing(
    input_text: str, _3d_mol: bool = False, descriptors: bool = False
) -> dict:
    """Preprocess user input text suitable for the COCONUT database submission.

    Args:
        input_text (str): The input text representing a chemical compound in Mol format.
        _3d_mol (bool, optional): Flag indicating whether to generate 3D coordinates for the molecule. Defaults to False.
        descriptors (bool, optional): Flag indicating whether to generate COCONUT descriptors for the molecule. Defaults to False.

    Returns:
        dict: A dictionary containing COCONUT preprocessed data with representations, descriptors, and errors.

    Raises:
        InvalidInputException: If the input SMILES string is invalid.
    """
    try:
        # Preprocess input text
        input_text = input_text.replace(" ", "+").replace("\\\\", "\\")

        # Original molecule
        original_mol = parse_input(input_text, "rdkit", False)

        original_mol_block = get_mol_block(input_text)
        original_mol_hash = get_molecule_hash(original_mol)
        original_representations = get_representations(original_mol)

        # Standardized molecule
        standardized_mol_block = standardizer.standardize_molblock(original_mol_block)
        standardized_SMILES = Chem.MolToSmiles(
            Chem.MolFromMolBlock(standardized_mol_block), kekuleSmiles=True
        )
        standardized_mol = parse_input(standardized_SMILES, "rdkit", False)
        standardized_representations = get_representations(standardized_mol)

        # Parent molecule
        parent_canonical_smiles = original_mol_hash["Canonical_SMILES"]
        parent_mol_block = get_mol_block(parent_canonical_smiles)
        rdkitParentMol = parse_input(parent_canonical_smiles, "rdkit", False)
        parent_representations = get_representations(rdkitParentMol)

        # Compute descriptors if requested
        if descriptors:
            original_descriptors = get_COCONUT_descriptors(input_text, "rdkit")
            standardized_descriptors = get_COCONUT_descriptors(
                standardized_SMILES, "rdkit"
            )
            parent_descriptors = get_COCONUT_descriptors(
                parent_canonical_smiles, "rdkit"
            )
        else:
            original_descriptors = {"descriptors": "Not computed, enable for computing"}
            standardized_descriptors = {
                "descriptors": "Not computed, enable for computing"
            }
            parent_descriptors = {"descriptors": "Not computed, enable for computing"}

        # Compute 3D conformers if requested
        if _3d_mol:
            original_3d_mol_block = rdkitmodules.get_3d_conformers(original_mol)
            standardized_3d_mol_block = rdkitmodules.get_3d_conformers(standardized_mol)
            parent_3D_mol_block = rdkitmodules.get_3d_conformers(rdkitParentMol)
        else:
            original_3d_mol_block = "Not computed, enable for computing"
            standardized_3d_mol_block = "Not computed, enable for computing"
            parent_3D_mol_block = "Not computed, enable for computing"

        # Construct and return the COCONUT preprocessed data
        return {
            "original": {
                "representations": {
                    "2D_MOL": original_mol_block,
                    "3D_MOL": original_3d_mol_block,
                    "canonical_smiles": original_mol_hash["Isomeric_SMILES"],
                    **original_representations,
                },
                "has_stereo": rdkitmodules.has_stereochemistry(original_mol),
                "descriptors": original_descriptors,
                "errors": checker.check_molblock(original_mol_block),
            },
            "standardized": {
                "representations": {
                    "2D_MOL": original_mol_block,
                    "3D_MOL": standardized_3d_mol_block,
                    "canonical_smiles": standardized_SMILES,
                    **standardized_representations,
                },
                "has_stereo": rdkitmodules.has_stereochemistry(standardized_mol),
                "descriptors": standardized_descriptors,
                "errors": checker.check_molblock(standardized_mol_block),
            },
            "parent": {
                "representations": {
                    "2D_MOL": parent_mol_block,
                    "3D_MOL": parent_3D_mol_block,
                    "canonical_smiles": parent_canonical_smiles,
                    **parent_representations,
                },
                "has_stereo": rdkitmodules.has_stereochemistry(rdkitParentMol),
                "descriptors": parent_descriptors,
            },
        }
    except InvalidInputException:
        raise InvalidInputException(f"Invalid input SMILES: {input_text}")
