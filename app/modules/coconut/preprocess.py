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


def get_parent_smiles(molecule: Chem.Mol) -> str:
    """
    Retrieves the parent SMILES string for a given SMILES string.

    Args:
        molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.

    Returns:
        str: The parent SMILES string for the given SMILES string.

    This function uses the RDKit and ChEMBL standardizer libraries to standardize the input structure
    and retrieve the parent molecule. The parent molecule represents the core molecular structure.

    If the input SMILES string is invalid or cannot be processed, the function returns an empty string.
    """

    if molecule:
        mol_block = Chem.MolToMolBlock(molecule)
        standarized = standardizer.standardize_molblock(mol_block)
        parent, _ = standardizer.get_parent_molblock(standarized)
        parent_mol = Chem.MolFromMolBlock(parent)

        if parent_mol:
            [a.SetAtomMapNum(0) for i,a in enumerate(parent_mol.GetAtoms())]
            parent_smiles = Chem.MolToSmiles(
                parent_mol, isomericSmiles=False, kekuleSmiles=True
            )
            parent_canonical_mol = Chem.MolFromSmiles(Chem.CanonSmiles(parent_smiles))

            if parent_canonical_mol:
                new_parent_smiles = Chem.MolToSmiles(
                    parent_canonical_mol, isomericSmiles=False, kekuleSmiles=True
                )
                return new_parent_smiles

    return "Error Check input SMILES"


def get_smiles(molecule: Chem.Mol, isomeric: bool = True) -> str:
    """
    Retrieves the SMILES string (Isomeric or Canonical) for a given RDKit molecule object.

    Args:
        molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.
        isomeric (bool, optional): Whether to retrieve the Isomeric SMILES (True) or the Canonical SMILES (False).
                                   Defaults to True.

    Returns:
        str: The Isomeric or Canonical SMILES string for the given molecule.
    """
    if molecule:
        [a.SetAtomMapNum(0) for i,a in enumerate(molecule.GetAtoms())]
        initial_smiles = Chem.MolToSmiles(
            molecule, isomericSmiles=isomeric, kekuleSmiles=True
        )
        canonical_mol = Chem.MolFromSmiles(Chem.CanonSmiles(initial_smiles))

        if canonical_mol:
            new_smiles = Chem.MolToSmiles(
                canonical_mol, isomericSmiles=isomeric, kekuleSmiles=True
            )
            return new_smiles

    return "Error Check input SMILES"


def get_standardized_smiles(standardized_mol_block: str) -> str:
    """
    Get the standardized SMILES representation of a molecule.

    This function takes a standardized molecular structure represented as a MolBlock and generates the corresponding
    standardized SMILES representation.

    Args:
        standardized_mol_block (str): The standardized molecular structure in MolBlock format.

    Returns:
        str: The standardized SMILES representation of the molecule.
    """
    mol = Chem.MolFromMolBlock(standardized_mol_block)
    [a.SetAtomMapNum(0) for i,a in enumerate(mol.GetAtoms())]
    standardized_smiles = Chem.MolToSmiles(
        mol, kekuleSmiles=True
    )
    canonical_mol = Chem.MolFromSmiles(Chem.CanonSmiles(standardized_smiles))
    if canonical_mol:
        new_smiles = Chem.MolToSmiles(
            canonical_mol, isomericSmiles=True, kekuleSmiles=True
        )
        return new_smiles

    return "Error Check input SMILES"


def get_molecule_hash(molecule: Chem.Mol) -> dict:
    """Return various molecule hashes for the provided SMILES.

    Args:
        molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.

    Returns:
        dict: Dictionary containing Formula, Isomeric SMILES, and Canonical SMILES.
    """
    if molecule:
        Formula = Chem.rdMolDescriptors.CalcMolFormula(molecule)
        Isomeric_SMILES = get_smiles(molecule, isomeric=True)
        Canonical_SMILES = get_smiles(molecule, isomeric=False)
        Parent_SMILES = get_parent_smiles(molecule)
        return {
            "Formula": Formula,
            "Isomeric_SMILES": Isomeric_SMILES,
            "Canonical_SMILES": Canonical_SMILES,
            "Parent_SMILES": Parent_SMILES,
        }
    else:
        return {"Error": "Check input SMILES"}


def get_representations(molecule: Chem.Mol) -> dict:
    """Return COCONUT representations for the provided SMILES.

    Args:
        molecule (Chem.Mol): An RDKit molecule object representing the molecular structure.

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
        standardized_SMILES = get_standardized_smiles(standardized_mol_block)
        standardized_mol = parse_input(standardized_SMILES, "rdkit", False)
        standardized_representations = get_representations(standardized_mol)

        # Parent molecule
        parent_canonical_smiles = original_mol_hash["Parent_SMILES"]
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
