from __future__ import annotations

from chembl_structure_pipeline import standardizer
from rdkit import Chem

from app.exception_handlers import InvalidInputException
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from app.modules.toolkits.cdk_wrapper import get_CDK_SDG_mol
from app.modules.toolkits.openbabel_wrapper import get_ob_mol


def parse_input(input: str, framework: str = "rdkit", standardize: bool = False):
    """Parse and check if the input is valid.

    Args:
        input (str): Input string.

    Returns:
        Mol or None: Valid molecule object or None if an error occurs.
            If an error occurs during SMILES parsing, an error message is returned.
    """

    # auto detect the format
    format = "SMILES"

    if format == "SMILES":
        return parse_SMILES(input, framework, standardize)


def parse_SMILES(smiles: str, framework: str = "rdkit", standardize: bool = False):
    """Check whether the input SMILES string is valid.

    If not, attempt to standardize
    the molecule using the ChEMBL standardization pipeline.

    Args:
        smiles (str): Input SMILES string.
        framework (str): Framework to use for parsing. Default is "rdkit".
        standardize (bool): Whether to standardize the molecule. Default is False.

    Returns:
        Chem.Mol or None: Valid molecule object or None if an error occurs.
            If an error occurs during SMILES parsing, an error message is returned.
    """
    try:
        smiles = smiles.replace(" ", "+")
        if framework == "rdkit":
            if "R" in smiles:
                mol = get_CDK_IAtomContainer(smiles)
                mol_str = get_CDK_SDG_mol(mol)
                mol = Chem.MolFromMolBlock(mol_str)
            else:
                mol = Chem.MolFromSmiles(smiles)
            if standardize:
                mol_block = Chem.MolToMolBlock(mol)
                standardized_mol = standardizer.standardize_molblock(mol_block)
                mol = Chem.MolFromMolBlock(standardized_mol)
        elif framework == "cdk":
            mol = get_CDK_IAtomContainer(smiles)
        elif framework == "openbabel":
            mol = get_ob_mol(smiles)
        else:
            raise ValueError(f"Invalid framework specified: {framework}")

        if mol:
            return mol
        else:
            mol = get_CDK_IAtomContainer(smiles)
            mol_str = get_CDK_SDG_mol(mol)
            return Chem.MolFromMolBlock(mol_str)
    except Exception:
        raise InvalidInputException(name="smiles", value=smiles)
