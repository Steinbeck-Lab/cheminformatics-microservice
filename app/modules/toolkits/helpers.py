from chembl_structure_pipeline import standardizer
from rdkit import Chem
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer, get_CDK_SDG_mol
from app.modules.toolkits.openbabel_wrapper import get_ob_mol
from app.exception_handlers import InvalidInputException


def parse_input(input: str, framework: str = "rdkit", standardize: bool = False):
    """
    Parse and check if the input is valid.

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
    """
    Check whether the input SMILES string is valid. If not, attempt to standardize
    the molecule using the ChEMBL standardization pipeline.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        Chem.Mol or None: Valid molecule object or None if an error occurs.
            If an error occurs during SMILES parsing, an error message is returned.
    """
    try:
        smiles = smiles.replace(" ", "+")
        if framework == "rdkit":
            if smiles.__contains__("R"):
                mol = get_CDK_IAtomContainer(smiles)
                mol_str = get_CDK_SDG_mol(mol)
                mol = Chem.MolFromMolBlock(mol_str)
            else:
                mol = Chem.MolFromSmiles(smiles)
            if standardize:
                mol = standardizer.standardize_molblock(mol)
        elif framework == "cdk":
            mol = get_CDK_IAtomContainer(smiles)
        elif framework == "openbabel":
            mol = get_ob_mol(smiles)
        if mol:
            return mol
        else:
            raise InvalidInputException(name="smiles", value=smiles)
    except Exception:
        raise InvalidInputException(name="smiles", value=smiles)
