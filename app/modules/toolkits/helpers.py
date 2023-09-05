from chembl_structure_pipeline import standardizer
from rdkit import Chem
from app.modules.toolkits.cdk_wrapper import getCDKIAtomContainer, getCDKSDGMol
from app.modules.toolkits.openbabel_wrapper import getOBMol
from app.exception_handlers import InvalidInputException


def parseInput(input: str, framework: str = "rdkit", standardize: bool = False):
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
        return parseSMILES(input, framework, standardize)


def parseSMILES(smiles: str, framework: str = "rdkit", standardize: bool = False):
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
                mol = getCDKIAtomContainer(smiles)
                mol_str = getCDKSDGMol(mol)
                mol = Chem.MolFromMolBlock(mol_str)
            else:
                mol = Chem.MolFromSmiles(smiles)
            if standardize:
                mol = standardizer.standardize_molblock(mol)
        elif framework == "cdk":
            mol = getCDKIAtomContainer(smiles)
        elif framework == "openbabel":
            mol = getOBMol(smiles)
        if mol:
            return mol
        else:
            raise InvalidInputException(name="smiles", value=smiles)
    except Exception:
        raise InvalidInputException(name="smiles", value=smiles)
