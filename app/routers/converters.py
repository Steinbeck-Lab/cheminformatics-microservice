import selfies as sf

from fastapi import APIRouter
from fastapi.responses import Response
from rdkit import Chem
from typing import Literal
from STOUT import translate_forward, translate_reverse
from app.modules.toolkits.cdk_wrapper import (
    getCDKSDGMol,
    getCXSMILES,
    getCanonSMILES,
    getInChI,
)
from app.modules.toolkits.rdkit_wrapper import (
    get3Dconformers,
    get2Dmol,
    getRDKitCXSMILES,
)
from app.modules.toolkits.openbabel_wrapper import (
    getOBMol,
    getOBCanonicalSMILES,
    getOBInChI,
)

router = APIRouter(
    prefix="/convert",
    tags=["convert"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def converters_index():
    return {"module": "converters", "message": "Successful", "status": 200}


@router.get("/mol2D")
async def Create2D_Coordinates(
    smiles: str, toolkit: Literal["cdk", "rdkit", "openbabel"]
):
    """
    Generates 2D Coordinates using the CDK Structure diagram generator/RDKit/Open Babel and returns the mol block.

    Parameters:
    - **SMILES**: required (str): The SMILES string.
    - **toolkit** (str, optional): The toolkit to use for generating 2D coordinates.
        - Supported values: "cdk" (default), "rdkit", "openbabel".

    Returns:
    - molblock (str): The generated mol block with 2D coordinates as a plain text response.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    if smiles:
        if toolkit:
            if toolkit == "cdk":
                return Response(
                    content=getCDKSDGMol(smiles).replace("$$$$\n", ""),
                    media_type="text/plain",
                )
            elif toolkit == "rdkit":
                return Response(
                    content=get2Dmol(smiles),
                    media_type="text/plain",
                )
            else:
                return Response(
                    content=getOBMol(smiles),
                    media_type="text/plain",
                )
        else:
            return "Error reading SMILES string, check again."


@router.get("/mol3D")
async def Create3D_Coordinates(smiles: str, toolkit: Literal["rdkit", "openbabel"]):
    """
    Generates a random 3D conformer from SMILES using the specified molecule toolkit.

    Parameters:
    - **SMILES**: required (str): The SMILES representation of the molecule.
    - **toolkit**: optional (str): The molecule toolkit to use.
        - Supported values: "rdkit" (default) & "openbabel".


    Returns:
    - molblock (str): The generated mol block with 3D coordinates as a plain text response.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    if smiles:
        if toolkit:
            if toolkit == "rdkit":
                return Response(
                    content=get3Dconformers(smiles, depict=False),
                    media_type="text/plain",
                )
            elif toolkit == "openbabel":
                return Response(
                    content=getOBMol(smiles, threeD=True),
                    media_type="text/plain",
                )

    else:
        return "Error reading SMILES string check again."


@router.get("/smiles")
async def IUPACname_or_SELFIES_to_SMILES(
    input_text: str, representation: Literal["iupac", "selfies"]
):
    """
    Generate SMILES from a given IUPAC name or a SELFIES representation.

    Parameters:
    - **input_text**: required (str): The input text containing either the IUPAC name or SELFIES representation.
    - **representation**: optional (str): The representation type of the input text.
        - Supported values: "iupac" (default) & "selfies".

    Returns:
    - If representation is "iupac": The generated SMILES string corresponding to the given IUPAC name.
    - If representation is "selfies": The generated SMILES string corresponding to the given SELFIES representation.

    Notes:
    - The IUPAC name should follow the standard IUPAC naming conventions for organic compounds.
    - SELFIES (Self-Referencing Embedded Strings) is a concise yet expressive chemical string notation.

    Example Usage:
    - To generate SMILES from an IUPAC name: /smiles?input_text=benzene&representation=iupac
    - To generate SMILES from a SELFIES representation: /smiles?input_text=[C][C][C]&representation=selfies
    """
    if input_text:
        if representation == "iupac":
            return translate_reverse(input_text)
        elif representation == "selfies":
            selfies_d = sf.decoder(input_text)
            return selfies_d


@router.get("/canonicalsmiles")
async def SMILES_Canonicalise(
    smiles: str, toolkit: Literal["cdk", "rdkit", "openbabel"]
):
    """
    Canonicalizes a given SMILES string according to the allowed toolkits.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to be canonicalized.
    - **toolkit**: optional (str): The toolkit to use for canonicalization.
        - Supported values: "cdk" (default), "rdkit" & "openbabel".

    Returns:
    - SMILES (str): The canonicalized SMILES string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.

    """

    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if toolkit:
            if toolkit == "cdk":
                return str(getCanonSMILES(smiles))
            elif toolkit == "rdkit":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return Chem.MolToSmiles(mol, kekuleSmiles=True)
            elif toolkit == "openbabel":
                return getOBCanonicalSMILES(smiles)

            else:
                return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/cxsmiles")
async def SMILES_to_CXSMILES(smiles: str, toolkit: Literal["cdk", "rdkit"]):
    """
    Convert SMILES to CXSMILES. For more informations:
    - https://docs.chemaxon.com/display/docs/chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.md

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default) & "rdkit".

    Returns:
    - CXSMILES (str): The converted CXSMILES string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.

    Note:
    - CXSMILES is a Chemaxon Extended SMILES which is used for storing special features of the molecules after the SMILES string.
    """

    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if toolkit:
            if toolkit == "cdk":
                cxsmiles = getCXSMILES(smiles)
                return cxsmiles
            else:
                return getRDKitCXSMILES(smiles)
    else:
        return "Error reading SMILES string check again."


@router.get("/inchi")
async def SMILES_to_InChI(smiles: str, toolkit: Literal["cdk", "rdkit", "openbabel"]):
    """
    Convert SMILES to InChI.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default), "openbabel" & "rdkit".

    Returns:
    - InChI (str): The resulting InChI string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.

    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")

    if smiles:
        if toolkit:
            if toolkit == "cdk":
                return str(getInChI(smiles))
            elif toolkit == "rdkit":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return Chem.inchi.MolToInchi(mol)
            elif toolkit == "openbabel":
                return getOBInChI(smiles)
        else:
            return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/inchikey")
async def SMILES_to_InChIKey(
    smiles: str, toolkit: Literal["cdk", "rdkit", "openbabel"]
):
    """
    Convert SMILES to InChI-Key.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default), "openbabel" & "rdkit".

    Returns:
    - InChI-Key (str): The resulting InChI-Key string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.

    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if toolkit:
            if toolkit == "cdk":
                return str(getInChI(smiles, InChIKey=True))
            elif toolkit == "rdkit":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return Chem.inchi.MolToInchiKey(mol)
            elif toolkit == "openbabel":
                return getOBInChI(smiles, InChIKey=True)
        else:
            return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/iupac")
async def SMILES_to_IUPACname(smiles: str):
    """
    Generates IUPAC name using STOUT package. For more information:
    - Rajan, K., Zielesny, A. & Steinbeck, C. STOUT: SMILES to IUPAC names using neural machine translation. J Cheminform 13, 34 (2021). https://doi.org/10.1186/s13321-021-00512-4

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.

    Returns:
    - IUPAC name (str): The resulting IUPAC name of the chemical compound.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.

    Note:
    - Here we are using STOUT v2.0 which is available at: https://github.com/Kohulan/Smiles-TO-iUpac-Translator

    Disclaimer:
    - Since STOUT is a deep learning model it does halucinate or may provide incorrect IUPAC names at times.

    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        iupac = translate_forward(smiles)
        return iupac


@router.get("/selfies")
async def encode_SELFIES(smiles: str):
    """
    Generates SELFIES string for a given SMILES string. For more information:
    - Krenn et al, SELFIES and the future of molecular string representations, Patterns, https://doi.org/10.1016/j.patter.2022.100588.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.

    Returns:
    - SELFIES (str): The resulting SELFIES of the chemical compound.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        selfies_e = sf.encoder(smiles)
        return selfies_e


@router.get("/formats")
async def SMILES_convert_to_Formats(
    smiles: str, toolkit: Literal["cdk", "rdkit", "openbabel"]
):
    """
    Convert SMILES to various molecular formats using different toolkits.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default), "openbabel" & "rdkit".

    Returns:
    - dict: A dictionary containing the converted data in various formats. The dictionary has the following keys:
        - "mol" (str): The generated 2D mol block of the molecule.
        - "canonicalsmiles" (str): The canonical SMILES representation of the molecule.
        - "inchi" (str): The InChI representation of the molecule.
        - "inchikey" (str): The InChIKey representation of the molecule.

    Note:
        - The returned dictionary may contain empty strings if conversion fails or the input SMILES string is invalid.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if toolkit:
            if toolkit == "cdk":
                response = {}
                response["mol"] = getCDKSDGMol(smiles).replace("$$$$\n", "")
                response["canonicalsmiles"] = str(getCanonSMILES(smiles))
                response["inchi"] = str(getInChI(smiles))
                response["inchikey"] = str(getInChI(smiles, InChIKey=True))
                return response

            elif toolkit == "rdkit":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    response = {}
                    response["mol"] = Chem.MolToMolBlock(mol)
                    response["canonicalsmiles"] = Chem.MolToSmiles(
                        mol, kekuleSmiles=True
                    )
                    response["inchi"] = Chem.inchi.MolToInchi(mol)
                    response["inchikey"] = Chem.inchi.MolToInchiKey(mol)
                    return response
            elif toolkit == "openbabel":
                response = {}
                response["mol"] = getOBMol(smiles)
                response["canonicalsmiles"] = getOBCanonicalSMILES(smiles)
                response["inchi"] = getOBInChI(smiles)
                response["inchikey"] = getOBInChI(smiles, InChIKey=True)
                return response
            else:
                return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."
