import selfies as sf

from fastapi import APIRouter
from fastapi.responses import Response
from rdkit import Chem
from typing import Optional
from STOUT import translate_forward, translate_reverse
from app.modules.toolkits.cdkmodules import (
    getCDKSDGMol,
    getCXSMILES,
    getCanonSMILES,
    getInChI,
)
from app.modules.toolkits.rdkitmodules import (
    get3Dconformers,
    get2Dmol,
    getRDKitCXSMILES,
)
from app.modules.toolkits.openbabelmodules import (
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
async def Create2D_Coordinates(smiles: str, generator: Optional[str] = "cdk"):
    """
    Generate 2D Coordinates using the CDK Structure diagram generator/Rdkit/Openbabel and return the mol block.

    - **SMILES**: required (query)
    """
    if smiles:
        if generator:
            if generator == "cdk":
                return Response(
                    content=getCDKSDGMol(smiles).replace("$$$$\n", ""),
                    media_type="text/plain",
                )
            elif generator == "rdkit":
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


@router.get("/mol3d")
async def Create3D_Coordinates(smiles: str, generator: Optional[str] = "rdkit"):
    """
    Generate a random 3D conformer from SMILES using RDKit/OpenBabel.
    CDK is not used for this purpose.

    - **SMILES**: required (query parameter)
    """
    if smiles:
        if generator:
            if generator == "rdkit":
                return Response(
                    content=get3Dconformers(smiles, depict=False),
                    media_type="text/plain",
                )
            elif generator == "openbabel":
                return Response(
                    content=getOBMol(smiles, threeD=True),
                    media_type="text/plain",
                )

    else:
        return "Error reading SMILES string check again."


@router.get("/canonicalsmiles")
async def SMILES_Canonicalise(smiles: str, generator: Optional[str] = "cdk"):
    """
    Canonicalise SMILES.

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if generator:
            if generator == "cdk":
                return str(getCanonSMILES(smiles))
            elif generator == "rdkit":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return Chem.MolToSmiles(mol, kekuleSmiles=True)
            elif generator == "openbabel":
                return getOBCanonicalSMILES(smiles)

            else:
                return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/inchi")
async def SMILES_to_InChI(smiles: str, generator: Optional[str] = "cdk"):
    """
    Convert SMILES to InChI

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")

    if smiles:
        if generator:
            if generator == "cdk":
                return str(getInChI(smiles))
            elif generator == "rdkit":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return Chem.inchi.MolToInchi(mol)
            elif generator == "openbabel":
                return getOBInChI(smiles)
        else:
            return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/inchikey")
async def SMILES_to_InChIKey(smiles: str, generator: Optional[str] = "cdk"):
    """
    Convert SMILES to InChIKey:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if generator:
            if generator == "cdk":
                return str(getInChI(smiles, InChIKey=True))
            elif generator == "rdkit":
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    return Chem.inchi.MolToInchiKey(mol)
            elif generator == "openbabel":
                return getOBInChI(smiles, InChIKey=True)
        else:
            return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/cxsmiles")
async def SMILES_to_CXSMILES(smiles: str, generator: Optional[str] = "cdk"):
    """
    Convert SMILES to CXSMILES:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if generator:
            if generator == "cdk":
                cxsmiles = getCXSMILES(smiles)
                return cxsmiles
            else:
                return getRDKitCXSMILES(smiles)
    else:
        return "Error reading SMILES string check again."


@router.get("/formats")
async def SMILES_convert_to_Formats(smiles: str, generator: Optional[str] = "cdk"):
    """
    Convert SMILES to mol block:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        if generator:
            if generator == "cdk":
                response = {}
                response["mol"] = getCDKSDGMol(smiles).replace("$$$$\n", "")
                response["canonicalsmiles"] = str(getCanonSMILES(smiles))
                response["inchi"] = str(getInChI(smiles))
                response["inchikey"] = str(getInChI(smiles, InChIKey=True))
                return response

            elif generator == "rdkit":
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
            elif generator == "openbabel":
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


@router.get("/iupac")
async def SMILES_to_IUPACname(smiles: str):
    """
    Generate IUPAC name using STOUT package:

    - **SMILES**: required (query)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        iupac = translate_forward(smiles)
        return iupac


@router.get("/smiles")
async def IUPACname_or_SELFIES_to_SMILES(
    input_text: str, representation: Optional[str] = "iupac"
):
    """
    Generate SMILES from IUPAC name or selfies:

    - **input_text**: required (query)
    - **representation**: optional
    """
    if input_text:
        if representation == "iupac":
            return translate_reverse(input_text)
        elif representation == "selfies":
            selfies_d = sf.decoder(input_text)
            return selfies_d


@router.get("/selfies")
async def encode_SELFIES(smiles: str):
    """
    Generate SELFIES from SMILES:

    - **SELFIES**: required (query)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        selfies_e = sf.encoder(smiles)
        return selfies_e
