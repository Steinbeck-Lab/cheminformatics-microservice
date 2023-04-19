import selfies as sf

from fastapi import APIRouter
from fastapi.responses import Response
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional
from STOUT import translate_forward, translate_reverse
from app.modules.cdkmodules import getCDKSDGMol
from app.modules.rdkitmodules import get3Dconformers

router = APIRouter(
    prefix="/convert",
    tags=["convert"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def converters_index():
    return {"module": "converters", "message": "Successful", "status": 200}


@router.get("/mol")
async def SMILES_Mol(smiles: str, generator: Optional[str] = "cdk"):
    """
    Convert SMILES to mol block:

    - **SMILES**: required (query parameter)
    - **generator**: optional (defaults: cdk)
    """
    if smiles:
        if generator:
            if generator == "cdk":
                return Response(content=getCDKSDGMol(smiles), media_type="text/plain")
            else:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    AllChem.Compute2DCoords(mol)
                    return Response(
                        content=Chem.MolToMolBlock(mol), media_type="text/plain"
                    )
                else:
                    return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/rdkit3d")
async def SMILES_Generate3DConformer(smiles: str):
    """
    Generate a random 3D conformer from SMILES using RDKit:

    - **SMILES**: required (query parameter)
    """
    if smiles:
        return Response(
            content=get3Dconformers(smiles, depict=False), media_type="text/plain"
        )
    else:
        return "Error reading SMILES string check again."


@router.get("/canonicalsmiles")
async def SMILES_Canonicalise(smiles: str):
    """
    Cannonicalise SMILES:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol)
        else:
            return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/inchi")
async def SMILES_to_InChI(smiles: str):
    """
    Convert SMILES to InChI:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.inchi.MolToInchi(mol)
        else:
            return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/inchikey")
async def SMILES_to_InChIKey(smiles: str):
    """
    Convert SMILES to InChIKey:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.inchi.MolToInchiKey(mol)
        else:
            return "Error reading SMILES string check again."
    else:
        return "Error reading SMILES string check again."


@router.get("/formats")
async def SMILES_convert_to_Formats(smiles: str):
    """
    Convert SMILES to mol block:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            response = {}
            response["mol"] = Chem.MolToMolBlock(mol)
            response["cannonicalsmiles"] = Chem.MolToSmiles(mol, kekuleSmiles=True)
            response["inchi"] = Chem.inchi.MolToInchi(mol)
            response["inchikey"] = Chem.inchi.MolToInchiKey(mol)
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
async def IUPACname_or_SELFIES_to_SMILES(iupac: Optional[str], selfies: Optional[str]):
    """
    Generate SMILES from IUPAC name or selfies:

    - **IUPAC**: optional
    - **SELFIES**: optional
    """
    if iupac:
        return translate_reverse(iupac)
    elif selfies:
        selfies_d = sf.decoder(selfies)
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
