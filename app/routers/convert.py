import selfies as sf

from fastapi import APIRouter
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Optional
from STOUT import translate_forward, translate_reverse
from app.modules.cdkmodules import getCDKSDGMol

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
async def smiles_mol(smiles: str, generator: Optional[str] = "cdk"):
    """
    Convert SMILES to mol block:

    - **smiles**: required (query parameter)
    - **generator**: optional (defaults: cdk)
    """
    if smiles:
        if generator:
            if generator == "cdk":
                return getCDKSDGMol(smiles)
            else:
                m = Chem.MolFromSmiles(smiles)
                AllChem.Compute2DCoords(m)
                return Chem.MolToMolBlock(m)
    else:
        return None


@router.get("/rdkit3d")
async def smiles_generate3dconformer(smiles: str):
    """
    Generate a random 3D conformer from SMILES:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(m)
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=0xF00D)
        AllChem.MMFFOptimizeMolecule(m)
        m = Chem.RemoveHs(m)
        return Chem.MolToMolBlock(m)
    else:
        return None


@router.get("/canonicalsmiles")
async def smiles_canonicalise(smiles: str):
    """
    Cannonicalise SMILES:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(m)
    else:
        return None


@router.get("/inchi")
async def smiles_inchi(smiles: str):
    """
    Convert SMILES to InChI:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        return Chem.inchi.MolToInchi(m)
    else:
        return None


@router.get("/inchikey")
async def smiles_inchikey(smiles: str):
    """
    Convert SMILES to InChIKey:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        return Chem.inchi.MolToInchiKey(m)
    else:
        return None


@router.get("/convert")
async def smiles_convert(smiles: str):
    """
    Convert SMILES to mol block:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        response = {}
        response["mol"] = Chem.MolToMolBlock(m)
        response["cannonicalsmiles"] = Chem.MolToSmiles(m)
        response["inchi"] = Chem.inchi.MolToInchi(m)
        response["inchikey"] = Chem.inchi.MolToInchiKey(m)
        return response
    else:
        return None


@router.get("/iupac")
async def smiles_iupac(smiles: str):
    """
    Generate IUPAC name using STOUT package:

    - **smiles**: required (query)
    """
    if smiles:
        iupac = translate_forward(smiles)
        return iupac


@router.get("/smiles")
async def iupac_smiles(iupac: Optional[str], selfies: Optional[str]):
    """
    Generate SMILES from IUPAC name or selfies:

    - **iupac**: optional
    - **selfies**: optional
    """
    if iupac:
        return translate_reverse(iupac)
    elif selfies:
        selfies_d = sf.decoder(selfies)
        return selfies_d


@router.get("/selfies")
async def encodeselfies(smiles: str):
    if smiles:
        selfies_e = sf.encoder(smiles)
        return selfies_e
