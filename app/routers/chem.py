import os
import requests
import selfies as sf
from fastapi import Request, APIRouter
from typing import Optional
from rdkit import Chem

from urllib.request import urlopen
from urllib.parse import urlsplit

# from ..database import db
# from fastapi_pagination import Page, add_pagination, paginate
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
)
from chembl_structure_pipeline import standardizer
from fastapi.responses import Response, JSONResponse
from rdkit.Chem.Scaffolds import MurckoScaffold
from STOUT import translate_forward, translate_reverse
from app.modules.npscorer import getnp_score
from app.modules.descriptor_calculator import GetBasicDescriptors
from app.modules.classyfire import classify, result
from app.modules.cdkmodules import getCDKSDGMol
from app.modules.depict import getRDKitDepiction, getCDKDepiction
from app.modules.decimermodules import getPredictedSegments

router = APIRouter(
    prefix="/chem",
    tags=["chem"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)

@router.get("/")
async def chem_index():
    return {
        "module"  : "chem",
        "message" : "Successful",
        "status"  : 200
    }


@router.get("/mol")
async def smiles_mol(smiles: str):
    """
    Convert smiles to mol block:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToMolBlock(m)
    else:
        return None
    
@router.get("/cannonicalsmiles")
async def smiles_mol(smiles: str):
    """
    Cannonicalise smiles:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(m)
    else:
        return None
    

@router.get("/inchi")
async def smiles_mol(smiles: str):
    """
    Convert smiles to InChI:

    - **smiles**: required (query parameter)
    """
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        return Chem.inchi.MolToInchi(m)
    else:
        return None
    
@router.get("/inchikey")
async def smiles_mol(smiles: str):
    """
    Convert smiles to InChIKey:

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
    Convert smiles to mol block:

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


@router.get("/stereoisomers")
async def smiles_stereoisomers(smiles: str):
    """
    Enumerate all possible stereoisomers based on the chiral centers in the given smiles:

    - **smiles**: required (query parameter)
    """
    m = Chem.MolFromSmiles(smiles)
    isomers = tuple(EnumerateStereoisomers(m))
    smilesArray = []
    for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
        smilesArray.append(smi)
    return smilesArray


@router.post("/standardise")
async def standardise_mol(request: Request):
    """
    Standardize molblock using the ChEMBL curation pipeline routine:

    - **mol**: required
    """
    body = await request.json()
    mol = body["mol"]
    if mol:
        standardised_mol = standardizer.standardize_molblock(mol)
        rdkit_mol = Chem.MolFromMolBlock(standardised_mol)
        smiles = Chem.MolToSmiles(rdkit_mol)
        response = {}
        response["standardised_mol"] = standardised_mol
        response["cannonical_smiles"] = smiles
        response["inchi"] = Chem.inchi.MolToInchi(rdkit_mol)
        response["inchikey"] = Chem.inchi.MolToInchiKey(rdkit_mol)
        core = MurckoScaffold.GetScaffoldForMol(rdkit_mol)
        response["murcko_scaffold"] = Chem.MolToSmiles(core)
        return response


@router.get("/descriptors")
async def smiles_descriptors(smiles: str):
    """
    Generate standard descriptors for the input molecules (smiles):

    - **smiles**: required (query)
    """
    if smiles:
        return GetBasicDescriptors(smiles)


@router.get("/iupac")
async def smiles_iupac(smiles: str):
    """
    Generate IUPAC name using STOUT package:

    - **smiles**: required (query)
    """
    if smiles:
        iupac = translate_forward(smiles)
        return iupac


@router.post("/smiles")
async def iupac_smiles(iupac: Optional[str], selfies: Optional[str]):
    """
    Generate smiles from IUPAC name or selfies:

    - **iupac**: optional
    - **selfies**: optional
    """
    if iupac:
        return translate_reverse(iupac)
    elif selfies:
        selfies_d = sf.decoder(selfies)
        return selfies_d


@router.get("/npscore")
async def nplikeliness_score(smiles: str):
    """
    Generate natural product likeliness score based on RDKit implementation
    
    - **smiles**: required (query)
    """
    if smiles:
        np_score = getnp_score(smiles)
        return np_score


@router.get("/selfies")
async def encodeselfies(smiles: str):
    if smiles:
        selfies_e = sf.encoder(smiles)
        return selfies_e


@router.get("/classyfire/classify")
async def classyfire_classify(smiles: str):
    if smiles:
        data = await classify(smiles)
        return data


@router.get("/classyfire/{id}/result")
async def classyfire_result(id: str):
    if id:
        data = await result(id)
        return data


@router.get("/cdk2d")
async def cdk2d_coordinates(smiles: str):
    if smiles:
        return getCDKSDGMol(smiles)


@router.get("/depict")
async def depick_molecule(
    smiles: str,
    generator: Optional[str] = "cdksdg",
    width: Optional[int] = 512,
    height: Optional[int] = 512,
    rotate: Optional[int] = 0,
):
    if generator:
        if generator == "cdksdg":
            return Response(
                content=getCDKDepiction(smiles, [width, height], rotate),
                media_type="image/svg+xml",
            )
        else:
            return Response(
                content=getRDKitDepiction(smiles, [width, height], rotate),
                media_type="image/svg+xml",
            )


@router.post("/process")
async def extract_chemicalinfo(request: Request):
    body = await request.json()
    image_path = body["path"]
    reference = body["reference"]
    split = urlsplit(image_path)
    filename = "/tmp/" + split.path.split("/")[-1]
    if "img" in body:
        imgDataURI = body["img"]
        if imgDataURI:
            response = urlopen(imgDataURI)
            with open(filename, "wb") as f:
                f.write(response.file.read())
                smiles = getPredictedSegments(filename)
                os.remove(filename)
                return JSONResponse(
                    content={"reference": reference, "smiles": smiles.split(".")}
                )
    else:
        response = requests.get(image_path)
        if response.status_code == 200:
            with open(filename, "wb") as f:
                f.write(response.content)
                smiles = getPredictedSegments(filename)
                os.remove(filename)
                return JSONResponse(
                    content={"reference": reference, "smiles": smiles.split(".")}
                )


# @app.get("/molecules/", response_model=List[schemas.Molecule])
# def read_molecules(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
#     molecules = crud.get_molecules(db, skip=skip, limit=limit)
#     return molecules
