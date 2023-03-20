from fastapi import Request, APIRouter
from typing import Optional
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
)
from chembl_structure_pipeline import standardizer
from fastapi.responses import Response, HTMLResponse
from rdkit.Chem.Scaffolds import MurckoScaffold
from app.modules.npscorer import getnp_score
from app.modules.descriptor_calculator import GetBasicDescriptors
from app.modules.classyfire import classify, result
from app.modules.cdkmodules import getCDKSDGMol
from app.modules.depict import getRDKitDepiction, getCDKDepiction
from app.modules.depict3D import get3DDepiction

router = APIRouter(
    prefix="/chem",
    tags=["chem"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def chem_index():
    return {"module": "chem", "message": "Successful", "status": 200}


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


@router.get("/npscore")
async def nplikeliness_score(smiles: str):
    """
    Generate natural product likeliness score based on RDKit implementation

    - **smiles**: required (query)
    """
    if smiles:
        np_score = getnp_score(smiles)
        return np_score


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
async def depict_molecule(
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


@router.get("/depict3D")
async def depict3D_molecule(
    smiles: str,
    width: Optional[int] = 512,
    height: Optional[int] = 512,
    style: Optional[str] = "stick",
):
    if smiles:
        content__ = get3DDepiction(smiles, [width, height], style)
        return HTMLResponse(content=content__.render(), status_code=200)


# @app.get("/molecules/", response_model=List[schemas.Molecule])
# def read_molecules(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
#     molecules = crud.get_molecules(db, skip=skip, limit=limit)
#     return molecules
