from fastapi import Request, APIRouter
from typing import Optional
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
)
from chembl_structure_pipeline import standardizer
from fastapi.responses import Response, HTMLResponse
from app.modules.npscorer import getNPScore
from app.modules.classyfire import classify, result
from app.modules.cdkmodules import getCDKSDGMol
from app.modules.depict import getRDKitDepiction, getCDKDepiction
from app.modules.rdkitmodules import get3Dconformers
from app.modules.coconutdescriptors import getCOCONUTDescriptors
import pandas as pd
from fastapi.templating import Jinja2Templates

router = APIRouter(
    prefix="/chem",
    tags=["chem"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)

templates = Jinja2Templates(directory="app/templates")


@router.get("/")
async def chem_index():
    return {"module": "chem", "message": "Successful", "status": 200}


@router.get("/stereoisomers")
async def SMILES_stereoisomers(smiles: str):
    """
    Enumerate all possible stereoisomers based on the chiral centers in the given smiles:

    - **smiles**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        isomers = tuple(EnumerateStereoisomers(mol))
        smilesArray = []
        for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
            smilesArray.append(smi)
        return smilesArray
    else:
        return "Error reading SMILES string, check again."


@router.post("/standardize")
async def standardize_mol(request: Request):
    """
    Standardize molblock using the ChEMBL curation pipeline routine:

    - **mol**: required
    """
    body = await request.json()
    mol = body["mol"]
    if mol:
        standardized_mol = standardizer.standardize_molblock(mol)
        rdkit_mol = Chem.MolFromMolBlock(standardized_mol)
        smiles = Chem.MolToSmiles(rdkit_mol, kekuleSmiles=True)
        response = {}
        response["standardized_mol"] = standardized_mol
        response["cannonical_smiles"] = smiles
        response["inchi"] = Chem.inchi.MolToInchi(rdkit_mol)
        response["inchikey"] = Chem.inchi.MolToInchiKey(rdkit_mol)
        return response
    else:
        return "Error reading SMILES string, check again."


@router.get("/descriptors")
async def SMILES_descriptors(
    smiles: str, format: Optional[str] = "json", toolkit: Optional[str] = "rdkit"
):
    """
    Generate standard descriptors for the input molecules (smiles):

    - **smiles**: required (query)
    """
    if smiles:
        if format == "html":
            data = getCOCONUTDescriptors(smiles, toolkit)
            if toolkit == "all":
                headers = ["Descriptor name", "RDKit Descriptors", "CDK Descriptors"]
                df = pd.DataFrame.from_dict(
                    data, orient="index", columns=headers[1:], dtype=object
                )
                df.insert(0, headers[0], df.index)
            else:
                headers = ["Descriptor name", "Values"]
                df = pd.DataFrame.from_dict(data, orient="index", columns=headers[1:])
                df.insert(0, headers[0], df.index)
            with open("app/templates/style.css", "r") as file:
                css_style = file.read()
            html_table = df.to_html(index=False)
            return Response(content=css_style + html_table, media_type="text/html")
        else:
            return getCOCONUTDescriptors(smiles, toolkit)


@router.get("/npscore")
async def NPlikeliness_score(smiles: str):
    """
    Generate natural product likeliness score based on RDKit implementation

    - **smiles**: required (query)
    """
    if smiles:
        np_score = getNPScore(smiles)
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
async def CDK2D_coordinates(smiles: str):
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return getCDKSDGMol(smiles)
        else:
            return "Error reading SMILES string, check again."


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


@router.get("/depict3D", response_class=HTMLResponse)
async def depict3D_molecule(
    request: Request,
    smiles: str,
):
    if smiles:
        content = {"request": request, "molecule": get3Dconformers(smiles)}
        return templates.TemplateResponse("mol.html", content)


# @app.get("/molecules/", response_model=List[schemas.Molecule])
# def read_molecules(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
#     molecules = crud.get_molecules(db, skip=skip, limit=limit)
#     return molecules
