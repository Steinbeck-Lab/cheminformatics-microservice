from fastapi import APIRouter, Depends, HTTPException

from rdkit import Chem

router = APIRouter(
    prefix="/chem",
    tags=["chem"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)

@router.get("/")
async def chem_index():
    return  {"message": "Hello Chem Router!"}


@router.get("/{smiles}/mol")
async def smiles_mol(smiles: str):
    if smiles:
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToMolBlock(m)
    else:
        return None