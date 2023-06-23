from fastapi import APIRouter
from app.modules.tools.surge import generateStructures
from app.modules.tools.sugarremoval import (
    getSugarInfo,
    removeLinearSugar,
    removeCircularSugar,
    removeLinearandCircularSugar,
)

router = APIRouter(
    prefix="/tools",
    tags=["tools"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def tools_index():
    return {"module": "tools", "message": "Successful", "status": 200}


@router.get("/generate-structures")
async def Generate_Structures(molecular_formula: str):
    """
    Generate structures using chemical structure generator based on the canonical generation path method.

    - **Molecular Formula**: required (query)
    """
    return generateStructures(molecular_formula)


@router.get("/sugar-information")
async def getsugarinformation(smiles: str):
    """
    Get information whether a molecule has circular or linear sugars.

    - **SMILES string**: required (query)
    """
    hasLinearSugar, hasCircularSugars = getSugarInfo(smiles)
    if hasLinearSugar and hasCircularSugars:
        return "The molecule contains Linear and Circular sugars"
    if hasLinearSugar and not hasCircularSugars:
        return "The molecule contains only Linear sugar"
    if hasCircularSugars and not hasLinearSugar:
        return "The molecule contains only Circular sugar"
    else:
        return "The molecule contains no sugar"


@router.get("/remove-linear-sugar")
async def removelinearsugars(smiles: str):
    """
    Detect and remove linear sugars from a given SMILES string.

    - **SMILES string**: required (query)
    """
    return removeLinearSugar(smiles)


@router.get("/remove-circular-sugar")
async def removecircularsugars(smiles: str):
    """
    Detect and remove circular sugars from a given SMILES string.

    - **SMILES string**: required (query)
    """
    return removeCircularSugar(smiles)


@router.get("/remove-linearandcircular-sugar")
async def removelinearandcircularsugars(smiles: str):
    """
    Detect and remove linear and circular sugars from a given SMILES string.

    - **SMILES string**: required (query)
    """
    return removeLinearandCircularSugar(smiles)
