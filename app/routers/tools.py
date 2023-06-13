from fastapi import APIRouter
from app.modules.tools.surge import generateStructures

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
    Generate structures using chemical structure generator based on the canonical generation path method

    - **Molecular Formula**: required (query)
    """
    return generateStructures(molecular_formula)
