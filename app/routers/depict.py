from fastapi import Request, APIRouter
from typing import Optional
from fastapi.responses import Response, HTMLResponse
from app.modules.depiction import getRDKitDepiction, getCDKDepiction
from app.modules.toolkits.rdkitmodules import get3Dconformers
from fastapi.templating import Jinja2Templates

templates = Jinja2Templates(directory="app/templates")

router = APIRouter(
    prefix="/depict",
    tags=["depict"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def depict_index():
    return {"module": "depict", "message": "Successful", "status": 200}


@router.get("/2D")
async def Depict2D_molecule(
    smiles: str,
    generator: Optional[str] = "cdksdg",
    width: Optional[int] = 512,
    height: Optional[int] = 512,
    rotate: Optional[int] = 0,
    CIP: Optional[bool] = False,
    unicolor: Optional[bool] = False,
):
    """
    Generate 2D Depictions using CDK or RDKit using given parameters.

    - **SMILES**: required (query)
    - **generator**: optional (defaults: cdk)
    - **width**: optional (defaults: 512)
    - **height**: optional (defaults: 512)
    - **rotate**: optional (defaults: 0)
    """
    if generator:
        if generator == "cdksdg":
            return Response(
                content=getCDKDepiction(smiles, [width, height], rotate, CIP, unicolor),
                media_type="image/svg+xml",
            )
        else:
            return Response(
                content=getRDKitDepiction(smiles, [width, height], rotate),
                media_type="image/svg+xml",
            )


@router.get("/3D", response_class=HTMLResponse)
async def Depict3D_Molecule(
    request: Request,
    smiles: str,
):
    """
    Generate 3D Depictions using RDKit.

    - **SMILES**: required (query)
    """
    if smiles:
        content = {"request": request, "molecule": get3Dconformers(smiles)}
        return templates.TemplateResponse("mol.html", content)
