from fastapi import Request, APIRouter, Query
from typing import Optional, Literal
from fastapi.responses import Response, HTMLResponse
from app.modules.depiction import getRDKitDepiction, getCDKDepiction
from app.modules.toolkits.rdkit_wrapper import get3Dconformers
from app.modules.toolkits.openbabel_wrapper import getOBMol
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
    toolkit: Literal["cdk", "rdkit"] = Query(
        default="rdkit", description="Cheminformatics toolkit used in the backend"
    ),
    width: Optional[int] = 512,
    height: Optional[int] = 512,
    rotate: Optional[int] = 0,
    CIP: Optional[bool] = False,
    unicolor: Optional[bool] = False,
):
    """
    Generates a 2D depiction of a molecule using CDK or RDKit with the given parameters.

    Parameters:
    - **SMILES**: required (query): The SMILES representation of the molecule. [required]
    - **toolkit**: (str, optional): The toolkit to use for the depiction. Defaults to "cdk".
        - Supported values: "cdk"/ "rdkit" (default), "cdk".
    - **width**: (int, optional): The width of the generated image in pixels. Defaults to 512.
    - **height**: (int, optional): The height of the generated image in pixels. Defaults to 512.
    - **rotate**: (int, optional): The rotation angle of the molecule in degrees. Defaults to 0.
    - CIP (bool, optional): Whether to include Cahn-Ingold-Prelog (CIP) stereochemistry information. Defaults to False.
    - unicolor (bool, optional): Whether to use a single color for the molecule. Defaults to False.

    Returns:
        Response: An HTTP response containing the generated image in SVG+xml format.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    Note:
        - The `smiles` parameter is required and must be provided as a query parameter.
        - The `toolkit` parameter determines the backend library to use for molecule depiction.
          Currently supported options are "cdksdg" (CDK with SDG) and RDKit (default).
        - The `width` and `height` parameters control the dimensions of the generated image.
        - The `rotate` parameter specifies the rotation angle of the molecule in degrees.
        - The `CIP` parameter controls whether Cahn-Ingold-Prelog (CIP) stereochemistry information should be included / not.
        - The `unicolor` parameter determines whether a single color is used for the molecule.

    """
    if toolkit:
        if toolkit == "cdk":
            return Response(
                content=getCDKDepiction(smiles, [width, height], rotate, CIP, unicolor),
                media_type="image/svg+xml",
            )
        elif toolkit == "rdkit":
            return Response(
                content=getRDKitDepiction(smiles, [width, height], rotate),
                media_type="image/svg+xml",
            )
        else:
            return "Check SMILES string or Toolkit configuration."


@router.get("/3D", response_class=HTMLResponse)
async def Depict3D_Molecule(
    request: Request,
    smiles: str,
    toolkit: Literal["rdkit", "openbabel"] = Query(
        default="rdkit", description="Cheminformatics toolkit used in the backend"
    ),
):
    """
    Generate 3D depictions of molecules using OpenBabel or RDKit.

    Parameters:
    - **SMILES**: required (str): The SMILES string representing the molecule to depict.
    - **toolkit**: optional (str): The molecule toolkit to use. Default is "rdkit".
          - Supported values: "rdkit"/ "openbabel" (default), "rdkit".
    Returns:
    - If toolkit is "openbabel", returns a TemplateResponse with the molecule depiction generated using OpenBabel.
    - If toolkit is "rdkit", returns a TemplateResponse with the 3D conformers of the molecule generated using RDKit.
    - If toolkit is neither "openbabel" nor "rdkit", returns a string indicating that the SMILES string or the toolkit configuration should be checked.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    Note:
    - The function expects a GET request to the "/3D" endpoint.
    - The generated depictions are rendered using the "mol.html" template found under templates directory.

    """
    if smiles:
        if toolkit == "openbabel":
            content = {
                "request": request,
                "molecule": getOBMol(smiles, threeD=True, depict=True),
            }
            return templates.TemplateResponse("mol.html", content)
        elif toolkit == "rdkit":
            content = {"request": request, "molecule": get3Dconformers(smiles)}
            return templates.TemplateResponse("mol.html", content)
        else:
            return "Check SMILES string or Toolkit configuration."
