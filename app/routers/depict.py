from __future__ import annotations

from typing import Literal
from typing import Optional

from fastapi import FastAPI
from fastapi import APIRouter
from fastapi import HTTPException
from fastapi import Query
from fastapi import Request
from fastapi import status
from fastapi.responses import Response
from fastapi.templating import Jinja2Templates
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

from app.modules.depiction import get_cdk_depiction
from app.modules.depiction import get_rdkit_depiction
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.openbabel_wrapper import get_ob_mol
from app.modules.toolkits.rdkit_wrapper import get_3d_conformers
from app.schemas import HealthCheck
from app.schemas.depict_schema import Depict2DResponse
from app.schemas.depict_schema import Depict3DResponse
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse
from app.schemas.error import NotFoundModel

templates = Jinja2Templates(directory="app/templates")
# Create the Limiter instance
limiter = Limiter(key_func=get_remote_address)

# Initialize FastAPI app
app = FastAPI()

# Add the middleware to handle rate limit exceeded errors
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

router = APIRouter(
    prefix="/depict",
    tags=["depict"],
    dependencies=[],
    responses={
        200: {"description": "OK"},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)


@router.get("/", include_in_schema=False)
@router.get(
    "/health",
    tags=["healthcheck"],
    summary="Perform a Health Check on Depict Module",
    response_description="Return HTTP Status Code 200 (OK)",
    status_code=status.HTTP_200_OK,
    response_model=HealthCheck,
    include_in_schema=False,
)
def get_health() -> HealthCheck:
    """## Perform a Health Check.

    Endpoint to perform a health check on. This endpoint can primarily be used by Docker
    to ensure a robust container orchestration and management are in place. Other
    services that rely on the proper functioning of the API service will not deploy if this
    endpoint returns any other HTTP status code except 200 (OK).
    Returns:
        HealthCheck: Returns a JSON response with the health status
    """
    return HealthCheck(status="OK")


@router.get(
    "/2D",
    summary="Generates a 2D depiction of a molecule",
    responses={
        200: {
            "description": "Successful response",
            "model": Depict2DResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def depict_2d_molecule(
    smiles: str = Query(
        title="SMILES",
        description="SMILES string to be converted",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    toolkit: Literal["cdk", "rdkit"] = Query(
        default="rdkit",
        description="Cheminformatics toolkit used in the backend",
    ),
    width: Optional[int] = Query(
        512,
        title="Width",
        description="The width of the generated image in pixels.",
    ),
    height: Optional[int] = Query(
        512,
        title="Height",
        description="The height of the generated image in pixels.",
    ),
    rotate: Optional[int] = Query(
        0,
        title="Rotate",
        description="The rotation angle of the molecule in degrees.",
    ),
    CIP: bool = Query(
        False,
        title="CIP",
        description="Whether to include Cahn-Ingold-Prelog (CIP) stereochemistry information.",
    ),
    unicolor: bool = Query(
        False,
        title="Unicolor",
        description="Whether to use a single colour for the molecule.",
    ),
    highlight: Optional[str] = Query(
        None,
        title="Substructure",
        description="SMARTS pattern to highlight atoms/bonds.",
    ),
    atomIds: Optional[str] = Query(
        None,
        title="Atom Indices",
        description="Comma-separated list of atom indices to highlight (0-based indexing).",
    ),
):
    """Generates a 2D depiction of a molecule using CDK or RDKit with the given.

    parameters.

    Parameters:
    - **SMILES**: required (query): The SMILES representation of the molecule. [required]
    - **toolkit**: (str, optional): The toolkit to use for the depiction. Defaults to "cdk".
        - Supported values: "cdk"/ "rdkit" (default), "cdk".
    - **width**: (int, optional): The width of the generated image in pixels. Defaults to 512.
    - **height**: (int, optional): The height of the generated image in pixels. Defaults to 512.
    - **rotate**: (int, optional): The rotation angle of the molecule in degrees. Defaults to 0.
    - CIP (bool, optional): Whether to include Cahn-Ingold-Prelog (CIP) stereochemistry information. Defaults to False.
    - unicolor (bool, optional): Whether to use a single colour for the molecule. Defaults to False.
    - highlight (Optional[str], optional): SMARTS pattern to highlight atoms/bonds. Defaults to "COSN".

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
        - The `unicolor` parameter determines whether a single colour is used for the molecule.
    """
    try:
        # Parse atom indices if provided
        highlight_atoms = None
        if atomIds:
            try:
                highlight_atoms = [
                    int(x.strip()) for x in atomIds.split(",") if x.strip().isdigit()
                ]
            except ValueError:
                raise HTTPException(
                    status_code=422,
                    detail="Invalid atomIds format. Please provide comma-separated integers.",
                )

        if toolkit == "cdk":
            mol = parse_input(smiles, "cdk", False)
            depiction = get_cdk_depiction(
                mol,
                [width, height],
                rotate,
                CIP=CIP,
                unicolor=unicolor,
                highlight=highlight,
                highlight_atoms=highlight_atoms,
            )
        elif toolkit == "rdkit":
            mol = parse_input(smiles, "rdkit", False)
            depiction = get_rdkit_depiction(
                mol,
                [width, height],
                rotate,
                unicolor=unicolor,
                highlight=highlight,
                highlight_atoms=highlight_atoms,
            )
        else:
            raise HTTPException(
                status_code=422,
                detail="Error reading SMILES string, please check again.",
            )
        return Response(content=depiction, media_type="image/svg+xml")
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/3D",
    summary="Generates a 3D depiction of a molecule",
    responses={
        200: {
            "description": "Successful response",
            "model": Depict3DResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
@limiter.limit("25/minute")
async def depict_3d_molecule(
    request: Request,
    smiles: str = Query(
        title="SMILES",
        description="SMILES string to be converted",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    toolkit: Literal["rdkit", "openbabel"] = Query(
        default="openbabel",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Generate 3D depictions of molecules using OpenBabel or RDKit.

    Parameters:
    - **SMILES**: required (str): The SMILES string representing the molecule to depict.
    - **toolkit**: optional (str): The molecule toolkit to use. The default is "rdkit".
          - Supported values: "rdkit"/ "openbabel" (default), "rdkit".
    Returns:
    - If the toolkit is "openbabel", returns a TemplateResponse with the molecule depiction generated using OpenBabel.
    - If the toolkit is "rdkit", returns a TemplateResponse with the 3D conformers of the molecule generated using RDKit.
    - If the toolkit is neither "openbabel" nor "rdkit", returns a string indicating that the SMILES string or the toolkit configuration should be checked.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    Note:
    - The function expects a GET request to the "/3D" endpoint.
    - The generated depictions are rendered using the "mol.html" template found under the templates directory.
    """
    try:
        if toolkit == "openbabel":
            content = {
                "request": request,
                "molecule": get_ob_mol(smiles, threeD=True, depict=True),
            }
        elif toolkit == "rdkit":
            mol = parse_input(smiles, "rdkit", False)
            content = {"request": request, "molecule": get_3d_conformers(mol)}
        else:
            raise HTTPException(
                status_code=422,
                detail="Error reading SMILES string, please check again.",
            )
        return templates.TemplateResponse("mol.html", content)
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
