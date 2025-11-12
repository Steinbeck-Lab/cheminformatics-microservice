from __future__ import annotations

import logging
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

from app.modules.depiction import get_rdkit_depiction
from app.modules.depiction_enhanced import get_cdk_depiction
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.openbabel_wrapper import get_ob_mol
from app.modules.toolkits.rdkit_wrapper import get_3d_conformers
from app.schemas import HealthCheck
from app.schemas.depict_schema import Depict2DResponse
from app.schemas.depict_schema import Depict3DResponse
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse
from app.schemas.error import NotFoundModel

# Configure logging
logger = logging.getLogger(__name__)

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
    showAtomNumbers: bool = Query(
        False,
        title="Show Atom Numbers",
        description="Whether to show atom numbers on the molecular depiction.",
    ),
    hydrogen_display: Literal[
        "Provided", "Minimal", "Explicit", "Stereo", "Smart"
    ] = Query(
        default="Smart",
        title="Hydrogen Display",
        description=(
            "Control how hydrogen atoms are displayed in the molecular depiction. "
            "Options: 'Provided' (as-is), 'Minimal' (suppress all H), "
            "'Explicit' (show all H), 'Stereo' (show stereo-relevant H only), "
            "'Smart' (intelligent H placement, recommended for chiral molecules)."
        ),
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
    - **CIP** (bool, optional): Whether to include Cahn-Ingold-Prelog (CIP) stereochemistry information. Defaults to False.
    - **highlight** (Optional[str], optional): SMARTS pattern to highlight atoms/bonds. Defaults to None.
    - **atomIds** (Optional[str], optional): Comma-separated atom indices to highlight. Defaults to None.
    - **showAtomNumbers** (bool, optional): Whether to show atom numbers on the molecular depiction. Defaults to False.
    - **hydrogen_display** (str, optional): Control hydrogen display mode. Defaults to "Smart".
        - "Provided": Keep hydrogens as-is (no changes)
        - "Minimal": Suppress all hydrogens (clean depiction)
        - "Explicit": Show all hydrogens explicitly
        - "Stereo": Show only stereo-relevant hydrogens (chiral centers, E/Z bonds)
        - "Smart": Intelligent hydrogen placement (recommended for chiral molecules)


    Returns:
        Response: An HTTP response containing the generated image in SVG+xml format.

    Raises:
    - HTTPException (422): If the SMILES string is invalid or parameters are incorrect.

    Note:
        - The `smiles` parameter is required and must be provided as a query parameter.
        - The `toolkit` parameter determines the backend library to use for molecule depiction.
          Currently supported options are "cdk" and "rdkit" (default).
        - The `width` and `height` parameters control the dimensions of the generated image.
        - The `rotate` parameter specifies the rotation angle of the molecule in degrees.
        - The `CIP` parameter controls whether Cahn-Ingold-Prelog (CIP) stereochemistry information should be included.
        - The `unicolor` parameter determines whether a single colour is used for the molecule.
        - The `hydrogen_display` parameter is crucial for molecules with stereochemistry:
          * Use "Smart" or "Stereo" for chiral molecules to show relevant hydrogens
          * Use "Minimal" for clean, publication-quality depictions
          * Use "Explicit" for educational materials showing all atoms
        - The `showAtomNumbers` parameter controls whether atom numbers are displayed on the depiction.
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
            from app.modules.depiction import (
                get_cdk_depiction as get_cdk_depiction_basic,
            )

            mol = parse_input(smiles, "cdk", False)
            depiction = get_cdk_depiction_basic(
                mol,
                [width, height],
                rotate,
                CIP=CIP,
                unicolor=unicolor,
                highlight=highlight,
                highlight_atoms=highlight_atoms,
                showAtomNumbers=showAtomNumbers,
                hydrogen_display=hydrogen_display,
            )
        elif toolkit == "rdkit":
            mol = parse_input(smiles, "rdkit", False)
            depiction = get_rdkit_depiction(
                mol,
                [width, height],
                rotate,
                CIP=CIP,
                unicolor=unicolor,
                highlight=highlight,
                highlight_atoms=highlight_atoms,
                showAtomNumbers=showAtomNumbers,
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


@router.get(
    "/2D_enhanced",
    summary="Generates a 2D depiction with Phase 1, 2, and 3 enhancements",
    responses={
        200: {
            "description": "Successful response - SVG image",
            "content": {"image/svg+xml": {}},
        },
    },
)
async def depict_2d_molecule_enhanced(
    smiles: str = Query(
        title="SMILES or CXSMILES",
        description="SMILES or CXSMILES string to be converted. Supports CXSMILES extensions for highlighting.",
        openapi_examples={
            "example1": {
                "summary": "Simple: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "CXSMILES: Benzene with highlighting",
                "value": "c1ccccc1 |ha:0,1,2|",
            },
            "example3": {
                "summary": "Coordination complex with dative bonds",
                "value": "[Co][N+]([O-])(=O)",
            },
            "example4": {
                "summary": "Chiral molecule with stereo hydrogens",
                "value": "C[C@H](N)C(=O)O",
            },
        },
    ),
    toolkit: Literal["cdk", "rdkit"] = Query(
        default="cdk",
        description="Cheminformatics toolkit. CDK recommended for all Phase 1-3 features.",
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
        title="CIP Stereochemistry",
        description="Whether to include Cahn-Ingold-Prelog (CIP) stereochemistry annotations (R/S, E/Z).",
    ),
    unicolor: bool = Query(
        False,
        title="Unicolor",
        description="Whether to use black and white coloring only (deprecated, use style parameter).",
    ),
    highlight: Optional[str] = Query(
        None,
        title="SMARTS Highlight",
        description="SMARTS pattern to highlight atoms/bonds.",
    ),
    atomIds: Optional[str] = Query(
        None,
        title="Atom Indices",
        description="Comma-separated list of atom indices to highlight (0-based indexing).",
    ),
    showAtomNumbers: bool = Query(
        False,
        title="Show Atom Numbers",
        description="Whether to show atom numbers on the molecular depiction.",
    ),
    hydrogen_display: Literal[
        "Provided", "Minimal", "Explicit", "Stereo", "Smart"
    ] = Query(
        default="Smart",
        title="Hydrogen Display",
        description=(
            "Control how hydrogen atoms are displayed. "
            "'Smart' recommended for chiral molecules."
        ),
    ),
    # ========== PHASE 1 PARAMETERS ==========
    abbreviate: Literal["off", "groups", "reagents", "on"] = Query(
        default="reagents",
        title="Abbreviations",
        description=(
            "Chemical abbreviation mode. "
            "'off': No abbreviations. "
            "'groups': Only functional groups (Ph, Me, Et, Boc, etc.). "
            "'reagents': Common reagents (THF, DMF, NaOH, etc.). "
            "'on': Both groups and reagents."
        ),
    ),
    dative: Literal["never", "metals", "always"] = Query(
        default="metals",
        title="Dative Bonds",
        description=(
            "Dative (coordinate) bond perception mode. "
            "'never': No dative bond perception. "
            "'metals': Only metal-ligand dative bonds (default). "
            "'always': Perceive all dative bonds including B, O."
        ),
    ),
    multicenter: Literal[
        "provided", "dative", "dashed", "dashed_neutral", "hidden", "hidden_neutral"
    ] = Query(
        default="provided",
        title="Multicenter Bonds",
        description=(
            "Multicenter bond display style for η-complexes. "
            "'provided': As-is from input. "
            "'dative': Arrow notation. "
            "'dashed': Dashed lines (preserve charges). "
            "'dashed_neutral': Dashed lines (neutralize charges). "
            "'hidden': Hide bonds (preserve charges). "
            "'hidden_neutral': Hide bonds (neutralize charges)."
        ),
    ),
    # ========== PHASE 2 PARAMETERS ==========
    annotate: Literal[
        "none", "number", "bondnumber", "mapidx", "atomvalue", "colmap", "cip"
    ] = Query(
        default="none",
        title="Annotations",
        description=(
            "Annotation mode for molecular depiction. "
            "'none': No annotations. "
            "'number': Show atom numbers (0-based). "
            "'bondnumber': Show bond numbers. "
            "'mapidx': Show atom mapping numbers (for reactions). "
            "'atomvalue': Show atom values/properties. "
            "'colmap': Color-code atoms by mapping number. "
            "'cip': Show CIP stereochemistry labels (R/S, E/Z)."
        ),
    ),
    style: Literal["cow", "cob", "cot", "bow", "bot", "wob", "nob"] = Query(
        default="cow",
        title="Style Preset",
        description=(
            "Color scheme preset. "
            "'cow': Color on White (default). "
            "'cob': Color on Black. "
            "'cot': Color on Transparent. "
            "'bow': Black on White. "
            "'bot': Black on Transparent. "
            "'wob': White on Black. "
            "'nob': Neon on Black."
        ),
    ),
    donuts: bool = Query(
        False,
        title="Aromatic Donuts",
        description="Whether to use circle-in-ring (donut) display for aromatic rings.",
    ),
    arrow: Literal["", "forward", "equ", "ngo", "ret", "res"] = Query(
        default="",
        title="Reaction Arrow",
        description=(
            "Reaction arrow type (for reactions only). "
            "'': Default forward arrow. "
            "'forward': Forward arrow (→). "
            "'equ': Equilibrium/bidirectional (⇌). "
            "'ngo': No-go/blocked (⇏). "
            "'ret': Retrosynthetic (⇒). "
            "'res': Resonance (↔)."
        ),
    ),
    alignrxnmap: bool = Query(
        True,
        title="Align Reaction Mapping",
        description="Whether to align reaction mapped atoms (for reactions only).",
    ),
    showtitle: bool = Query(
        False,
        title="Show Title",
        description="Whether to display molecule/reaction title in depiction.",
    ),
    bgcolor: Optional[str] = Query(
        None,
        title="Background Color",
        description="Custom background color as hex string (e.g., '#FFFFFF') or 'default'.",
    ),
    fgcolor: Optional[str] = Query(
        None,
        title="Foreground Color",
        description="Custom foreground/annotation color as hex string (e.g., '#000000') or 'default'.",
    ),
    # ========== PHASE 3 PARAMETERS ==========
    zoom: float = Query(
        1.3,
        title="Zoom Level",
        description="Zoom level for depiction (0.1 to 5.0). Default: 1.3.",
        ge=0.1,
        le=5.0,
    ),
    ratio: float = Query(
        1.0,
        title="Stroke Ratio",
        description="Bond thickness/stroke ratio (0.5 to 2.0). Default: 1.0.",
        ge=0.5,
        le=2.0,
    ),
    flip: bool = Query(
        False,
        title="Flip Structure",
        description="Whether to horizontally flip the molecular structure.",
    ),
    anon: bool = Query(
        False,
        title="Anonymous Display",
        description="Whether to use anonymous atom display (IUPAC recommendations).",
    ),
    smalim: int = Query(
        100,
        title="SMARTS Hit Limit",
        description="Maximum number of SMARTS pattern matches to highlight (1 to 1000).",
        ge=1,
        le=1000,
    ),
    svgunits: Literal["px", "mm", "cm", "in"] = Query(
        default="px",
        title="SVG Units",
        description="SVG coordinate units. 'px': pixels, 'mm': millimeters, 'cm': centimeters, 'in': inches.",
    ),
    perceive_radicals: bool = Query(
        False,
        title="Perceive Radicals",
        description="Whether to detect and mark unpaired electrons/radicals.",
    ),
    apply_mdl_highlighting: bool = Query(
        True,
        title="Apply MDL HILITE",
        description="Whether to apply MDL V3000 HILITE highlighting from molecule properties.",
    ),
):
    """Generate 2D molecular depictions with comprehensive Phase 1, 2, and 3 enhancements.

    This endpoint provides the most advanced molecular depiction capabilities including:

    **Phase 1 - Core Functionality:**
    - **CXSMILES**: Extended SMILES with highlighting (|ha:0,1| for atoms, |hb:0,1| for bonds)
    - **Abbreviations**: Automatic recognition of common groups (Ph, Me, Boc) and reagents (THF, DMF)
    - **Dative Bonds**: Coordinate bond detection with arrow notation (N→B, N→Metal)
    - **Multicenter Bonds**: η-complex handling for organometallic compounds (ferrocene, etc.)

    **Phase 2 - Advanced Features:**
    - **Annotations**: 7 annotation modes including CIP stereochemistry labels
    - **Style Presets**: 7 predefined color schemes for different contexts
    - **Aromatic Display**: Circle-in-ring (donut) representation for aromatic systems
    - **Reaction Arrows**: 5 different arrow types for chemical reactions
    - **Custom Colors**: Override background and foreground colors

    **Phase 3 - Professional Controls:**
    - **Radical Perception**: Detect and mark unpaired electrons
    - **MDL HILITE**: Support for MDL/SDF V3000 highlighting
    - **Zoom & Stroke**: Fine control over image scale and line thickness
    - **Structure Manipulation**: Flip structures horizontally
    - **Anonymous Display**: IUPAC-recommended atom display
    - **SMARTS Limiting**: Control number of pattern matches highlighted
    - **SVG Units**: Export in various measurement units

    **Usage Examples:**
    - Basic: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` (caffeine)
    - CXSMILES: `c1ccccc1 |ha:0,1,2|` (highlighted benzene)
    - Abbreviations: `C1CCOC1` with `abbreviate=reagents` → shows "THF"
    - Coordination: `[Co][N+]([O-])(=O)` with `dative=metals` → shows arrows
    - Styling: Use `style=cob` for color on black background
    - Annotations: Use `annotate=cip` to show R/S, E/Z labels
    - Donuts: Use `donuts=true` for aromatic circles

    **Returns**: SVG image of the molecular structure with all requested enhancements
    """
    try:
        # Parse atom indices for highlighting
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

        # Generate depiction based on toolkit
        if toolkit == "cdk":
            # CDK's SmilesParser automatically handles CXSMILES (e.g., "CCO |ha:0,1|")
            # The highlighting is extracted internally via extract_cxsmiles_highlighting()
            mol = parse_input(smiles, "cdk", False)
            depiction = get_cdk_depiction(
                mol,
                molSize=(width, height),
                rotate=rotate,
                kekulize=not donuts,  # Don't kekulize if using donuts
                CIP=CIP,
                unicolor=unicolor,
                highlight=highlight,
                highlight_atoms=highlight_atoms,
                showAtomNumbers=showAtomNumbers,
                hydrogen_display=hydrogen_display,
                # Phase 1
                abbreviate=abbreviate,
                dative=dative,
                multicenter=multicenter,
                # Phase 2
                annotate=annotate,
                style=style,
                donuts=donuts,
                arrow=arrow,
                alignrxnmap=alignrxnmap,
                showtitle=showtitle,
                bgcolor=bgcolor,
                fgcolor=fgcolor,
                # Phase 3
                zoom=zoom,
                ratio=ratio,
                flip=flip,
                anon=anon,
                smalim=smalim,
                svgunits=svgunits,
                perceive_radicals=perceive_radicals,
                apply_mdl_highlighting=apply_mdl_highlighting,
            )
        elif toolkit == "rdkit":
            # RDKit doesn't support Phase 1, 2, 3 features
            logger.warning(
                "RDKit toolkit only supports basic depiction. Use CDK for Phase 1-3 features."
            )

            if abbreviate != "off":
                logger.warning("Abbreviations not supported with RDKit toolkit")
            if dative != "never":
                logger.warning("Dative bonds not supported with RDKit toolkit")
            if multicenter != "provided":
                logger.warning("Multicenter bonds not supported with RDKit toolkit")
            if "|" in smiles:
                logger.warning(
                    "CXSMILES highlighting not fully supported with RDKit toolkit"
                )
            if annotate != "none":
                logger.warning("Annotations not supported with RDKit toolkit")
            if style != "cow":
                logger.warning("Style presets not supported with RDKit toolkit")
            if donuts:
                logger.warning("Aromatic donuts not supported with RDKit toolkit")

            mol = parse_input(smiles, "rdkit", False)
            depiction = get_rdkit_depiction(
                mol,
                [width, height],
                rotate,
                CIP=CIP,
                unicolor=unicolor,
                highlight=highlight,
                highlight_atoms=highlight_atoms,
                showAtomNumbers=showAtomNumbers,
            )
        else:
            raise HTTPException(
                status_code=422,
                detail="Invalid toolkit. Choose 'cdk' or 'rdkit'.",
            )

        return Response(content=depiction, media_type="image/svg+xml")

    except Exception as e:
        logger.error(f"Error generating enhanced depiction: {e}")
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/2D/info",
    summary="Get information about all Phase 1, 2, and 3 features",
)
async def get_features_info():
    """Get comprehensive information about all Phase 1, 2, and 3 depiction features.

    Returns detailed documentation about CXSMILES, abbreviations, dative bonds, multicenter bonds,
    annotations, style presets, aromatic display, reaction arrows, and advanced controls.
    """
    return {
        "version": "3.0.0",
        "phases": "Phase 1, 2, and 3 - Complete Feature Set",
        "features": {
            # ========== PHASE 1 ==========
            "phase1": {
                "cxsmiles": {
                    "description": "Extended SMILES with highlighting and coordinates",
                    "syntax": {
                        "highlighted_atoms": "|ha:0,1,2|",
                        "highlighted_bonds": "|hb:0,1|",
                        "coordinates_2d": "|c:x1,y1,x2,y2,...|",
                        "coordinates_3d": "|C:x1,y1,z1,x2,y2,z2,...|",
                    },
                    "example": "c1ccccc1 |ha:0,1,2|",
                },
                "abbreviations": {
                    "description": "Automatic chemical abbreviation system",
                    "modes": {
                        "off": "No abbreviations",
                        "groups": "Only functional groups (Ph, Me, Et, Boc, Fmoc, etc.)",
                        "reagents": "Common reagents (THF, DMF, DCM, NaOH, LiAlH4, etc.)",
                        "on": "Both groups and reagents",
                    },
                    "examples": {
                        "phenyl": "c1ccccc1 → Ph",
                        "thf": "C1CCOC1 → THF",
                        "sodium_hydroxide": "[Na+].[OH-] → NaOH",
                    },
                    "total_abbreviations": "198 reagents + 60 groups",
                },
                "dative_bonds": {
                    "description": "Coordinate bond perception with arrow notation",
                    "modes": {
                        "never": "No dative bond perception",
                        "metals": "Only metal-ligand bonds (default)",
                        "always": "All dative bonds including B, O",
                    },
                    "examples": [
                        "[NH3]BF3 (ammonia borane)",
                        "[Co][N+]([O-])(=O) (nitro complex)",
                        "O=[N+]([O-])O (nitrate)",
                    ],
                    "arrow_notation": "Donor→Acceptor",
                },
                "multicenter_bonds": {
                    "description": "η-complex and π-bonding display",
                    "styles": {
                        "provided": "As-is from input",
                        "dative": "Arrow notation from ring to metal",
                        "dashed": "Dashed lines (preserve charges)",
                        "dashed_neutral": "Dashed lines (neutralize charges)",
                        "hidden": "Hide bonds (preserve charges)",
                        "hidden_neutral": "Hide bonds (neutralize charges)",
                    },
                    "examples": [
                        "η5-Cyclopentadienyl (ferrocene)",
                        "η6-Benzene (chromium complex)",
                        "η3-Allyl (palladium complex)",
                    ],
                },
            },
            # ========== PHASE 2 ==========
            "phase2": {
                "annotations": {
                    "description": "Comprehensive annotation system for molecular structures",
                    "modes": {
                        "none": "No annotations",
                        "number": "Show atom numbers (0-based indexing)",
                        "bondnumber": "Show bond numbers",
                        "mapidx": "Show atom mapping numbers (for reactions)",
                        "atomvalue": "Show atom values/properties",
                        "colmap": "Color-code atoms by mapping number",
                        "cip": "Show CIP stereochemistry labels (R/S, E/Z)",
                    },
                    "use_cases": [
                        "Education: Show atom numbers for teaching",
                        "Reactions: Use mapidx for reaction mapping",
                        "Stereochemistry: Use cip for R/S, E/Z labels",
                    ],
                },
                "style_presets": {
                    "description": "Predefined color schemes for different contexts",
                    "presets": {
                        "cow": "Color on White (default) - Standard publication style",
                        "cob": "Color on Black - For dark backgrounds",
                        "cot": "Color on Transparent - For overlays",
                        "bow": "Black on White - High contrast, printer-friendly",
                        "bot": "Black on Transparent - Simple overlays",
                        "wob": "White on Black - Inverted high contrast",
                        "nob": "Neon on Black - Eye-catching presentations",
                    },
                    "use_cases": [
                        "Publications: cow or bow",
                        "Presentations: cob or nob",
                        "Web overlays: cot or bot",
                    ],
                },
                "aromatic_display": {
                    "description": "Circle-in-ring (donut) representation for aromatic systems",
                    "modes": {
                        "false": "Alternating double bonds (Kekulé structure)",
                        "true": "Circle inside aromatic rings (donut representation)",
                    },
                    "examples": [
                        "Benzene: 6-membered ring with inner circle",
                        "Naphthalene: Fused rings with circles",
                        "Pyridine: Nitrogen-containing aromatic",
                    ],
                },
                "reaction_arrows": {
                    "description": "Different arrow types for chemical reactions",
                    "types": {
                        "forward": "Standard forward arrow (→)",
                        "equ": "Equilibrium/reversible arrow (⇌)",
                        "ngo": "No-go/blocked arrow (⇏)",
                        "ret": "Retrosynthetic arrow (⇒)",
                        "res": "Resonance arrow (↔)",
                    },
                    "use_cases": [
                        "Standard reactions: forward",
                        "Equilibria: equ",
                        "Retrosynthesis: ret",
                        "Resonance structures: res",
                    ],
                },
                "custom_colors": {
                    "description": "Override default colors with custom hex values",
                    "parameters": {
                        "bgcolor": "Background color (e.g., '#FFFFFF')",
                        "fgcolor": "Foreground/annotation color (e.g., '#000000')",
                    },
                },
            },
            # ========== PHASE 3 ==========
            "phase3": {
                "radical_perception": {
                    "description": "Automatic detection and marking of unpaired electrons",
                    "examples": [
                        "Methyl radical: [CH3]",
                        "Oxygen radical: [O]",
                        "Nitric oxide: [N]=O",
                    ],
                    "display": "Dots or markers on atoms with unpaired electrons",
                },
                "mdl_hilite": {
                    "description": "Support for MDL/SDF V3000 HILITE properties",
                    "format": "M  V30 HIGHLIGHT ATOMS=(count index1 index2 ...)",
                    "source": "Automatically parsed from molecule properties",
                    "use_case": "Import highlighting from chemical databases",
                },
                "advanced_controls": {
                    "zoom": {
                        "description": "Zoom level for depiction",
                        "range": "0.1 to 5.0",
                        "default": 1.3,
                        "use_case": "Adjust molecule size in image",
                    },
                    "ratio": {
                        "description": "Bond thickness/stroke ratio",
                        "range": "0.5 to 2.0",
                        "default": 1.0,
                        "use_case": "Thicker bonds for presentations, thinner for publications",
                    },
                    "flip": {
                        "description": "Horizontally flip molecular structure",
                        "use_case": "Mirror image of molecule",
                    },
                    "anon": {
                        "description": "Anonymous atom display (IUPAC recommendations)",
                        "use_case": "Hide atom labels for simplified view",
                    },
                    "smalim": {
                        "description": "SMARTS hit limit",
                        "range": "1 to 1000",
                        "default": 100,
                        "use_case": "Control maximum number of highlighted matches",
                    },
                    "svgunits": {
                        "description": "SVG coordinate units",
                        "options": ["px", "mm", "cm", "in"],
                        "default": "px",
                        "use_case": "Export for print (mm/cm/in) or web (px)",
                    },
                },
            },
        },
        "usage_tips": [
            "Use CXSMILES for precise highlighting: 'CC(O)C |ha:1|'",
            "Combine abbreviations with highlighting for clean depictions",
            "CDK toolkit required for all Phase 1-3 features",
            "Use 'dative=metals' as default for coordination chemistry",
            "Apply 'multicenter=dative' for clearer organometallic structures",
            "Use 'style=cow' for publications, 'style=cob' for presentations",
            "Enable 'donuts=true' for cleaner aromatic ring display",
            "Use 'annotate=cip' to show stereochemistry labels",
            "Adjust 'zoom' and 'ratio' for optimal depiction size and clarity",
            "Use 'svgunits=mm' for print-ready output",
        ],
        "endpoint": "/depict/2D_enhanced",
        "recommended_settings": {
            "publication": {
                "style": "cow",
                "donuts": False,
                "zoom": 1.3,
                "ratio": 0.8,
                "annotate": "none",
            },
            "presentation": {
                "style": "cob",
                "donuts": True,
                "zoom": 1.5,
                "ratio": 1.2,
                "annotate": "cip",
            },
            "teaching": {
                "style": "cow",
                "donuts": False,
                "zoom": 1.5,
                "ratio": 1.0,
                "annotate": "number",
            },
            "coordination_chemistry": {
                "style": "cow",
                "dative": "metals",
                "multicenter": "dative",
                "abbreviate": "off",
                "zoom": 1.3,
            },
        },
    }
