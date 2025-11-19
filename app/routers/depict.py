from __future__ import annotations

import logging
from typing import Literal
from typing import Optional

from fastapi import APIRouter
from fastapi import HTTPException
from fastapi import Query
from fastapi import Request
from fastapi import status
from fastapi.responses import Response
from fastapi.templating import Jinja2Templates

from app.modules.depiction import get_rdkit_depiction
from app.modules.depiction_enhanced import get_cdk_depiction
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.openbabel_wrapper import get_ob_mol
from app.modules.toolkits.rdkit_wrapper import get_3d_conformers
from app.schemas import HealthCheck
from app.schemas.depict_schema import Depict2DResponse, Depict3DResponse
from app.schemas.error import BadRequestModel, ErrorResponse, NotFoundModel

# Use the shared limiter instance
from app.limiter import limiter

# Configure logging
logger = logging.getLogger(__name__)

templates = Jinja2Templates(directory="app/templates")


# Removed local FastAPI app instance and limiter/exception handler setup in favour of shared one.

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
    summary="Advanced 2D molecular depiction with comprehensive rendering options",
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
        True,
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
    # ========== CHEMICAL STRUCTURE ENHANCEMENTS ==========
    abbreviate: Literal["off", "groups", "reagents", "on"] = Query(
        default="off",
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
    # ========== VISUAL STYLE AND ANNOTATIONS ==========
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
    title: Optional[str] = Query(
        None,
        title="Title",
        description=(
            "Optional title to display when showtitle=true. "
            "If not provided, title will be extracted from SMILES string if present. "
            "Example: ?smiles=CCO&title=Ethanol&showtitle=true"
        ),
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
    # ========== ADVANCED RENDERING CONTROLS ==========
    zoom: float = Query(
        1.0,
        title="Zoom Level",
        description="Zoom level for depiction (0.1 to 5.0). Default: 1.0.",
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
    """Generate advanced 2D molecular depictions with comprehensive customization options.

    This endpoint provides professional-grade molecular structure rendering with extensive
    control over visualization, suitable for publications, presentations, and chemical databases.

    ## Key Features

    ### Extended SMILES Support (CXSMILES)
    Supports CXSMILES extensions for precise control over highlighting and atom coordinates:
    - Atom highlighting: `c1ccccc1 |ha:0,1,2|` highlights atoms 0, 1, and 2
    - Bond highlighting: `CCO |hb:0,1|` highlights bonds 0 and 1
    - 2D coordinates: `|c:x1,y1,x2,y2,...|`
    - 3D coordinates: `|C:x1,y1,z1,x2,y2,z2,...|`

    ### Chemical Abbreviations
    Automatic recognition and rendering of common chemical groups and reagents:
    - **groups**: Functional groups like Ph (phenyl), Me (methyl), Et (ethyl), Boc, Fmoc, etc.
    - **reagents**: Common reagents like THF, DMF, DCM, NaOH, LiAlH₄, etc.
    - Total library: 198 reagents + 60 functional groups
    - Example: `C1CCOC1` → "THF" when `abbreviate=reagents`

    ### Coordination Chemistry
    Advanced support for metal complexes and coordination compounds:
    - **Dative bonds**: Coordinate bonds shown with arrow notation (Donor→Acceptor)
      * Example: `[NH3]BF3` shows ammonia-borane with N→B arrow
    - **Multicenter bonds**: η-complexes (ferrocene, chromium-benzene) with multiple display styles
      * Dative arrows, dashed lines, or hidden bonds
      * Charge neutralization options

    ### Visual Styles and Annotations
    Seven predefined color schemes for different contexts:
    - **cow**: Color on White (default) - publications
    - **cob**: Color on Black - presentations with dark backgrounds
    - **bow**: Black on White - high contrast, printer-friendly
    - **nob**: Neon on Black - eye-catching presentations
    - Custom colors: Override with hex values via `bgcolor` and `fgcolor`

    Seven annotation modes:
    - **number**: Atom numbering (0-based indexing) for education
    - **cip**: CIP stereochemistry labels (R/S, E/Z)
    - **mapidx**: Reaction atom mapping
    - **colmap**: Color-coded atom mapping
    - **bondnumber**: Bond numbering
    - **atomvalue**: Atom properties display

    ### Aromatic Ring Display
    Choose between Kekulé structures (alternating double bonds) or circle-in-ring
    representation using the `donuts` parameter for cleaner aromatic visualization.

    ### Reaction Depiction
    Multiple arrow types for chemical reactions:
    - **forward** (→): Standard reaction
    - **equ** (⇌): Equilibrium/reversible
    - **ret** (⇒): Retrosynthetic analysis
    - **res** (↔): Resonance structures
    - **ngo** (⇏): No-go/blocked reaction

    ### Professional Controls
    - **Radical perception**: Automatic detection of unpaired electrons
    - **MDL HILITE**: Support for MDL/SDF V3000 highlighting from databases
    - **Zoom/Ratio**: Fine control over molecule size and bond thickness
    - **Flip**: Horizontal structure mirroring
    - **SVG units**: Export in pixels, millimeters, centimeters, or inches

    ## Usage Examples

    ```
    # Basic molecule
    ?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C

    # CXSMILES with highlighting
    ?smiles=c1ccccc1 |ha:0,1,2|

    # Abbreviated reagent
    ?smiles=C1CCOC1&abbreviate=reagents  # Shows "THF"

    # Coordination complex with dative bonds
    ?smiles=[Co][N+]([O-])(=O)&dative=metals

    # Publication-ready (black on white, no abbreviations)
    ?smiles=c1ccccc1&style=bow&abbreviate=off&zoom=1.3&ratio=0.8

    # Presentation style (color on black, larger, annotated)
    ?smiles=C[C@H](N)C(=O)O&style=cob&annotate=cip&zoom=1.5&ratio=1.2

    # Teaching (atom numbers, clear display)
    ?smiles=c1ccccc1&annotate=number&zoom=1.5&showAtomNumbers=true
    ```

    ## Toolkit

    This endpoint uses the **Chemistry Development Kit (CDK)** exclusively for rendering.
    CDK provides full support for all advanced features including abbreviations, dative bonds,
    multicenter bonds, CXSMILES extensions, annotations, style presets, and advanced rendering controls.

    For basic 2D depiction with RDKit support, use the `/depict/2D` endpoint instead.

    ## Returns
    SVG image of the molecular structure with all requested enhancements and customizations.
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

        # Generate depiction using CDK (only toolkit supported for enhanced features)
        # CDK's SmilesParser automatically handles CXSMILES (e.g., "CCO |ha:0,1|")
        # The highlighting is extracted internally via extract_cxsmiles_highlighting()
        mol = parse_input(smiles, "cdk", False)
        if title:
            try:
                from jpype import JClass

                cdk_base = "org.openscience.cdk"
                CDKConstants = JClass(cdk_base + ".CDKConstants")
                IReaction = JClass(cdk_base + ".interfaces.IReaction")
                IReactionSet = JClass(cdk_base + ".interfaces.IReactionSet")

                # Handle different molecule types
                if isinstance(mol, IReactionSet):
                    # Set title on all reactions in set
                    for rxn in mol.reactions():
                        rxn.setProperty(CDKConstants.TITLE, title)
                elif isinstance(mol, IReaction):
                    # Set title on single reaction
                    mol.setProperty(CDKConstants.TITLE, title)
                else:
                    # Set title on molecule
                    mol.setProperty(CDKConstants.TITLE, title)
            except Exception:
                # Silently ignore title setting errors
                pass
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
            # Chemical structure enhancements
            abbreviate=abbreviate,
            dative=dative,
            multicenter=multicenter,
            # Visual style and annotations
            annotate=annotate,
            style=style,
            donuts=donuts,
            arrow=arrow,
            alignrxnmap=alignrxnmap,
            showtitle=showtitle,
            bgcolor=bgcolor,
            fgcolor=fgcolor,
            # Advanced rendering controls
            zoom=zoom,
            ratio=ratio,
            flip=flip,
            anon=anon,
            smalim=smalim,
            svgunits=svgunits,
            perceive_radicals=perceive_radicals,
            apply_mdl_highlighting=apply_mdl_highlighting,
        )

        return Response(content=depiction, media_type="image/svg+xml")

    except Exception as e:
        logger.error(f"Error generating enhanced depiction: {e}")
        raise HTTPException(status_code=422, detail=str(e))
