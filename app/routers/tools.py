from __future__ import annotations

from fastapi import APIRouter
from fastapi import HTTPException
from fastapi import Query
from fastapi import status

from app.modules.toolkits.helpers import parse_input
from app.modules.tools.sugar_removal import extract_aglycone_and_sugars
from app.modules.tools.sugar_removal import get_sugar_info
from app.modules.tools.sugar_removal import preservation_modes_enum
from app.modules.tools.sugar_removal import remove_circular_sugars
from app.modules.tools.sugar_removal import remove_linear_and_circular_sugars
from app.modules.tools.sugar_removal import remove_linear_sugars
from app.modules.tools.surge import generate_structures_SURGE
from app.schemas import HealthCheck
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse
from app.schemas.error import NotFoundModel
from app.schemas.tools_schema import ExtractAglyconeAndSugarsResponse
from app.schemas.tools_schema import GenerateStructuresResponse
from app.schemas.tools_schema import GetCircularandLinearSugarResponse
from app.schemas.tools_schema import GetCircularSugarResponse
from app.schemas.tools_schema import GetLinearSugarResponse
from app.schemas.tools_schema import GetSugarInformationResponse

router = APIRouter(
    prefix="/tools",
    tags=["tools"],
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
    summary="Perform a Health Check on Tools Module",
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
    "/generate-structures",
    summary="Generates structures using the chemical structure generator",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateStructuresResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def generate_structures(
    molecular_formula: str = Query(
        title="Molecular Formula",
        description="Molecular Formula for the chemical structure to be generated",
        openapi_examples={
            "example1": {
                "summary": "Example: Heavy atom count 6",
                "value": "C6H6",
            },
            "example2": {"summary": "Example: Heavy atom count 8", "value": "C8H8"},
        },
    ),
):
    """Generates structures using the chemical structure generator based on the.

    canonical generation path method.

    For more information refer to:
    - McKay, B.D., Yirik, M.A. & Steinbeck, C. Surge: a fast open-source chemical graph generator. J Cheminform 14, 24 (2022). https://doi.org/10.1186/s13321-022-00604-9

    Parameters:
     - **Molecular_Formula**: required (str): The molecular formula of the compound.

    Returns:
    - Dict: A dictionary containing structure generation results with:
        - total_count: Total number of possible structures
        - generated_count: Number of structures actually generated
        - structures: List of SMILES strings (limited to 1000)
        - settings: Dictionary describing the surge settings used
        - formula: The input molecular formula
        - limit_applied: Whether a limit was applied to results

    Raises:
    - HTTPException: If there was an error generating the structures.

    Example:
    - GET /generate-structures?molecular_formula=C4H10

    Note:
    - The maximum allowable count of heavy atoms is restricted to 10 to mitigate excessive utilization of this service.
    - Results are limited to the first 1000 structures when the total count exceeds this limit.
    """
    try:
        result = generate_structures_SURGE(molecular_formula)
        if isinstance(result, str):
            # Error message returned
            raise HTTPException(status_code=400, detail=result)
        else:
            # Success - return the structured response
            return GenerateStructuresResponse(message="Success", output=result)
    except HTTPException:
        raise  # Re-raise HTTP exceptions
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/sugars-info",
    response_model=str,
    summary="Get information whether a given molecule has circular or linear sugars",
    responses={
        200: {
            "description": "Successful response",
            "model": GetSugarInformationResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def get_sugar_info_endpoint(
    smiles: str = Query(
        title="SMILES",
        description="SMILES: string representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: 1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]pentane-1,2,3,4,5-pentol",
                "value": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
            },
            "example2": {
                "summary": "Example: 5-[1-(3,4-dihydroxyphenyl)-3-(2,3,4,5,6,7-hexahydroxyheptoxycarbonyl)-6-hydroxy-7-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,2-dihydronaphthalene-2-carbonyl]oxy-3,4-dihydroxycyclohexene-1-carboxylic acid",
                "value": "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
            },
        },
    ),
    gly_bond: bool = Query(
        default=False,
        title="Detect only Circular Sugars with O-Glycosidic Bonds",
        description="Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.",
    ),
    oxygen_atoms: bool = Query(
        default=True,
        title="Detect only Circular Sugars with enough exocyclic Oxygen Atoms",
        description="Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.",
    ),
    oxygen_atoms_threshold: float = Query(
        default=0.5,
        minimum=0.0,
        maximum=1.0,
        title="Exocyclic Oxygen Atoms to Atoms in Ring Ratio Threshold",
        description="A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!",
    ),
    linear_sugars_in_rings: bool = Query(
        default=False,
        title="Detect Linear Sugars in Rings",
        description="Whether to consider linear sugars in rings. Default is False.",
    ),
    linear_sugars_min_size: int = Query(
        default=4,
        minimum=0,
        title="Linear Sugars Minimum Size",
        description="Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.",
    ),
    linear_sugars_max_size: int = Query(
        default=7,
        minimum=1,
        title="Linear Sugars Maximum Size",
        description="Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.",
    ),
    linear_acidic_sugars: bool = Query(
        default=False,
        title="Detect Linear Acidic Sugars",
        description="Whether to consider linear acidic sugars. Default is False.",
    ),
    spiro_sugars: bool = Query(
        default=False,
        title="Detect Spiro Sugars",
        description="Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.",
    ),
    keto_sugars: bool = Query(
        default=False,
        title="Detect Keto Sugars",
        description="Whether circular sugars with keto groups should be detected. Default is False.",
    ),
):
    """
    Get information on whether a given molecule has circular or linear sugars.

    For more information refer to:
    - Schaub, J., Zielesny, A., Steinbeck, C., Sorokina, M. Too sweet: cheminformatics for deglycosylation in natural products. J Cheminform 12, 67 (2020). https://doi.org/10.1186/s13321-020-00467-y.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)
    - **gly_bond**: (bool): Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.
    - **oxygen_atoms**: (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.
    - **oxygen_atoms_threshold**: (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
    - **linear_sugars_in_rings**: (bool): Whether to consider linear sugars in rings. Default is False.
    - **linear_sugars_min_size**: (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.
    - **linear_sugars_max_size**: (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.
    - **linear_acidic_sugars**: (bool): Whether to consider linear acidic sugars. Default is False.
    - **spiro_sugars**: (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.
    - **keto_sugars**: (bool): Whether circular sugars with keto groups should be detected. Default is False.

    Returns:
    - str: A message indicating the type of sugars present in the molecule.

    Note:
    The function `get_sugar_info` is used internally to determine the presence of linear and circular sugars in the molecule.

    The returned message indicates the types of sugars present in the molecule:
        - If both linear and circular sugars are present, it returns "The molecule contains Linear and Circular sugars."
        - If only linear sugar is present, it returns "The molecule contains only Linear sugar."
        - If only circular sugars are present, it returns "The molecule contains only Circular sugar."
        - If no sugars are found, it returns "The molecule contains no sugar."
    """
    try:
        mol = parse_input(smiles, "cdk", False)
        _has_linear_sugar, _has_circular_sugars = get_sugar_info(
            molecule=mol,
            gly_bond=gly_bond,
            oxygen_atoms=oxygen_atoms,
            oxygen_atoms_threshold=oxygen_atoms_threshold,
            linear_sugars_in_rings=linear_sugars_in_rings,
            linear_sugars_min_size=linear_sugars_min_size,
            linear_sugars_max_size=linear_sugars_max_size,
            linear_acidic_sugars=linear_acidic_sugars,
            spiro_sugars=spiro_sugars,
            keto_sugars=keto_sugars,
        )
        if _has_linear_sugar and _has_circular_sugars:
            return "The molecule contains Linear and Circular sugars"
        if _has_linear_sugar and not _has_circular_sugars:
            return "The molecule contains only Linear sugar"
        if _has_circular_sugars and not _has_linear_sugar:
            return "The molecule contains only Circular sugar"
        else:
            return "The molecule contains no sugar"
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/remove-linear-sugars",
    summary="Detect and remove linear sugars",
    responses={
        200: {"description": "Successful response", "model": GetLinearSugarResponse},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def remove_linear_sugars_endpoint(
    smiles: str = Query(
        title="SMILES",
        description="SMILES: string representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: 1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]pentane-1,2,3,4,5-pentol",
                "value": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
            },
            "example2": {
                "summary": "Example: 5-[1-(3,4-dihydroxyphenyl)-3-(2,3,4,5,6,7-hexahydroxyheptoxycarbonyl)-6-hydroxy-7-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,2-dihydronaphthalene-2-carbonyl]oxy-3,4-dihydroxycyclohexene-1-carboxylic acid",
                "value": "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
            },
        },
    ),
    only_terminal: bool = Query(
        default=True,
        title="Remove Only Terminal Sugars",
        description="Whether only terminal sugars should be removed. Default is True.",
    ),
    preservation_mode: int = Query(
        default=2,
        minimum=1,
        maximum=3,
        title="Preservation Mode",
        description="Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).",
    ),
    preservation_threshold: int = Query(
        default=5,
        minimum=0,
        title="Preservation Mode Threshold",
        description="Threshold value for the selected preservation mode. Default is 5 (heavy atoms).",
    ),
    linear_sugars_in_rings: bool = Query(
        default=False,
        title="Detect Linear Sugars in Rings",
        description="Whether to consider linear sugars in rings. Default is False.",
    ),
    linear_sugars_min_size: int = Query(
        default=4,
        minimum=0,
        title="Linear Sugars Minimum Size",
        description="Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.",
    ),
    linear_sugars_max_size: int = Query(
        default=7,
        minimum=1,
        title="Linear Sugars Maximum Size",
        description="Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.",
    ),
    linear_acidic_sugars: bool = Query(
        default=False,
        title="Detect Linear Acidic Sugars",
        description="Whether to consider linear acidic sugars. Default is False.",
    ),
    mark_attach_points: bool = Query(
        default=False,
        title="Mark Attachment Points",
        description="Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.",
    ),
):
    """
    Detect and remove linear sugars from a given SMILES string using Sugar Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)
    - **only_terminal**: (bool): Whether only terminal linear sugars should be removed. Default is True.
    - **preservation_mode**: (int): Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).
    - **preservation_threshold**: (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
    - **linear_sugars_in_rings**: (bool): Whether to consider linear sugars in rings. Default is False.
    - **linear_sugars_min_size**: (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.
    - **linear_sugars_max_size**: (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.
    - **linear_acidic_sugars**: (bool): Whether to consider linear acidic sugars. Default is False.
    - **mark_attach_points**: (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.

    Returns:
    - str: The modified SMILES string with linear sugars removed.
    """
    try:
        mol = parse_input(smiles, "cdk", False)
        _removed_smiles = remove_linear_sugars(
            molecule=mol,
            only_terminal=only_terminal,
            preservation_mode=preservation_modes_enum(preservation_mode),
            preservation_threshold=preservation_threshold,
            linear_sugars_in_rings=linear_sugars_in_rings,
            linear_sugars_min_size=linear_sugars_min_size,
            linear_sugars_max_size=linear_sugars_max_size,
            linear_acidic_sugars=linear_acidic_sugars,
            mark_attach_points=mark_attach_points,
        )
        if _removed_smiles:
            return _removed_smiles
        else:
            raise HTTPException(
                status_code=422,
                detail="Error reading SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/remove-circular-sugars",
    summary="Detect and remove linear sugars",
    responses={
        200: {"description": "Successful response", "model": GetCircularSugarResponse},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def remove_circular_sugars_endpoint(
    smiles: str = Query(
        title="SMILES",
        description="SMILES: string representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: 1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]pentane-1,2,3,4,5-pentol",
                "value": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
            },
            "example2": {
                "summary": "Example: 5-[1-(3,4-dihydroxyphenyl)-3-(2,3,4,5,6,7-hexahydroxyheptoxycarbonyl)-6-hydroxy-7-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,2-dihydronaphthalene-2-carbonyl]oxy-3,4-dihydroxycyclohexene-1-carboxylic acid",
                "value": "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
            },
        },
    ),
    gly_bond: bool = Query(
        default=False,
        title="Detect only Circular Sugars with O-Glycosidic Bonds",
        description="Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.",
    ),
    only_terminal: bool = Query(
        default=True,
        title="Remove Only Terminal Sugars",
        description="Whether only terminal sugars should be removed. Default is True.",
    ),
    preservation_mode: int = Query(
        default=2,
        minimum=1,
        maximum=3,
        title="Preservation Mode",
        description="Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).",
    ),
    preservation_threshold: int = Query(
        default=5,
        minimum=0,
        title="Preservation Mode Threshold",
        description="Threshold value for the selected preservation mode. Default is 5 (heavy atoms).",
    ),
    oxygen_atoms: bool = Query(
        default=True,
        title="Detect only Circular Sugars with enough exocyclic Oxygen Atoms",
        description="Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.",
    ),
    oxygen_atoms_threshold: float = Query(
        default=0.5,
        minimum=0.0,
        maximum=1.0,
        title="Exocyclic Oxygen Atoms to Atoms in Ring Ratio Threshold",
        description="A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!",
    ),
    spiro_sugars: bool = Query(
        default=False,
        title="Detect Spiro Sugars",
        description="Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.",
    ),
    keto_sugars: bool = Query(
        default=False,
        title="Detect Keto Sugars",
        description="Whether circular sugars with keto groups should be detected. Default is False.",
    ),
    mark_attach_points: bool = Query(
        default=False,
        title="Mark Attachment Points",
        description="Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.",
    ),
):
    """
    Detect and remove circular sugars from a given SMILES string using Sugar Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)
    - **gly_bond**: (bool): Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.
    - **only_terminal**: (bool): Whether only terminal sugars should be removed. Default is True.
    - **preservation_mode**: (int): Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).
    - **preservation_threshold**: (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
    - **oxygen_atoms**: (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.
    - **oxygen_atoms_threshold**: (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
    - **spiro_sugars**: (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.
    - **keto_sugars**: (bool): Whether circular sugars with keto groups should be detected. Default is False.
    - **mark_attach_points**: (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.

    Returns:
    - str: The modified SMILES string with circular sugars removed.
    """
    try:
        mol = parse_input(smiles, "cdk", False)
        removed_smiles = remove_circular_sugars(
            molecule=mol,
            gly_bond=gly_bond,
            only_terminal=only_terminal,
            preservation_mode=preservation_modes_enum(preservation_mode),
            preservation_threshold=preservation_threshold,
            oxygen_atoms=oxygen_atoms,
            oxygen_atoms_threshold=oxygen_atoms_threshold,
            spiro_sugars=spiro_sugars,
            keto_sugars=keto_sugars,
            mark_attach_points=mark_attach_points,
        )
        if removed_smiles:
            return removed_smiles
        else:
            raise HTTPException(
                status_code=422,
                detail="Error processing SMILES string.",
            )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/remove-sugars",
    summary="Detect and remove linear and circular sugars",
    responses={
        200: {
            "description": "Successful response",
            "model": GetCircularandLinearSugarResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def remove_linear_and_circular_sugars_endpoint(
    smiles: str = Query(
        title="SMILES",
        description="SMILES: string representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: 1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]pentane-1,2,3,4,5-pentol",
                "value": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
            },
            "example2": {
                "summary": "Example: 5-[1-(3,4-dihydroxyphenyl)-3-(2,3,4,5,6,7-hexahydroxyheptoxycarbonyl)-6-hydroxy-7-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,2-dihydronaphthalene-2-carbonyl]oxy-3,4-dihydroxycyclohexene-1-carboxylic acid",
                "value": "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
            },
        },
    ),
    gly_bond: bool = Query(
        default=False,
        title="Detect only Circular Sugars with O-Glycosidic Bonds",
        description="Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.",
    ),
    only_terminal: bool = Query(
        default=True,
        title="Remove Only Terminal Sugars",
        description="Whether only terminal sugars should be removed. Default is True.",
    ),
    preservation_mode: int = Query(
        default=2,
        minimum=1,
        maximum=3,
        title="Preservation Mode",
        description="Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).",
    ),
    preservation_threshold: int = Query(
        default=5,
        minimum=0,
        title="Preservation Mode Threshold",
        description="Threshold value for the selected preservation mode. Default is 5 (heavy atoms).",
    ),
    oxygen_atoms: bool = Query(
        default=True,
        title="Detect only Circular Sugars with enough exocyclic Oxygen Atoms",
        description="Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.",
    ),
    oxygen_atoms_threshold: float = Query(
        default=0.5,
        minimum=0.0,
        maximum=1.0,
        title="Exocyclic Oxygen Atoms to Atoms in Ring Ratio Threshold",
        description="A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!",
    ),
    linear_sugars_in_rings: bool = Query(
        default=False,
        title="Detect Linear Sugars in Rings",
        description="Whether to consider linear sugars in rings. Default is False.",
    ),
    linear_sugars_min_size: int = Query(
        default=4,
        minimum=0,
        title="Linear Sugars Minimum Size",
        description="Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.",
    ),
    linear_sugars_max_size: int = Query(
        default=7,
        minimum=1,
        title="Linear Sugars Maximum Size",
        description="Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.",
    ),
    linear_acidic_sugars: bool = Query(
        default=False,
        title="Detect Linear Acidic Sugars",
        description="Whether to consider linear acidic sugars. Default is False.",
    ),
    spiro_sugars: bool = Query(
        default=False,
        title="Detect Spiro Sugars",
        description="Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.",
    ),
    keto_sugars: bool = Query(
        default=False,
        title="Detect Keto Sugars",
        description="Whether circular sugars with keto groups should be detected. Default is False.",
    ),
    mark_attach_points: bool = Query(
        default=False,
        title="Mark Attachment Points",
        description="Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.",
    ),
):
    """
    Detect and remove linear and circular sugars from a given SMILES string using Sugar Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)
    - **gly_bond**: (bool): Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.
    - **only_terminal**: (bool): Whether only terminal sugars should be removed. Default is True.
    - **preservation_mode**: (int): Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).
    - **preservation_threshold**: (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
    - **oxygen_atoms**: (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.
    - **oxygen_atoms_threshold**: (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
    - **linear_sugars_in_rings**: (bool): Whether to consider linear sugars in rings. Default is False.
    - **linear_sugars_min_size**: (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.
    - **linear_sugars_max_size**: (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.
    - **linear_acidic_sugars**: (bool): Whether to consider linear acidic sugars. Default is False.
    - **spiro_sugars**: (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.
    - **keto_sugars**: (bool): Whether circular sugars with keto groups should be detected. Default is False.
    - **mark_attach_points**: (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.

    Returns:
    - str: The modified SMILES string with linear and circular sugars removed.
    """
    try:
        mol = parse_input(smiles, "cdk", False)
        _removed_smiles = remove_linear_and_circular_sugars(
            molecule=mol,
            gly_bond=gly_bond,
            only_terminal=only_terminal,
            preservation_mode=preservation_modes_enum(preservation_mode),
            preservation_threshold=preservation_threshold,
            oxygen_atoms=oxygen_atoms,
            oxygen_atoms_threshold=oxygen_atoms_threshold,
            linear_sugars_in_rings=linear_sugars_in_rings,
            linear_sugars_min_size=linear_sugars_min_size,
            linear_sugars_max_size=linear_sugars_max_size,
            linear_acidic_sugars=linear_acidic_sugars,
            spiro_sugars=spiro_sugars,
            keto_sugars=keto_sugars,
            mark_attach_points=mark_attach_points,
        )
        if _removed_smiles:
            return _removed_smiles
        else:
            raise HTTPException(
                status_code=422,
                detail="Error processing SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/extract-aglycone-and-sugars",
    summary='Extracts the aglycone and the sugars from a given molecule. Returns a printed tuple (["<SMILES>", "<SMILES>", ...]). The first position is always the aglycone.',
    responses={
        200: {
            "description": "Successful response",
            "model": ExtractAglyconeAndSugarsResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def extract_aglycone_and_sugars_endpoint(
    smiles: str = Query(
        title="SMILES",
        description="SMILES: string representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: 1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]pentane-1,2,3,4,5-pentol",
                "value": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
            },
            "example2": {
                "summary": "Example: 5-[1-(3,4-dihydroxyphenyl)-3-(2,3,4,5,6,7-hexahydroxyheptoxycarbonyl)-6-hydroxy-7-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-1,2-dihydronaphthalene-2-carbonyl]oxy-3,4-dihydroxycyclohexene-1-carboxylic acid",
                "value": "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
            },
        },
    ),
    extract_circular_sugars: bool = Query(
        default=True,
        title="Extract Circular Sugars",
        description="Whether to extract circular sugars. Default is True.",
    ),
    extract_linear_sugars: bool = Query(
        default=True,
        title="Extract Linear Sugars",
        description="Whether to extract linear sugars. Default is True.",
    ),
    gly_bond: bool = Query(
        default=False,
        title="Detect only Circular Sugars with O-Glycosidic Bonds",
        description="Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.",
    ),
    only_terminal: bool = Query(
        default=True,
        title="Remove Only Terminal Sugars",
        description="Whether only terminal sugars should be removed. Default is True.",
    ),
    preservation_mode: int = Query(
        default=2,
        minimum=1,
        maximum=3,
        title="Preservation Mode",
        description="Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).",
    ),
    preservation_threshold: int = Query(
        default=5,
        minimum=0,
        title="Preservation Mode Threshold",
        description="Threshold value for the selected preservation mode. Default is 5 (heavy atoms).",
    ),
    oxygen_atoms: bool = Query(
        default=True,
        title="Detect only Circular Sugars with enough exocyclic Oxygen Atoms",
        description="Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.",
    ),
    oxygen_atoms_threshold: float = Query(
        default=0.5,
        minimum=0.0,
        maximum=1.0,
        title="Exocyclic Oxygen Atoms to Atoms in Ring Ratio Threshold",
        description="A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!",
    ),
    linear_sugars_in_rings: bool = Query(
        default=False,
        title="Detect Linear Sugars in Rings",
        description="Whether to consider linear sugars in rings. Default is False.",
    ),
    linear_sugars_min_size: int = Query(
        default=4,
        minimum=0,
        title="Linear Sugars Minimum Size",
        description="Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.",
    ),
    linear_sugars_max_size: int = Query(
        default=7,
        minimum=1,
        title="Linear Sugars Maximum Size",
        description="Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.",
    ),
    linear_acidic_sugars: bool = Query(
        default=False,
        title="Detect Linear Acidic Sugars",
        description="Whether to consider linear acidic sugars. Default is False.",
    ),
    spiro_sugars: bool = Query(
        default=False,
        title="Detect Spiro Sugars",
        description="Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.",
    ),
    keto_sugars: bool = Query(
        default=False,
        title="Detect Keto Sugars",
        description="Whether circular sugars with keto groups should be detected. Default is False.",
    ),
    mark_attach_points: bool = Query(
        default=False,
        title="Mark Attachment Points",
        description="Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.",
    ),
    post_process_sugars: bool = Query(
        default=False,
        title="Post-process Sugars",
        description="Whether the extracted sugar moieties should be post-processed, i.e. bond splitting (O-glycosidic, ether, ester, peroxide) to separate the individual sugars, before being output. Default is False.",
    ),
    limit_post_process_by_size: bool = Query(
        default=False,
        title="Limit Post-processing by Size",
        description="Whether the post-processing of extracted sugar moieties should be limited to structures bigger than a defined size (see preservation mode (threshold)) to preserve smaller modifications. Default is False.",
    ),
):
    """
    Extracts the aglycone and sugars from a given SMILES string using Sugar Detection Utility.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)
    - **extract_circular_sugars**: (bool): Whether to extract circular sugars. Default is True.
    - **extract_linear_sugars**: (bool): Whether to extract linear sugars. Default is False.
    - **gly_bond**: (bool): Whether to consider only circular sugars with glycosidic bonds in the analysis. Default is False.
    - **only_terminal**: (bool): Whether only terminal sugars should be removed. Default is True.
    - **preservation_mode**: (int): Mode to determine which disconnected structures to preserve. All (1): Preserve all disconnected structures (note: this might lead to no circular sugar moieties being detected, depending on the other settings). Heavy atom count (2): Remove disconnected structures that do not have enough heavy atoms. Molecular weight (3): Remove disconnected structures that do not have a sufficient molecular weight. Default is heavy atom count (2).
    - **preservation_threshold**: (int): Threshold value for the selected preservation mode. Default is 5 (heavy atoms).
    - **oxygen_atoms**: (bool): Whether to consider only circular sugars with a sufficient number of exocyclic oxygen atoms in the analysis (see oxygen_atoms_threshold). Default is True.
    - **oxygen_atoms_threshold**: (float): A number giving the minimum attached exocyclic oxygen atoms to atom number in the ring ratio a circular sugar needs to have to be considered in the analysis. Default is 0.5 (a 6-membered ring needs at least 3 attached exocyclic oxygen atoms). Must be positive!
    - **linear_sugars_in_rings**: (bool): Whether to consider linear sugars in rings. Default is False.
    - **linear_sugars_min_size**: (int): Minimum size of linear sugars to consider. Default is 4. Must be positive and higher than or equal to 0 and also smaller than the linear sugar candidate maximum size.
    - **linear_sugars_max_size**: (int): Maximum size of linear sugars to consider. Default is 7. Must be positive and higher than or equal to 1 and also higher than the linear sugar candidate minimum size.
    - **linear_acidic_sugars**: (bool): Whether to consider linear acidic sugars. Default is False.
    - **spiro_sugars**: (bool): Whether spiro rings (rings that share one atom with another cycle) should be included in the circular sugar detection. Default is False.
    - **keto_sugars**: (bool): Whether circular sugars with keto groups should be detected. Default is False.
    - **mark_attach_points**: (bool): Whether to mark the attachment points of removed sugars with a dummy atom. Default is False.
    - **post_process_sugars**: (bool): Whether the extracted sugar moieties should be post-processed, i.e. bond splitting (O-glycosidic, ether, ester, peroxide) to separate the individual sugars, before being output. Default is False.
    - **limit_post_process_by_size**: (bool): Whether the post-processing of extracted sugar moieties should be limited to structures bigger than a defined size (see preservation mode (threshold)) to preserve smaller modifications. Default is False.

    Returns:
    - tuple: The SMILES representations of the aglycone and sugars. The first one is always the aglycone.
    """
    try:
        mol = parse_input(smiles, "cdk", False)
        _aglycone_and_sugars_tuple = extract_aglycone_and_sugars(
            molecule=mol,
            extract_circular_sugars=extract_circular_sugars,
            extract_linear_sugars=extract_linear_sugars,
            gly_bond=gly_bond,
            only_terminal=only_terminal,
            preservation_mode=preservation_modes_enum(preservation_mode),
            preservation_threshold=preservation_threshold,
            oxygen_atoms=oxygen_atoms,
            oxygen_atoms_threshold=oxygen_atoms_threshold,
            linear_sugars_in_rings=linear_sugars_in_rings,
            linear_sugars_min_size=linear_sugars_min_size,
            linear_sugars_max_size=linear_sugars_max_size,
            linear_acidic_sugars=linear_acidic_sugars,
            spiro_sugars=spiro_sugars,
            keto_sugars=keto_sugars,
            mark_attach_points=mark_attach_points,
            post_process_sugars=post_process_sugars,
            limit_post_process_by_size=limit_post_process_by_size,
        )
        if _aglycone_and_sugars_tuple:
            return list(_aglycone_and_sugars_tuple)
        else:
            raise HTTPException(
                status_code=422,
                detail="Error processing SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
