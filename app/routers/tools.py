from fastapi import APIRouter, status
from app.modules.tools.surge import generateStructures
from app.modules.tools.sugarremoval import (
    getSugarInfo,
    removeLinearSugar,
    removeCircularSugar,
    removeLinearandCircularSugar,
)
from app.schemas import HealthCheck

router = APIRouter(
    prefix="/tools",
    tags=["tools"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
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
    """
    ## Perform a Health Check
    Endpoint to perform a healthcheck on. This endpoint can primarily be used Docker
    to ensure a robust container orchestration and management is in place. Other
    services which rely on proper functioning of the API service will not deploy if this
    endpoint returns any other HTTP status code except 200 (OK).
    Returns:
        HealthCheck: Returns a JSON response with the health status
    """
    return HealthCheck(status="OK")


@router.get("/generate-structures")
async def Generate_Structures(molecular_formula: str):
    """
    Generates structures using the chemical structure generator based on the canonical generation path method.
    For more information refer to:
    - McKay, B.D., Yirik, M.A. & Steinbeck, C. Surge: a fast open-source chemical graph generator. J Cheminform 14, 24 (2022). https://doi.org/10.1186/s13321-022-00604-9


    Parameters:
     - **Molecular_Formula**: required (str): The molecular formula of the compound.

    Returns:
    - List[Structure]: A list of generated structures.

    Raises:
    - HTTPException: If there was an error generating the structures.

    Example:
    - GET /generate-structures?molecular_formula=C4H10

    Note:
    - The maximum allowable count of heavy atoms is restricted to 10 to mitigate excessive utilization of this service.

    """

    return generateStructures(molecular_formula)


@router.get("/sugar-information")
async def getsugarinformation(smiles: str):
    """
    Get information whether a given molecule has circular or linear sugars.
    For more information refer to:
    - Schaub, J., Zielesny, A., Steinbeck, C. et al. Too sweet: cheminformatics for deglycosylation in natural products. J Cheminform 12, 67 (2020). https://doi.org/10.1186/s13321-020-00467-y.

    Parameters:
    - **SMILES string**: (str): SMILES string of the molecule (required, query parameter)

    Returns:
    - str: A message indicating the type of sugars present in the molecule.

    Note:
    The function `getSugarInfo` is used internally to determine the presence of linear and circular sugars in the molecule.

    The returned message indicates the types of sugars present in the molecule:
        - If both linear and circular sugars are present, it returns "The molecule contains Linear and Circular sugars."
        - If only linear sugar is present, it returns "The molecule contains only Linear sugar."
        - If only circular sugars are present, it returns "The molecule contains only Circular sugar."
        - If no sugars are found, it returns "The molecule contains no sugar."

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
    Detect and remove linear sugars from a given SMILES string using Sugar Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES string of the molecule (required, query parameter)

    Returns:
    - str: The modified SMILES string with linear sugars removed.

    """
    return removeLinearSugar(smiles)


@router.get("/remove-circular-sugar")
async def removecircularsugars(smiles: str):
    """
    Detect and remove circular sugars from a given SMILES string using Sugar Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES string of the molecule (required, query parameter)

    Returns:
    - str: The modified SMILES string with circular sugars removed.

    """
    return removeCircularSugar(smiles)


@router.get("/remove-linearandcircular-sugar")
async def removelinearandcircularsugars(smiles: str):
    """
    Detect and remove linear and circular sugars from a given SMILES string using Sugar Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES string of the molecule (required, query parameter)

    Returns:
    - str: The modified SMILES string with linear and circular sugars removed.

    """
    return removeLinearandCircularSugar(smiles)
