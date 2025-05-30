from __future__ import annotations
import os
import uuid
from PyPDF2 import PdfReader
from fastapi import (
    UploadFile,
    File,
    Form,
)
from fastapi import APIRouter
from fastapi import HTTPException
from fastapi import Query
from fastapi import status

from app.modules.toolkits.helpers import parse_input
from app.modules.tools.sugar_removal import get_sugar_info
from app.modules.tools.sugar_removal import remove_circular_sugar
from app.modules.tools.sugar_removal import remove_linear_and_circular_sugar
from app.modules.tools.sugar_removal import remove_linear_sugar
from app.modules.tools.surge import generate_structures_SURGE
from app.schemas import HealthCheck
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse
from app.schemas.error import NotFoundModel
from app.schemas.tools_schema import GenerateStructuresResponse
from app.schemas.tools_schema import GetCircularandLinearSugarResponse
from app.schemas.tools_schema import GetCircularSugarResponse
from app.schemas.tools_schema import GetLinearSugarResponse
from app.schemas.tools_schema import GetSugarInformationResponse
from app.config import PDF_DIR
from app.modules.docling_wrapper import (
    get_converted_document,
    extract_from_docling_document,
    combine_to_paragraph,
)

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
    operation_id="generate_structures",
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
    - List[str]: A list of generated structures.

    Raises:
    - HTTPException: If there was an error generating the structures.

    Example:
    - GET /generate-structures?molecular_formula=C4H10

    Note:
    - The maximum allowable count of heavy atoms is restricted to 10 to mitigate excessive utilization of this service.
    """
    try:
        structures = generate_structures_SURGE(molecular_formula)
        if structures:
            return structures
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/sugars-info",
    operation_id="get_sugar_information",
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
async def get_sugar_information(
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
):
    """Get information on whether a given molecule has circular or linear.

    sugars.

    For more information refer to:
    - Schaub, J., Zielesny, A., Steinbeck, C. et al. Too sweet: cheminformatics for deglycosylation in natural products. J Cheminform 12, 67 (2020). https://doi.org/10.1186/s13321-020-00467-y.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)

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
    mol = parse_input(smiles, "cdk", False)
    try:
        hasLinearSugar, hasCircularSugars = get_sugar_info(mol)
        if hasLinearSugar and hasCircularSugars:
            return "The molecule contains Linear and Circular sugars"
        if hasLinearSugar and not hasCircularSugars:
            return "The molecule contains only Linear sugar"
        if hasCircularSugars and not hasLinearSugar:
            return "The molecule contains only Circular sugar"
        else:
            return "The molecule contains no sugar"
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/remove-linear-sugars",
    operation_id="remove_linear_sugars",
    summary="Detect and remove linear sugars",
    responses={
        200: {"description": "Successful response", "model": GetLinearSugarResponse},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def remove_linear_sugars(
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
):
    """Detect and remove linear sugars from a given SMILES string using Sugar.

    Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)

    Returns:
    - str: The modified SMILES string with linear sugars removed.
    """
    mol = parse_input(smiles, "cdk", False)
    try:
        removed_smiles = remove_linear_sugar(mol)
        if removed_smiles:
            return removed_smiles
        else:
            raise HTTPException(
                status_code=422,
                detail="Error reading SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/remove-circular-sugars",
    operation_id="remove_circular_sugars",
    summary="Detect and remove linear sugars",
    responses={
        200: {"description": "Successful response", "model": GetCircularSugarResponse},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def remove_circular_sugars(
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
):
    """Detect and remove circular sugars from a given SMILES string using Sugar.

    Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)

    Returns:
    - str: The modified SMILES string with circular sugars removed.
    """
    mol = parse_input(smiles, "cdk", False)
    try:
        removed_smiles = remove_circular_sugar(mol)
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
    operation_id="remove_linear_and_circular_sugars",
    summary="Detect and remove linear sugars",
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
async def remove_linear_and_circular_sugars(
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
):
    """Detect and remove linear and circular sugars from a given SMILES string.

    using Sugar Removal Utility.

    Parameters:
    - **SMILES string**: (str): SMILES: string representation of the molecule (required, query parameter)

    Returns:
    - str: The modified SMILES string with linear and circular sugars removed.
    """
    mol = parse_input(smiles, "cdk", False)
    try:
        removed_smiles = remove_linear_and_circular_sugar(mol)
        if removed_smiles:
            return removed_smiles
        else:
            raise HTTPException(
                status_code=422,
                detail="Error processing SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.post(
    "/extract_json",
    summary="Extract JSON data from PDF",
    response_description="Return extracted document structure as JSON",
    status_code=status.HTTP_200_OK,
)
async def extract_pdf_json(
    pdf_file: UploadFile = File(...),
    pages: int = Form(1, description="Number of pages to process"),
):
    """
    Upload a PDF file and extract its content as a structured JSON document.

    Args:
        pdf_file: The PDF file to process
        pages: Number of pages to process (default: 1)

    Returns:
        JSON: The extracted document structure
    """
    if not pdf_file.filename.lower().endswith(".pdf"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Uploaded file must be a PDF",
        )

    try:
        # Generate unique filename to prevent collisions
        unique_id = str(uuid.uuid4())[:8]
        original_filename = pdf_file.filename
        filename_parts = os.path.splitext(original_filename)
        safe_filename = f"{filename_parts[0]}_{unique_id}{filename_parts[1]}"

        # Create the full file path
        file_path = os.path.join(PDF_DIR, safe_filename)

        # Save the uploaded file
        with open(file_path, "wb") as buffer:
            content = await pdf_file.read()
            buffer.write(content)

        # Process the PDF file
        json_data = get_converted_document(file_path, number_of_pages=3)

        # Keep the file for future reference
        # (you can implement a cleanup strategy if needed)

        return json_data

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error processing PDF: {str(e)}",
        )


@router.post(
    "/extract_text",
    summary="Extract combined text from PDF",
    response_description="Return combined text extracted from PDF",
    status_code=status.HTTP_200_OK,
)
async def extract_pdf_text(
    pdf_file: UploadFile = File(...),
    pages: int = Form(1, description="Number of pages to process"),
):
    """
    Upload a PDF file and extract its content as a combined text paragraph.

    Args:
        pdf_file: The PDF file to process
        pages: Number of pages to process (default: 1)

    Returns:
        JSON: Object containing the combined text
    """
    if not pdf_file.filename.lower().endswith(".pdf"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Uploaded file must be a PDF",
        )

    try:
        # Use original filename without adding unique IDs
        original_filename = pdf_file.filename
        safe_filename = original_filename.replace(" ", "_")

        # Create the full file path
        file_path = os.path.join(PDF_DIR, safe_filename)

        # Check if file already exists
        if os.path.exists(file_path):
            # Use existing file
            pass
        else:
            # Save the uploaded file with original name
            with open(file_path, "wb") as buffer:
                content = await pdf_file.read()
                buffer.write(content)

        # Process the PDF file
        json_data = get_converted_document(file_path, number_of_pages=pages)
        result = extract_from_docling_document(json_data)
        combined_text = combine_to_paragraph(result)

        # Check if the extracted text has more than 10 words
        word_count = len(combined_text.split())
        if word_count <= 10:
            reader = PdfReader(file_path)
            if len(reader.pages) > 0:
                first_page_text = reader.pages[0].extract_text()
                combined_text = first_page_text.strip()

        return {"text": combined_text, "pdf_filename": safe_filename}

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error processing PDF: {str(e)}",
        )
