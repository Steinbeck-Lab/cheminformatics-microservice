import requests
import uuid
from fastapi.responses import JSONResponse
from urllib.request import urlopen
from urllib.parse import urlsplit
from fastapi import Body, APIRouter, status, HTTPException
from app.modules.decimer import getPredictedSegmentsFromFile
from app.schemas import HealthCheck
from app.schemas.error import ErrorResponse, BadRequestModel, NotFoundModel
from app.schemas.ocsr_schema import ExtractChemicalInfoResponse
from fastapi import UploadFile

router = APIRouter(
    prefix="/ocsr",
    tags=["ocsr"],
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
    summary="Perform a Health Check on OCSR Module",
    response_description="Return HTTP Status Code 200 (OK)",
    status_code=status.HTTP_200_OK,
    response_model=HealthCheck,
    include_in_schema=False,
)
def get_health() -> HealthCheck:
    """
    ## Perform a Health Check
    Endpoint to perform a health check on. This endpoint can primarily be used by Docker
    to ensure a robust container orchestration and management is in place. Other
    services that rely on the proper functioning of the API service will not deploy if this
    endpoint returns any other HTTP status code except 200 (OK).
    Returns:
        HealthCheck: Returns a JSON response with the health status
    """
    return HealthCheck(status="OK")


@router.post(
    "/process",
    summary="Detect, segment and convert a chemical structure depiction into a SMILES string using the DECIMER",
    responses={
        200: {
            "description": "Successful response",
            "model": ExtractChemicalInfoResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def Extract_ChemicalInfo_From_File(
    path: str = Body(
        None,
        embed=True,
        description="URL or local file path to the chemical structure depiction image.",
        openapi_examples={
            "example1": {
                "summary": "Cheminformatics - Article example image",
                "value": "https://static-content.springer.com/image/art%3A10.1186%2Fs13321-023-00744-6/MediaObjects/13321_2023_744_Figa_HTML.png",
            },
            "example2": {
                "summary": "Cheminformatics - Article example image",
                "value": "https://static-content.springer.com/image/art%3A10.1186%2Fs13321-023-00743-7/MediaObjects/13321_2023_743_Figa_HTML.png",
            },
        },
    ),
    reference: str = Body(None, embed=True, description="Reference information."),
    img: str = Body(
        None,
        embed=True,
        description="Bytes content of the chemical structure depiction image.",
    ),
):
    """
    Detect, segment and convert a chemical structure depiction into a SMILES string using the DECIMER modules.

    Parameters:
    - **path**: optional if img is provided (str): URL or local file path to the chemical structure depiction image.
    - **reference**: optional (str): URL or local file path to the chemical structure depiction image.
    - **img**: optional if a valid path is provided (str): URL or local file path to the chemical structure depiction image.

    Returns:
    - JSONResponse: A JSON response containing the extracted SMILES and the reference (if provided).

    Raises:
    - HTTPException: If the 'path' parameter is not provided or if it is an invalid URL or file path.
    - HTTPException: If the 'img' parameter is provided, but the content is not accessible.
    """

    if img:
        try:
            filename = "/tmp/" + str(uuid.uuid4())
            response = urlopen(img)
            smiles = getPredictedSegmentsFromFile(response.file.read(), filename)
            return JSONResponse(
                content={"reference": reference, "smiles": smiles.split(".")}
            )
        except Exception as e:
            raise HTTPException(
                status_code=400, detail="Error accessing image content: " + str(e)
            )
    else:
        try:
            split = urlsplit(path)
            filename = "/tmp/" + split.path.split("/")[-1]
            response = requests.get(path)
            if response.status_code == 200:
                smiles = getPredictedSegmentsFromFile(response.content, filename)
            return JSONResponse(
                content={"reference": reference, "smiles": smiles.split(".")}
            )
        except Exception as e:
            raise HTTPException(
                status_code=400, detail="Invalid URL or file path: " + str(e)
            )


@router.post(
    "/process-upload",
    summary="Detect, segment and convert a chemical structure depiction in the uploaded file into a SMILES string using the DECIMER",
    responses={
        200: {
            "description": "Successful response",
            "model": ExtractChemicalInfoResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def Extract_ChemicalInfo(file: UploadFile):
    """
    Detect, segment and convert a chemical structure depiction in the uploaded image file into a SMILES string using the DECIMER modules.

    Parameters:
    - **file**: required (File):

    Returns:
    - JSONResponse: A JSON response containing the extracted SMILES and the reference (if provided).

    Raises:
    - HTTPException: If the 'path' parameter is not provided or if it is an invalid URL or file path.
    - HTTPException: If the 'img' parameter is provided, but the URL is not accessible.
    """

    if file:
        filename = file.filename
        try:
            contents = file.file.read()
            try:
                smiles = getPredictedSegmentsFromFile(contents, filename)
                return JSONResponse(
                    content={"reference": None, "smiles": smiles.split(".")}
                )
            except Exception as e:
                raise HTTPException(
                    status_code=400, detail="Error accessing image URL: " + str(e)
                )
        except Exception as e:
            raise HTTPException(
                status_code=400, detail="Error accessing image URL: " + str(e)
            )
        finally:
            file.file.close()
