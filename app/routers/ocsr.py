import os
import requests

from fastapi.responses import JSONResponse
from urllib.request import urlopen
from urllib.parse import urlsplit
from fastapi import Body, APIRouter, status, HTTPException
from typing import Dict
from app.modules.decimer import getPredictedSegments
from app.schemas import HealthCheck
from app.schemas.pydanticmodels import ErrorResponse

router = APIRouter(
    prefix="/ocsr",
    tags=["ocsr"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
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
    Endpoint to perform a healthcheck on. This endpoint can primarily be used Docker
    to ensure a robust container orchestration and management is in place. Other
    services which rely on proper functioning of the API service will not deploy if this
    endpoint returns any other HTTP status code except 200 (OK).
    Returns:
        HealthCheck: Returns a JSON response with the health status
    """
    return HealthCheck(status="OK")


@router.post(
    "/process",
    response_model=Dict,
    summary="Detect, segment and convert a chemical structure depiction into a SMILES string using the DECIMER",
    responses={400: {"model": ErrorResponse}},
)
async def Extract_ChemicalInfo(
    path: str = Body(
        None,
        embed=True,
        description="URL or local file path to the chemical structure depiction image.",
    ),
    reference: str = Body(None, embed=True, description="Reference information."),
    img: str = Body(
        None, embed=True, description="URL of the chemical structure depiction image."
    ),
):
    """
    Detect, segment and convert a chemical structure depiction into a SMILES string using the DECIMER modules.

    Parameters:
    - **Images**: required (str): URL or local file path to the chemical structure depiction image.

    Returns:
    - JSONResponse: A JSON response containing the extracted SMILES and the reference (if provided).

    Raises:
    - HTTPException: If the 'path' parameter is not provided or if it is an invalid URL or file path.
    - HTTPException: If the 'img' parameter is provided, but the URL is not accessible.
    """

    split = urlsplit(path)
    filename = "/tmp/" + split.path.split("/")[-1]
    if img:
        try:
            response = urlopen(img)
            with open(filename, "wb") as f:
                f.write(response.file.read())
                smiles = getPredictedSegments(filename)
                os.remove(filename)
                return JSONResponse(
                    content={"reference": reference, "smiles": smiles.split(".")}
                )
        except Exception as e:
            raise HTTPException(
                status_code=400, detail="Error accessing image URL: " + str(e)
            )
    else:
        try:
            response = requests.get(path)
            if response.status_code == 200:
                with open(filename, "wb") as f:
                    f.write(response.content)
                    smiles = getPredictedSegments(filename)
                    os.remove(filename)
                    return JSONResponse(
                        content={"reference": reference, "smiles": smiles.split(".")}
                    )
            else:
                raise HTTPException(
                    status_code=400, detail="Error accessing provided URL"
                )
        except Exception as e:
            raise HTTPException(
                status_code=400, detail="Invalid URL or file path: " + str(e)
            )
