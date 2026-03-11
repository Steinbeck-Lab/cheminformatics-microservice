from __future__ import annotations

import ipaddress
import logging
import uuid
from typing import Annotated
from urllib.parse import urlsplit
from urllib.request import urlopen

import httpx
from fastapi import APIRouter
from fastapi import Body
from fastapi import File
from fastapi import HTTPException
from fastapi import status
from fastapi import UploadFile
from fastapi.responses import JSONResponse

from app.modules.decimer import get_predicted_segments_from_file
from app.schemas import HealthCheck
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse
from app.schemas.error import NotFoundModel
from app.schemas.ocsr_schema import ExtractChemicalInfoResponse

logger = logging.getLogger(__name__)

ALLOWED_SCHEMES = {"http", "https"}


def _validate_url(url: str) -> None:
    """Validate that a URL is safe to fetch (no SSRF to internal networks)."""
    parsed = urlsplit(url)
    if parsed.scheme not in ALLOWED_SCHEMES:
        raise ValueError(f"URL scheme '{parsed.scheme}' is not allowed")
    hostname = parsed.hostname
    if not hostname:
        raise ValueError("URL has no hostname")
    try:
        addr = ipaddress.ip_address(hostname)
        if addr.is_private or addr.is_loopback or addr.is_link_local or addr.is_reserved:
            raise ValueError("URLs pointing to internal/private networks are not allowed")
    except ValueError as exc:
        if "internal" in str(exc) or "not allowed" in str(exc):
            raise
        # hostname is a DNS name — block known internal names
        lower = hostname.lower()
        if lower in ("localhost",) or lower.endswith(".local"):
            raise ValueError(
                "URLs pointing to internal/private networks are not allowed"
            )


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
    """## Perform a Health Check.

    Endpoint to perform a health check on. This endpoint can primarily be used by Docker
    to ensure a robust container orchestration and management are in place. Other
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
def Extract_ChemicalInfo_From_File(
    path: str = Body(
        None,
        embed=True,
        description="Local or Remote path to the image file",
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
    reference: str = Body(
        None,
        embed=True,
        description="User-defined reference information for tracking",
    ),
    img: str = Body(
        None,
        embed=True,
        description="Image: Bytes content of the chemical structure depiction image",
    ),
    hand_drawn: bool = Body(
        False,
        embed=True,
        description="Use hand-drawn model for prediction",
    ),
):
    """Detect, segment and convert a chemical structure depiction into a SMILES.

    string using the DECIMER modules.

    Parameters:
    - **path**: optional if img is provided (str): Local or Remote path to the image file.
    - **reference**: optional (str): User-defined reference information for tracking.
    - **img**: optional if a valid path is provided (str): Image: Bytes content of the chemical structure depiction image.
    - **hand_drawn**: optional (bool): Use hand-drawn model for prediction. Defaults to False.

    Returns:
    - JSONResponse: A JSON response containing the extracted SMILES and the reference (if provided).

    Raises:
    - HTTPException: If the 'path' parameter is not provided or if it is an invalid URL or file path.
    - HTTPException: If the 'img' parameter is provided, but the content is not accessible.
    """

    if img:
        try:
            _validate_url(img)
            filename = "/tmp/" + str(uuid.uuid4())
            response = urlopen(img)
            smiles = get_predicted_segments_from_file(
                response.file.read(),
                filename,
                hand_drawn=hand_drawn,
            )
            return JSONResponse(
                content={"reference": reference, "smiles": smiles.split(".")},
            )
        except HTTPException:
            raise
        except Exception:
            logger.exception("Error processing image from img parameter")
            raise HTTPException(
                status_code=422,
                detail="Error processing the image",
            )
    else:
        try:
            _validate_url(path)
            filename = "/tmp/" + str(uuid.uuid4())
            response = httpx.get(path)
            if response.status_code == 200:
                smiles = get_predicted_segments_from_file(
                    response.content,
                    filename,
                    hand_drawn=hand_drawn,
                )
            return JSONResponse(
                content={"reference": reference, "smiles": smiles.split(".")},
            )
        except HTTPException:
            raise
        except (httpx.HTTPError, ValueError, IOError):
            logger.exception("Error processing image from URL path")
            raise HTTPException(
                status_code=422,
                detail="Error processing the image",
            )
        except Exception:
            logger.exception("Unexpected error processing image from URL path")
            raise HTTPException(
                status_code=422,
                detail="Error processing the image",
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
def extract_chemicalinfo_from_upload(
    file: Annotated[UploadFile, File(description="Chemical structure depiction image")],
    hand_drawn: bool = Body(False, description="Use hand-drawn model for prediction"),
):
    """Detect, segment and convert a chemical structure depiction in the.

    uploaded image file into a SMILES string using the DECIMER modules.

    Parameters:
    - **file**: required (File): Chemical structure depiction image
    - **hand_drawn**: optional (bool): Use hand-drawn model for prediction. Defaults to False.

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
                smiles = get_predicted_segments_from_file(
                    contents, filename, hand_drawn=hand_drawn
                )
                return JSONResponse(
                    content={"reference": None, "smiles": smiles.split(".")},
                )
            except Exception:
                logger.exception("Error processing uploaded image")
                raise HTTPException(
                    status_code=422,
                    detail="Error processing the image",
                )
        except HTTPException:
            raise
        except Exception:
            logger.exception("Error reading uploaded image file")
            raise HTTPException(
                status_code=422,
                detail="Error processing the image",
            )
        finally:
            file.file.close()
