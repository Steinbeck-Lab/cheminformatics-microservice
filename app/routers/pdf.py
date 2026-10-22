from __future__ import annotations

import logging
import os

from fastapi import APIRouter
from fastapi import File
from fastapi import Form
from fastapi import HTTPException
from fastapi import Request
from fastapi import status
from fastapi import UploadFile

from app.limiter import limiter
from app.modules.docling_wrapper import ConverterBusyError
from app.modules.docling_wrapper import extract_text
from app.modules.docling_wrapper import get_converted_document
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse

logger = logging.getLogger(__name__)

MAX_PDF_SIZE = 10 * 1024 * 1024
MAX_PDF_PAGES = 20

PDF_FILE_PARAM = File(
    ..., description=f"PDF file, {MAX_PDF_SIZE // (1024 * 1024)} MB or smaller"
)
PAGES_PARAM = Form(1, ge=1, le=MAX_PDF_PAGES, description="Number of pages to process")

# Same /tools prefix as the tools router; kept apart so main.py can leave it out
# of the lite image, which ships without docling.
router = APIRouter(
    prefix="/tools",
    tags=["tools"],
    responses={
        400: {"description": "Bad Request", "model": BadRequestModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)


def _convert_upload(pdf_file: UploadFile, pages: int) -> dict:
    """Validate an uploaded PDF and convert its first pages with Docling."""
    content = pdf_file.file.read(MAX_PDF_SIZE + 1)
    if len(content) > MAX_PDF_SIZE:
        raise HTTPException(
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            detail=f"PDF must be {MAX_PDF_SIZE // (1024 * 1024)} MB or smaller",
        )
    if not content.startswith(b"%PDF-"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Uploaded file must be a PDF",
        )
    try:
        return get_converted_document(content, number_of_pages=pages)
    except ConverterBusyError:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="PDF converter is busy, please retry later",
        )
    except Exception:
        logger.exception("Error converting PDF")
        raise HTTPException(status_code=422, detail="Error processing the PDF")


@router.post(
    "/extract_json",
    summary="Extract the document structure of a PDF as JSON using Docling",
    response_description="Return the Docling document as JSON",
    status_code=status.HTTP_200_OK,
)
@limiter.limit("5/minute")
def extract_pdf_json(
    request: Request,
    pdf_file: UploadFile = PDF_FILE_PARAM,
    pages: int = PAGES_PARAM,
):
    """## Extract the document structure of a PDF using Docling.

    Parameters:
    - **pdf_file**: required (File): The PDF file to process.
    - **pages**: optional (int): Number of pages to process from the start. Defaults to 1.

    Returns:
    - JSON: The Docling document structure.
    """
    return _convert_upload(pdf_file, pages)


@router.post(
    "/extract_text",
    summary="Extract the title, abstract and introduction of a PDF as text",
    response_description="Return the combined text extracted from the PDF",
    status_code=status.HTTP_200_OK,
)
@limiter.limit("5/minute")
def extract_pdf_text(
    request: Request,
    pdf_file: UploadFile = PDF_FILE_PARAM,
    pages: int = PAGES_PARAM,
):
    """## Extract the title, abstract and introduction of a PDF as one paragraph.

    Falls back to all text on the first page when 10 words or fewer are found.

    Parameters:
    - **pdf_file**: required (File): The PDF file to process.
    - **pages**: optional (int): Number of pages to process from the start. Defaults to 1.

    Returns:
    - JSON: The combined text and the name of the uploaded file.
    """
    doc = _convert_upload(pdf_file, pages)
    filename = os.path.basename(pdf_file.filename or "").replace(" ", "_")
    return {"text": extract_text(doc), "pdf_filename": filename}
