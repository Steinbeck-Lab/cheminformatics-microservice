import os
import requests

from fastapi.responses import JSONResponse
from urllib.request import urlopen
from urllib.parse import urlsplit
from fastapi import Body, APIRouter
from typing_extensions import Annotated
from app.modules.decimer import getPredictedSegments

router = APIRouter(
    prefix="/ocsr",
    tags=["ocsr"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def ocsr_index():
    return {"module": "ocsr", "message": "Successful", "status": 200}


@router.post("/process")
async def Extract_ChemicalInfo(
    path: Annotated[str, Body(embed=True)] = None,
    reference: Annotated[str, Body(embed=True)] = None,
    img: Annotated[str, Body(embed=True)] = None,
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
        response = urlopen(img)
        with open(filename, "wb") as f:
            f.write(response.file.read())
            smiles = getPredictedSegments(filename)
            os.remove(filename)
            return JSONResponse(
                content={"reference": reference, "smiles": smiles.split(".")}
            )
    else:
        response = requests.get(path)
        if response.status_code == 200:
            with open(filename, "wb") as f:
                f.write(response.content)
                smiles = getPredictedSegments(filename)
                os.remove(filename)
                return JSONResponse(
                    content={"reference": reference, "smiles": smiles.split(".")}
                )
