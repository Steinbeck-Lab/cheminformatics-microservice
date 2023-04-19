import os
import requests

from fastapi.responses import JSONResponse
from urllib.request import urlopen
from urllib.parse import urlsplit
from fastapi import Request, APIRouter
from app.modules.decimermodules import getPredictedSegments

router = APIRouter(
    prefix="/decimer",
    tags=["decimer"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def DECIMER_Index():
    return {"module": "decimer", "message": "Successful", "status": 200}


@router.post("/process")
async def Extract_ChemicalInfo(request: Request):
    """
    Extract chemical structure depictions and convert them into SMILES using DECIMER:

    - **Images**: required (query)
    """
    body = await request.json()
    image_path = body["path"]
    reference = body["reference"]
    split = urlsplit(image_path)
    filename = "/tmp/" + split.path.split("/")[-1]
    if "img" in body:
        imgDataURI = body["img"]
        if imgDataURI:
            response = urlopen(imgDataURI)
            with open(filename, "wb") as f:
                f.write(response.file.read())
                smiles = getPredictedSegments(filename)
                os.remove(filename)
                return JSONResponse(
                    content={"reference": reference, "smiles": smiles.split(".")}
                )
    else:
        response = requests.get(image_path)
        if response.status_code == 200:
            with open(filename, "wb") as f:
                f.write(response.content)
                smiles = getPredictedSegments(filename)
                os.remove(filename)
                return JSONResponse(
                    content={"reference": reference, "smiles": smiles.split(".")}
                )
