import os
import requests

from fastapi.responses import JSONResponse
from urllib.request import urlopen
from urllib.parse import urlsplit
from fastapi import Body, APIRouter
from typing_extensions import Annotated
from app.modules.decimermodules import getPredictedSegments

router = APIRouter(
    prefix="/decimer",
    tags=["decimer"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def decimer_index():
    return {"module": "decimer", "message": "Successful", "status": 200}


@router.post("/process")
async def extract_chemicalinfo(path: Annotated[str, Body(embed=True)], reference: Annotated[str, Body(embed=True)], img: Annotated[str, Body(embed=True)]):
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
