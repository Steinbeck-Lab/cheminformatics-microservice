from fastapi import APIRouter, Depends, HTTPException

router = APIRouter(
    prefix="/converters",
    tags=["converters"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def chem_index():
    return {"message": "Hello Converters Router!"}
