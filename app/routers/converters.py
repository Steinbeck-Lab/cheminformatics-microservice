from fastapi import APIRouter

router = APIRouter(
    prefix="/converters",
    tags=["converters"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


@router.get("/")
async def converters_index():
    return {"module": "converters", "message": "Successful", "status": 200}
