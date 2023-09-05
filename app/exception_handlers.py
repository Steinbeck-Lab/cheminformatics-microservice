from fastapi.responses import JSONResponse
from fastapi import Request


class InvalidInputException(Exception):
    def __init__(self, name: str, value: str):
        self.name = name
        self.value = value


async def input_exception_handler(request: Request, exc: InvalidInputException):
    return JSONResponse(
        status_code=422,
        content={"detail": f"Error reading {exc.name}, check again: {exc.value}"},
    )
