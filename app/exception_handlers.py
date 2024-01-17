from __future__ import annotations

from fastapi import Request
from fastapi.responses import JSONResponse


class InvalidInputException(Exception):
    def __init__(self, name: str, value: str):
        self.name = name
        self.value = value


async def input_exception_handler(request: Request, exc: InvalidInputException):
    """Custom exception handler for InvalidInputException.

    Args:
        request (Request): The FastAPI Request object.
        exc (InvalidInputException): The InvalidInputException instance.

    Returns:
        JSONResponse: A JSON response containing error details.
    """
    return JSONResponse(
        status_code=422,
        content={"detail": f"Error reading {exc.name}, check again: {exc.value}"},
    )
