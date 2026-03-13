from __future__ import annotations

import logging

from fastapi import Request
from fastapi.responses import JSONResponse

logger = logging.getLogger(__name__)


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
    logger.warning("Invalid input for %s: %s", exc.name, exc.value)
    return JSONResponse(
        status_code=422,
        content={"detail": f"Invalid input for '{exc.name}'"},
    )
