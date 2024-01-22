from __future__ import annotations

from pydantic import BaseModel


class ErrorResponse(BaseModel):
    """Represents an error response.

    Attributes:
        detail (str): A description of the error. This attribute contains a human-readable description of the error.
    """

    detail: str


class BadRequestModel(BaseModel):
    """Represents a Bad Request response in an API.

    Args:
        detail (str): A detailed message describing the nature of the bad request.
    """

    error: str = "Bad Request"
    detail: str


class NotFoundModel(BaseModel):
    """Represents a Not Found response in an API.

    Args:
        detail (str): A detailed message indicating that the requested resource was not found.
    """

    error: str = "Not Found"
    detail: str
