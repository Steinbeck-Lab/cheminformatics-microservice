from pydantic import BaseModel


class ErrorResponse(BaseModel):
    """
    Represents an error response.

    Attributes:
        detail (str): A description of the error. This attribute contains a human-readable description of the error.
    """

    detail: str
