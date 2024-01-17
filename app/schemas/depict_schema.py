from __future__ import annotations

from pydantic import BaseModel


class Depict2DResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SVG code block of the depicted image
    """

    message: str = "Success"
    output: str

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {"input": "CCCOC", "message": "Success", "output": "SVG string"},
            ],
        }


class Depict3DResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): HTML code block of the depicted image
    """

    message: str = "Success"
    output: str

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {"input": "CCCOC", "message": "Success", "output": "HTML string"},
            ],
        }
