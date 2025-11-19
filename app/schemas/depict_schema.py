from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class Depict2DResponse(BaseModel):
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SVG code block of the depicted image
    """

    message: str = "Success"
    output: str

    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
                {"input": "CCCOC", "message": "Success", "output": "SVG string"},
            ],
        }
    )


class Depict3DResponse(BaseModel):
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): HTML code block of the depicted image
    """

    message: str = "Success"
    output: str

    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
                {"input": "CCCOC", "message": "Success", "output": "HTML string"},
            ],
        }
    )
