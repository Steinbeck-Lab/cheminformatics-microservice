from __future__ import annotations

from typing import List

from pydantic import BaseModel


class GenerateStructuresResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (List[str]): A list of generated structures.

    """

    message: str = "Success"
    output: List[str]

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "C4H8",
                    "message": "Success",
                    "output": ["CC(C)C", "CCCC"],
                },
            ],
        }


class GetSugarInformationResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): Information on containing sugar

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
                {
                    "input": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
                    "message": "Success",
                    "output": "The molecule contains Linear and Circular sugars",
                },
            ],
        }


class GetLinearSugarResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SMILES without linear sugar

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
                {
                    "input": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
                    "message": "Success",
                    "output": "C(C1C(C(C(CO1)O)O)O)O",
                },
            ],
        }


class GetCircularSugarResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SMILES without circular sugar

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
                {
                    "input": "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
                    "message": "Success",
                    "output": "C(C(C(C(C(C1C(C(C(C(CO)O1)O)O)O)O)O)O)O)O",
                },
            ],
        }


class GetCircularandLinearSugarResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SMILES without circular and Linear sugar

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
                {
                    "input": "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
                    "message": "Success",
                    "output": "C1=C(C=C(C(=C1)O)O)C2C3=C(C=C(C=O)C2C(=O)OC4CC(=CC(C4O)O)C(=O)O)C=C(C(=C3)O)O",
                },
            ],
        }
