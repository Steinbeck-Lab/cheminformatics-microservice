from __future__ import annotations

from typing import List, Dict, Union

from pydantic import BaseModel


class GenerateStructuresResponse(BaseModel):
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (Dict): A dictionary containing structure generation results with:
            - total_count: Total number of possible structures
            - generated_count: Number of structures actually generated
            - structures: List of SMILES strings (limited)
            - settings: Dictionary describing the surge settings used
            - formula: The input molecular formula
            - limit_applied: Whether a limit was applied to results
    """

    message: str = "Success"
    output: Dict[str, Union[int, List[str], Dict[str, str], str, bool]]

    class Config:
        """Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "message": "Success",
                    "output": {
                        "total_count": 24000,
                        "generated_count": 1000,
                        "structures": ["CC(C)C", "CCCC"],
                        "settings": {
                            "-P": "Require planarity",
                            "-T": "Disallow triple bonds",
                            "-B1,2,3,4,5,7,9": "Avoid substructures: no triple bonds in small rings, Bredt's rule violations, cumulative double bonds, forbidden topologies",
                            "-t0": "No rings of length 3 allowed",
                            "-f0": "No cycles of length 4 allowed",
                            "-S": "Output in SMILES format",
                        },
                        "formula": "C10H16",
                        "limit_applied": True,
                    },
                },
            ],
        }


class GetSugarInformationResponse(BaseModel):
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): Information on containing sugar
    """

    message: str = "Success"
    output: str

    class Config:
        """Pydantic model configuration.

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
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SMILES without linear sugar
    """

    message: str = "Success"
    output: str

    class Config:
        """Pydantic model configuration.

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
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SMILES without circular sugar
    """

    message: str = "Success"
    output: str

    class Config:
        """Pydantic model configuration.

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
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): SMILES without circular and Linear sugar
    """

    message: str = "Success"
    output: str

    class Config:
        """Pydantic model configuration.

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
