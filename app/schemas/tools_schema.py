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
        output (str): A message indicating the type of sugars present in the molecule, either "The molecule contains Linear and Circular sugars", "The molecule contains only Linear sugar", "The molecule contains only Circular sugar", or "The molecule contains no sugar".
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
                    "input": "CC(=O)N[C@H]1[C@@H](O[C@@H]([C@H](O)[C@@H](O)C=O)[C@H](O)CO)O[C@H](CO)[C@H](O)[C@@H]1O",
                    "message": "Success",
                    "output": "The molecule contains Linear and Circular sugars",
                },
            ],
        }


class GetLinearSugarResponse(BaseModel):
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): The aglycone SMILES string or "No Linear sugar found".
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
                    "input": "CC1=CC(O)=C2C(=O)C3=C(OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)OC[C@@H](O)[C@@](O)(OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O)[C@H](O)[C@@H](O)C=O)C=CC=C3C(=O)C2=C1",
                    "message": "Success",
                    "output": "CC1=CC(=C2C(=C1)C(=O)C3=CC=CC(=C3C2=O)O)O",
                },
            ],
        }


class GetCircularSugarResponse(BaseModel):
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): The aglycone SMILES string or "No Circular sugar found".
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
                    "input": "C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
                    "message": "Success",
                    "output": "C=CC1C(C[C@H]2C3=C(CCN2)C4=C(C=CC=C4)N3)C(=CO[C@H]1O)C(=O)O",
                },
            ],
        }


class GetCircularandLinearSugarResponse(BaseModel):
    """A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): The aglycone SMILES string or "No Linear or Circular sugars found".
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
                    "input": "O=C(O)C1=C[C@@H](O)[C@@H](O)[C@H](OC(=O)[C@@H]2C(C(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)CO)=CC3=CC(O)=C(O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)C=C3[C@H]2C2=CC=C(O)C(O)=C2)C1",
                    "message": "Success",
                    "output": "C1=C(C=C(C(=C1)O)O)[C@@H]2C3=C(C=C([C@H]2C(=O)O[C@@H]4CC(=C[C@H]([C@H]4O)O)C(=O)O)C(=O)O)C=C(C(=C3)O)O",
                },
            ],
        }


class ExtractAglyconeAndSugarsResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): The SMILES representations of the aglycone and sugars as a printed list (["<SMILES>", "<SMILES>", ...]). The first one is always the aglycone. The list has a variable length dependening on how many sugar moieties were found.
    """

    message: str = "Success"
    output: List[str]

    class Config:
        """Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
                    "message": "Success",
                    "output": '["C=CC1C(C[C@H]2C3=C(CCN2)C4=C(C=CC=C4)N3)C(=CO[C@H]1O)C(=O)O","C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O)O1)O)O)O)O"]',
                },
            ],
        }


class GetAglyconeAndSugarIndicesResponse(BaseModel):
    """
    A Pydantic model representing a successful response.

    Attributes:
        message (str): A message indicating the success status (default: "Success").
        output (str): The atom indices of the aglycone and sugars in the given molecule as a printed list of indices lists. The first set of indices is always the aglycone. The list has a variable length dependening on how many sugar moieties were found.
    """

    message: str = "Success"
    output: List[List[int]]

    class Config:
        """Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CCCCC/C=C/C=C/[C@@H](O)C/C=C/C=C/C(=O)OC1C(O)[C@H](C2=C(O)C=C(O)C=C2CO)O[C@H](CO)[C@H]1O[C@@H]1OC(CO)[C@H](O)[C@H](O)C1O[C@@H]1OC(CO)[C@H](O)[C@H](O)C1O",
                    "message": "Success",
                    "output": "[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38], [38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60]]",
                },
            ],
        }
