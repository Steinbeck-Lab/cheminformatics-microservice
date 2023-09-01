from pydantic import BaseModel, Field
from typing import Dict, Any, List


class GenerateStereoisomersResponse(BaseModel):
    """
    Represents a response containing enumerated stereo isomers for a molecule.

    Properties:
    - stereoisomers (List[str]): A list of stereoisomer SMILES strings.
    """

    stereoisomers: list = Field(
        ...,
        title="Stereoisomers",
        description="A list of stereoisomer SMILES strings.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "message": "Success",
                    "output": '["Cn1c(=O)c2c(ncn2C)n(C)c1=O"]',
                }
            ]
        }


class GenerateDescriptorsResponse(BaseModel):
    """
    Represents a successful response containing descriptors for a molecule.

    Properties:
    - descriptors (dict): A dictionary of descriptors and their values.
    """

    descriptors: dict = Field(
        ...,
        title="Descriptors",
        description="A dictionary of descriptors and their values.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "message": "Success",
                    "output": "a defined set of moleculare descriptors calculated",
                }
            ]
        }


class GenerateMultipleDescriptorsResponse(BaseModel):
    """
    Represents a response containing multiple descriptors for a list of SMILES strings.

    Properties:
    - descriptors (dict): A dictionary with each SMILES as the key and the corresponding descriptors as the value.
    """

    descriptors: Dict[str, Any] = Field(
        ...,
        title="Descriptors",
        description="A dictionary with each SMILES as the key and the corresponding descriptors as the value.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C, CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
                    "message": "Success",
                    "output": '{"CN1C=NC2=C1C(=O)N(C(=O)N2C)C": {"A set of calculated molecular descriptors"},"CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1": {"A set of calculated molecular descriptors"}}',
                }
            ]
        }


class GenerateHOSECodeResponse(BaseModel):
    """
    Represents a response containing HOSE codes for a molecule.

    Properties:
    - hose_codes (List[str]): A list of HOSE codes for each atom in the molecule.
    """

    hose_codes: List[str] = Field(
        ...,
        title="HOSE Codes",
        description="A list of HOSE codes for each atom in the molecule.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "message": "Success",
                    "output": '[  "C-4;N(//)",  "N-3;*C*CC(//)",  "C-3;*N*N(//)",  "N-2;*C*C(//)",  "C-3;*C*N*N(//)",  "C-3;*C*C*N(//)",  "C-3;=O*C*N(//)",  "O-1;=C(//)",  "N-3;*C*CC(//)",  "C-3;=O*N*N(//)",  "O-1;=C(//)",  "N-3;*C*CC(//)",  "C-4;N(//)"]',
                }
            ]
        }


class GenerateStandardizeResponse(BaseModel):
    """
    Represents a response containing standardized molecule information.

    Properties:
    - standardized_mol (str): The standardized molblock of the molecule.
    - canonical_smiles (str): The canonical SMILES representation of the molecule.
    - inchi (str): The InChI representation of the molecule.
    - inchikey (str): The InChI-Key of the molecule.
    """

    standardized_mol: str = Field(
        ...,
        title="Standardized Molblock",
        description="The standardized molblock of the molecule.",
    )
    canonical_smiles: str = Field(
        ...,
        title="Canonical SMILES",
        description="The canonical SMILES representation of the molecule.",
    )
    inchi: str = Field(
        ...,
        title="InChI",
        description="The InChI representation of the molecule.",
    )
    inchikey: str = Field(
        ...,
        title="InChI-Key",
        description="The InChI-Key of the molecule.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": """
          CDK     09012310592D

          1  0  0  0  0  0  0  0  0  0999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        M  END""",
                    "message": "Success",
                    "output": "Standardized Molblock, Canonical SMILES, InChI, InChI-Key",
                }
            ]
        }


class NPlikelinessScoreResponse(BaseModel):
    """
    Represents a response containing the natural product likeness score.

    Properties:
    - np_score (float): The calculated natural product likeness score.
    """

    np_score: float = Field(
        ...,
        title="Natural Product Likeness Score",
        description="The calculated natural product likeness score.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "message": "Success",
                    "output": "-1.09",
                }
            ]
        }


class TanimotoSimilarityResponse(BaseModel):
    """
    Represents the Tanimoto similarity index for a pair of SMILES strings.

    Properties:
    - similarity (float): The Tanimoto similarity index.
    """

    similarity: float = Field(
        ...,
        title="Tanimoto Similarity",
        description="The Tanimoto similarity index as a floating-point value.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1,CC",
                    "message": "Success",
                    "output": "0.024390243902439025",
                }
            ]
        }


class TanimotoMatrixResponse(BaseModel):
    """
    Response model for Tanimoto similarity matrix.

    Attributes:
    - similarity_matrix (List[List[float]]): A 2D list representing the Tanimoto similarity matrix.
      Each inner list corresponds to a row in the matrix, and each value inside the inner list
      represents the similarity score between two molecules.
    """

    similarity_matrix: list = Field(
        ...,
        title="Tanimoto Similarities",
        description="The Tanimoto similarities as a 2D list representing the Tanimoto similarity matrix.",
    )

    class Config:
        """
        Pydantic model configuration.

        JSON Schema Extra:
        - Includes examples of the response structure.
        """

        json_schema_extra = {
            "examples": [
                {
                    "input": "CCC,CC,CCC",
                    "message": "Success",
                    "output": "[[1.0, 0.2, 1.0], [0.2, 1.0, 0.2], [1.0, 0.2, 1.0]]",
                }
            ]
        }
