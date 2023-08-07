from __future__ import annotations
from pydantic import BaseModel
from typing import List


class SMILESValidationResult(BaseModel):
    """
    Represents the result of validating a SMILES string.

    Attributes:
        smi (str): The SMILES string.
        messages (tuple): A tuple of messages indicating validation results.
    """

    smi: str
    messages: List[str]


class SMILESStandardizedResult(BaseModel):
    """
    Represents the original and standardized versions of a molecule.

    Attributes:
        original (SMILESValidationResult): The original SMILES validation result.
        standardized (SMILESValidationResult): The standardized SMILES validation result.
    """

    original: SMILESValidationResult
    standardized: SMILESValidationResult


class StandardizedResult(BaseModel):
    standardized_mol: str
    canonical_smiles: str
    inchi: str
    inchikey: str
