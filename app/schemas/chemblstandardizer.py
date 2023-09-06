from __future__ import annotations
from pydantic import BaseModel


class SMILESValidationResult(BaseModel):
    """
    Represents the result of validating a SMILES string.

    Attributes:
        smi (str): The SMILES string.
        messages (tuple): A tuple of messages indicating validation results.
    """

    smi: str
    messages: tuple


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
    """
    Represents a standardized result for a molecule.

    Attributes:
        standardized_mol (str): The standardized molecule representation.
        canonical_smiles (str): The canonical SMILES notation for the molecule.
        inchi (str): The InChI (International Chemical Identifier) for the molecule.
        inchikey (str): The InChIKey for the molecule.
    """

    standardized_mol: str
    canonical_smiles: str
    inchi: str
    inchikey: str
