from pydantic import BaseModel


class ErrorResponse(BaseModel):
    """
    Represents an error response.

    Attributes:
        detail (str): A description of the error.
    """

    detail: str


class StandardizeRequest(BaseModel):
    """
    Represents a request to standardize a molecular structure.

    Attributes:
        molblock (str): The molecular structure in MDL Molfile (molblock) format.
    """

    molblock: str


class StandardizeResponse(BaseModel):
    """
    Represents a response after standardizing a molecular structure.

    Attributes:
        standardized_mol (str): The standardized molecular structure in MDL Molfile format.
        canonical_smiles (str): The canonical SMILES representation of the molecule.
        inchi (str): The InChI (International Chemical Identifier) of the molecule.
        inchikey (str): The InChIKey of the molecule.
    """

    standardized_mol: str
    canonical_smiles: str
    inchi: str
    inchikey: str


class SMILESValidationResult(BaseModel):
    """
    Represents the result of validating a SMILES string.

    Attributes:
        smi (str): The SMILES string.
        messages (tuple): A tuple of messages indicating validation results.
    """

    smi: str
    messages: tuple


class StandardizedResult(BaseModel):
    """
    Represents the original and standardized versions of a molecule.

    Attributes:
        original (SMILESValidationResult): The original SMILES validation result.
        standardized (SMILESValidationResult): The standardized SMILES validation result.
    """

    original: SMILESValidationResult
    standardized: SMILESValidationResult


class NPlikelinessScoreResponse(BaseModel):
    """
    Represents a response containing the NP likeness score of a molecule.

    Attributes:
        np_score (float): The NP likeness score of the molecule.
    """

    np_score: float
