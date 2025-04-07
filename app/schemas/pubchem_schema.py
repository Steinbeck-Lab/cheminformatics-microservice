from __future__ import annotations
from typing import Optional

from pydantic import BaseModel


class PubChemResponse(BaseModel):
    """Response model for the SMILES endpoint."""

    input: str
    canonical_smiles: Optional[str] = None
    input_type: Optional[str] = None
    success: bool
