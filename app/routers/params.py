from __future__ import annotations

from fastapi import Query

AUTO_DETECT_PARAM = Query(
    default=False,
    description=(
        "When true, auto-detect input format (SMILES, InChI, molblock, etc.) "
        "instead of assuming SMILES."
    ),
)
