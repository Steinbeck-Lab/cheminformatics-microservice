# Coding Conventions

**Analysis Date:** 2026-03-12

## Naming Patterns

**Files:**
- Python modules: `lowercase_with_underscores.py` (e.g., `rdkit_wrapper.py`, `cdk_wrapper.py`)
- Router files: `noun_plural.py` (e.g., `chem.py`, `converters.py`, `depict.py`, `tools.py`)
- Test files: `test_*.py` or `*_test.py` prefixed with `test_` (e.g., `test_chem.py`, `test_converters.py`)
- Schema files: `noun_schema.py` or descriptive names (e.g., `chem_schema.py`, `error.py`, `classyfire.py`)

**Functions:**
- Python functions: `snake_case` (e.g., `parse_input()`, `get_tanimoto_similarity()`, `check_RO5_violations()`)
- Private/internal functions: same convention with leading underscore if needed
- Wrapper/utility functions named by action: `get_*()`, `check_*()`, `parse_*()`, `ensure_*()` (e.g., `get_CDK_HOSE_codes()`, `ensure_2d()`)

**Variables:**
- Local/module variables: `snake_case` (e.g., `test_smiles`, `logger`, `client`)
- Constants: `UPPER_SNAKE_CASE` (e.g., `DEFAULT_TIMEOUT`, `MAX_RETRIES`)
- Instance attributes: `snake_case` (e.g., `self.timeout`, `self.cache_size`)

**Types:**
- Classes: `PascalCase` (e.g., `InvalidInputException`, `PubChemClient`, `FilterCatalog`)
- Type hints: Python 3.10+ syntax using `from __future__ import annotations` at top of file
- Enums/Literals: Used for constrained values (e.g., `framework: Literal["rdkit", "cdk", "openbabel"]`)

**API/Routes:**
- Endpoint paths: `/lowercase/with/slashes` or `/snake_case_parameters` (e.g., `/chem/`, `/depict/2D`, `/convert/mol2D`)
- Query parameters: `camelCase` (e.g., `?smiles=`, `?toolkit=`, `?framework=`, `?standardize=`)

## Code Style

**Formatting:**
- Formatter: `black` (v22.6.0)
- Line length: No strict limit enforced (E501 ignored by flake8)
- Indentation: 4 spaces (Python standard)

**Linting:**
- Tool: `flake8` (v7.0.0)
- Config location: `.pre-commit-config.yaml`
- Ignores: `E402` (module-level import not at top), `E501` (line length), `W503` (line break before binary operator), `E203` (whitespace before colon)
- Per-file ignores: `__init__.py:F401` (unused imports in init files allowed)

**Docstrings:**
- Tool: `pydocstringformatter` (v0.7.3)
- Format: Google-style docstrings with `Args:`, `Returns:`, `Raises:` sections
- Example from `rdkit_wrapper.py`:
```python
def check_RO5_violations(molecule: any) -> int:
    """Check the molecule for violations of Lipinski's Rule of Five.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        int: Number of Lipinski Rule violations.
    """
```

## Import Organization

**Order:**
1. Future imports: `from __future__ import annotations`
2. Standard library: `import os`, `import io`, `from typing import ...`
3. Third-party imports: `from rdkit import Chem`, `import pandas as pd`
4. Local/app imports: `from app.modules import ...`, `from app.routers import ...`

**Tool:** `reorder-python-imports` (v3.12.0) enforces import order in pre-commit hooks

**Path Aliases:**
- Absolute imports used throughout (e.g., `from app.routers import chem`)
- No relative imports (`from . import`) observed in production code
- App root is `app/` directory

## Error Handling

**Custom Exceptions:**
- Location: `app/exception_handlers.py`
- `InvalidInputException(name: str, value: str)` - raised when input parsing fails
- Imported as: `from app.exception_handlers import InvalidInputException`
- Raised in: `app/modules/toolkits/helpers.py:parse_SMILES()` when parsing fails

**Exception Pattern:**
```python
try:
    # parsing logic
except Exception:
    raise InvalidInputException(name="smiles", value=smiles)
```

**HTTP Response Mapping:**
- `InvalidInputException` → HTTP 422 (Unprocessable Entity)
- Generic exceptions → HTTP 500 (handled via FastAPI error handlers)
- Validation errors → HTTP 400 or 422 depending on context

**Logging:**
- Module-level logger: `logger = logging.getLogger(__name__)`
- Log levels used: `logger.exception()` for error cases, `logger.debug()` for debug info
- Example: `logger.exception("Error processing SMILES")`

## Logging

**Framework:** Python's `logging` module (standard library)

**Patterns:**
- Module-level logger instantiation: `logger = logging.getLogger(__name__)` (one per module)
- Exception logging: `logger.exception("Description")` - logs full traceback
- Debug logging: `logger.debug("message: %s", variable)` - for verbose output
- Used in: routers (`chem.py`, `tools.py`) and critical modules

**When to Log:**
- Exceptions: Always use `logger.exception()` at error boundaries
- Debug info: Use `logger.debug()` for intermediate values (e.g., `logger.debug("Sugar removal result: %s", removed_smiles)`)
- Do NOT log errors that are caught and handled gracefully

## Comments

**When to Comment:**
- Complex algorithms or non-obvious logic (e.g., ChEMBL standardization pipeline)
- Integration points between toolkits (RDKit, CDK, OpenBabel fallbacks)
- Workarounds for library limitations

**JSDoc/TSDoc:**
- Python: Use docstrings, no JSDoc
- TypeScript/JavaScript: Not widely used in frontend (React components documented via PropTypes/TypeScript)
- Docstring format: Google-style with `Args:`, `Returns:`, `Raises:` sections

**Comment Style:**
- Descriptive comments explain "why" not "what" (code is self-documenting)
- Example from `test_converters.py`: `# Check if SMILES contains at least carbon`
- Commented-out code: Acceptable temporarily but should not persist (e.g., lines 103-108 in `test_converters.py`)

## Function Design

**Size:**
- Functions kept relatively focused (most under 50 lines)
- Toolkit wrappers encapsulate specific operations (e.g., `get_tanimoto_similarity()`, `check_RO5_violations()`)

**Parameters:**
- Named parameters preferred over positional
- Type hints required for all parameters (Python 3.10+ with `from __future__ import annotations`)
- Optional parameters use `Optional[Type]` with defaults
- Example: `def parse_input(input: str, framework: str = "rdkit", standardize: bool = False):`

**Return Values:**
- Explicit return types required
- Can return `None` for optional outcomes
- Toolkit functions return toolkit-native objects (RDKit `Mol`, CDK `IAtomContainer`)
- API endpoints return Pydantic models (`BaseModel` subclasses)

## Module Design

**Exports:**
- Public functions/classes exported directly from modules
- No explicit `__all__` observed; all public functions importable
- Modules organized by functionality: `toolkits/`, `modules/`, `routers/`, `schemas/`

**Barrel Files:**
- Not used extensively (no `__init__.py` aggregating exports)
- `__init__.py` files mostly empty (F401 ignored)
- Import directly from source modules: `from app.routers import chem`

**Module Structure Pattern:**
- Wrapper modules (e.g., `rdkit_wrapper.py`, `cdk_wrapper.py`) export utility functions
- Helper modules (e.g., `helpers.py`) export dispatcher/parsing logic
- Router modules (e.g., `chem.py`) export FastAPI `APIRouter` instance and endpoint handlers
- Schema modules (e.g., `chem_schema.py`) export Pydantic `BaseModel` subclasses

---

*Convention analysis: 2026-03-12*
