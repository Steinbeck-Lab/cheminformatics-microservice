# Testing Patterns

**Analysis Date:** 2026-03-12

## Test Framework

**Runner:**
- Framework: `pytest` (latest version)
- Config: No explicit `pytest.ini` or `pyproject.toml` config; uses pytest defaults
- CI config: `.github/workflows/test.yml` runs tests on Python 3.10

**Assertion Library:**
- pytest's built-in assertions used throughout
- No separate assertion library imported

**Run Commands:**
```bash
python3 -m pytest                              # Run all tests
python3 -m pytest tests/test_chem.py           # Run a single test file
python3 -m pytest tests/test_chem.py::test_name # Run a single test function
python3 -m pytest --cov=./ --cov-report=xml    # With coverage (CI command)
python3 -m pytest tests/ -v                     # Verbose output
```

## Test File Organization

**Location:**
- All tests in `tests/` directory at project root
- Co-located pattern: One `tests/` directory, NOT per-module (all tests together)

**Naming:**
- Pattern: `test_*.py` (e.g., `test_chem.py`, `test_converters.py`, `test_depict.py`)
- Multiple test files per router: `test_cdk_*.py` for CDK-specific depiction tests (e.g., `test_cdk_annotations.py`, `test_cdk_dative_bonds.py`)
- Fixture files: included in test modules themselves (no separate `conftest.py`)

**Structure:**
```
tests/
├── test_chem.py                      # Main chemistry endpoints
├── test_converters.py                # Format conversion endpoints
├── test_depict.py                    # 2D/3D depiction endpoints
├── test_tools.py                     # Tools/utilities endpoints
├── test_ocsr.py                      # OCSR (Optical recognition) endpoints
├── test_helper.py                    # Helper function tests
├── test_invalid_functions.py         # Error handling tests
├── test_cdk_*.py                     # CDK-specific depiction features (15+ files)
├── test_filters.py                   # Filter/sanitization tests
├── test_cdk_cxsmiles_parser.py      # CXSMILES parsing tests
├── caffeine.png                      # Test fixtures (images, molecules)
├── segment_sample.png
└── ketcher.cdx
```

## Test Structure

**Suite Organization:**
```python
from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.mark.parametrize(
    "smiles, response_text, response_code",
    [
        ("CC", '["CC"]', 200),
        ("INVALID_INPUT", "", 422),
    ],
)
def test_smiles_to_stereo_isomers(smiles, response_text, response_code):
    response = client.get(f"/latest/chem/stereoisomers?smiles={smiles}")
    assert response.status_code == response_code
    if smiles != "INVALID_INPUT":
        assert response.text == response_text
```

**Patterns:**

1. **Fixtures:**
   - Test-level fixtures (no module or session scope by default)
   - Named descriptively: `test_smiles`, `invalid_smiles`, `molfile`
   - Decorator: `@pytest.fixture`
   - Autouse fixtures for teardown (e.g., `_reset_limiter()` in `test_converters.py`)

2. **Parametrization:**
   - Heavy use of `@pytest.mark.parametrize()` for matrix testing
   - Parameters: `(input, expected_output, response_code)`
   - Format: tuple list with values for each parameter set
   - Used to test multiple toolkits, input formats, and edge cases

3. **Assertions:**
   - Status codes: `assert response.status_code == expected_code`
   - Response bodies: `assert response.json() == expected_dict`
   - Text comparison: `assert response.text == expected_text`
   - Headers: `assert response.headers["content-type"] == "application/json"`
   - List membership: `assert data["summary"]["successful"] >= 1`

4. **Test Client:**
   - FastAPI `TestClient` from `fastapi.testclient`
   - Module-level: `client = TestClient(app)`
   - Used for GET/POST requests: `client.get()`, `client.post()`

## Mocking

**Framework:**
- `unittest.mock` from Python standard library
- `patch()` and `MagicMock()` for mocking external dependencies

**Patterns:**

1. **HTTP Mocking (External APIs):**
   - Location: `tests/test_pubchem_client.py`
   - Uses `@patch()` decorator to mock `requests` library
   - Creates `MagicMock()` response objects with `.json()`, `.status_code`, `.raise_for_status()`
   - Example:
   ```python
   from unittest.mock import patch, MagicMock

   with patch('requests.get') as mock_get:
       mock_response = MagicMock()
       mock_response.status_code = 200
       mock_response.json.return_value = {"PropertyTable": {...}}
       mock_get.return_value = mock_response
   ```

2. **Async Mocking:**
   - `AsyncMock` used for async functions (e.g., `test_classyfire.py`)
   - Example: `from unittest.mock import AsyncMock`

3. **Mock Client Class:**
   - Custom mock class created in tests when needed (e.g., `MockClient` in `test_pubchem_client.py`)
   - Mimics actual client interface for testing without external dependencies

**What to Mock:**
- External HTTP APIs (PubChemClient, ClassyFire)
- Slow operations (ML model inference, large computations)
- File system operations (CDX parsing, file uploads)

**What NOT to Mock:**
- Internal chemistry toolkit calls (RDKit, CDK, OpenBabel) - test integration with them
- FastAPI app and client (use real TestClient)
- Database queries (if any)
- Request/response models

## Fixtures and Factories

**Test Data:**

1. **Inline Molecules:**
   ```python
   CAFFEINE_MOL_BLOCK = """
     RDKit          2D

    14 15  0  0  0  0  0  0  0  0999 V2000
    ...
   M  END"""
   ```

2. **Test SMILES:**
   - Common test molecule: Caffeine `"CN1C=NC2=C1C(=O)N(C(=O)N2C)C"`
   - Simple molecules: `"CCO"` (ethanol), `"c1ccccc1"` (benzene)
   - Edge cases: `"C1#CC#C1"` (triple bonds), `"C[C@H](Cl)Br"` (stereo)

3. **Fixtures defined in test modules:**
   ```python
   @pytest.fixture
   def test_smiles():
       return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

   @pytest.fixture
   def invalid_smiles():
       return "invalid"
   ```

**Location:**
- Fixtures defined in test modules themselves (no shared `conftest.py`)
- Test data (images, structures) in `tests/` directory: `caffeine.png`, `segment_sample.png`, `ketcher.cdx`

**Factory Pattern:**
- Not extensively used
- Test utilities created on-the-fly within test functions
- Example: Mock response objects created with `MagicMock()` and `.return_value` assignments

## Coverage

**Requirements:** ~88% (as of March 2026, per MEMORY.md)

**Target Coverage Path:**
- Ongoing effort to push beyond 88%
- Focus on untested error paths and edge cases
- Coverage tests explicitly target uncovered lines (e.g., comments like `# --- Lines 303-311: /smiles endpoint error handling ---`)

**View Coverage:**
```bash
python3 -m pytest --cov=./ --cov-report=html    # Generate HTML report
python3 -m pytest --cov=./ --cov-report=term    # Terminal summary
python3 -m pytest --cov=./ --cov-report=xml     # XML for CI (Codecov)
```

**CI Integration:**
- Codecov token stored in GitHub Actions secrets (`CODECOV_TOKEN`)
- Coverage report uploaded automatically on push/PR
- XML coverage file generated: `coverage.xml`

## Test Types

**Unit Tests:**
- Scope: Individual functions and toolkit wrappers
- Approach: Test single toolkit (RDKit, CDK, or OpenBabel) in isolation
- Example: `test_helper.py` tests `parse_input()` with each framework separately
- Tool-specific tests: `test_cdk_*.py` files isolate CDK depiction features

**Integration Tests:**
- Scope: API endpoints with real toolkits (no mocking of chemistry)
- Approach: Use `TestClient` to make real HTTP calls to `app`
- Example: `test_converters.py` tests complete SMILES→InChI pipeline
- Multi-step workflows: Batch operations in `test_converters.py:test_batch_*`

**E2E Tests:**
- Framework: Not formally organized as "E2E"
- Approach: Full API endpoint testing with `TestClient` (closest to E2E)
- Complex workflows: Batch conversion, format transformations
- Image uploads: CDX file processing, OCSR on images (`test_ocsr.py`)

## Common Patterns

**Async Testing:**
- Test functions are synchronous (`def test_*()`), not async
- `TestClient` handles async FastAPI routes internally
- Async mocking done with `AsyncMock()` when needed

**Error Testing:**
```python
def test_invalid_canonical_smiles(invalid_smiles):
    with pytest.raises(InvalidInputException):
        get_ob_canonical_SMILES(invalid_smiles)
```

**Parametrized Error Testing:**
```python
@pytest.mark.parametrize(
    "smiles, response_code",
    [
        ("CC", 200),
        ("INVALID_INPUT", 422),
    ],
)
def test_endpoint(smiles, response_code):
    response = client.get(f"/endpoint?smiles={smiles}")
    assert response.status_code == response_code
```

**Status Code Assertions:**
- 200: Success
- 400: Bad request (validation error)
- 404: Not found
- 422: Unprocessable entity (invalid input)
- 500: Server error

**Ignore Warnings:**
```python
import warnings

@pytest.fixture(autouse=True)
def ignore_deprecation_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
```

**Rate Limiter Reset:**
```python
@pytest.fixture(autouse=True)
def _reset_limiter():
    """Reset rate limiter state before each test."""
    limiter.reset()
    yield
```

## Test Execution Details

**Python Version:** 3.10 (enforced in CI)

**Dependencies Installed:**
- `pytest` and `pytest-cov` for testing
- `flake8` for linting
- All application dependencies from `requirements.txt`
- Optional: DECIMER (deep learning OCSR model) via pip
- External: `surge` binary (structure generator)

**Pre-test Setup (from CI):**
1. Install dependencies
2. Install DECIMER deep learning model
3. Download and install `surge` binary to `/usr/bin`
4. Run flake8 linter

**Test Execution Order:**
- pytest discovers all `test_*.py` files
- Runs tests in file order by default
- Within files: tests run top-to-bottom
- Parametrized tests expand to multiple test cases

---

*Testing analysis: 2026-03-12*
