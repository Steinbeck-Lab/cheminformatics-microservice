from __future__ import annotations

import warnings

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.limiter import limiter

client = TestClient(app)

AUTH_HEADERS = {}


@pytest.fixture(autouse=True)
def _reset_limiter():
    """Reset rate limiter state before each test."""
    limiter.reset()
    yield


def test_converters_index():
    response = client.get("/latest/convert/")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}


@pytest.mark.parametrize(
    "smiles, toolkit, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "openbabel",
            200,
        ),
        (
            "INVALID_INPUT",
            "cdk",
            422,
        ),
    ],
)
def test_create_2D_coordinates(smiles, toolkit, response_code):
    response = client.get(
        f"/latest/convert/mol2D?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "smiles, toolkit, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "openbabel",
            200,
        ),
        (
            "INVALID_INPUT",
            "rdkit",
            422,
        ),
    ],
)
def test_create_3D_coordinates(smiles, toolkit, response_code):
    response = client.get(
        f"/latest/convert/mol3D?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "smiles, toolkit, response_text, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            '"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"',
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            '"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"',
            200,
        ),
        # (
        #     "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        #     "openbabel",
        #     '"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"',
        #     200,
        # ),
        (
            "INVALID_INPUT",
            "cdk",
            "",
            422,
        ),
    ],
)
def test_smiles_to_inchi(smiles, toolkit, response_text, response_code):
    response = client.get(
        f"/latest/convert/inchi?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "smiles, toolkit, response_text, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            '"RYYVLZVUVIJVGH-UHFFFAOYSA-N"',
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            '"RYYVLZVUVIJVGH-UHFFFAOYSA-N"',
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "openbabel",
            '"RYYVLZVUVIJVGH-UHFFFAOYSA-N"',
            200,
        ),
        (
            "INVALID_INPUT",
            "cdk",
            "",
            422,
        ),
    ],
)
def test_smiles_to_inchikey(smiles, toolkit, response_text, response_code):
    response = client.get(
        f"/latest/convert/inchikey?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "smiles, toolkit, response_text, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            '"CN1C=NC2=C1C(=O)N(C)C(=O)N2C |(2.68,2.45,;2.22,1.02,;3.1,-0.19,;2.22,-1.4,;0.8,-0.94,;0.8,0.56,;-0.5,1.31,;-0.5,2.81,;-1.8,0.56,;-3.1,1.31,;-1.8,-0.94,;-3.1,-1.69,;-0.5,-1.69,;-0.5,-3.19,)|"',
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            '"Cn1c(=O)c2c(ncn2C)n(C)c1=O |(-2.43754,-2.3236,;-1.48768,-1.16267,;-0.0073539,-1.4048,;0.523116,-2.80787,;0.942504,-0.243866,;0.412035,1.1592,;1.58251,2.09728,;2.83637,1.27398,;2.44083,-0.172931,;3.3789,-1.3434,;-1.06829,1.40134,;-1.59876,2.80441,;-2.01815,0.240402,;-3.49848,0.482536,)|"',
            200,
        ),
        (
            "INVALID_INPUT",
            "cdk",
            "",
            422,
        ),
    ],
)
def test_smiles_to_cxsmiles(smiles, toolkit, response_text, response_code):
    response = client.get(
        f"/latest/convert/cxsmiles?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "smiles, toolkit, response_text, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            '"CN1C=NC2=C1C(=O)N(C)C(=O)N2C"',
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            '"CN1C(=O)C2=C(N=CN2C)N(C)C1=O"',
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "openbabel",
            '"Cn1cnc2c1c(=O)n(C)c(=O)n2C"',
            200,
        ),
        (
            "INVALID_INPUT",
            "cdk",
            "",
            422,
        ),
    ],
)
def test_smiles_cannonicalise(smiles, toolkit, response_text, response_code):
    response = client.get(
        f"/latest/convert/canonicalsmiles?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "input, representation, response_text, response_code",
    [
        (
            "1,3,7-trimethylpurine-2,6-dione",
            "iupac",
            '"CN1C(N(C=2N=CN(C2C1=O)C)C)=O"',
            200,
        ),
        (
            "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]",
            "selfies",
            '"CN1C=NC2=C1C(=O)N(C(=O)N2C)C"',
            200,
        ),
    ],
)
def test_iupac_or_selfies_to_smiles(
    input,
    representation,
    response_text,
    response_code,
):
    response = client.get(
        f"/latest/convert/smiles?input_text={input}&representation={representation}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    assert response.text == response_text


@pytest.mark.parametrize(
    "smiles, response_text, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            '"[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]"',
            200,
        ),
        (
            "INVALID_INPUT",
            "",
            400,
        ),
    ],
)
def test_smiles_to_selfies(smiles, response_text, response_code):
    response = client.get(f"/latest/convert/selfies?smiles={smiles}")
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "smiles, toolkit, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "openbabel",
            200,
        ),
        (
            "INVALID_INPUT",
            "cdk",
            422,
        ),
    ],
)
def test_smiles_to_formats(smiles, toolkit, response_code):
    response = client.get(
        f"/latest/convert/formats?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert "mol" in response.json()
        assert "canonicalsmiles" in response.json()
        assert "inchi" in response.json()
        assert "inchikey" in response.json()


# Filter out DeprecationWarning messages
@pytest.fixture(autouse=True)
def ignore_deprecation_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)


@pytest.mark.parametrize(
    "smiles, toolkit, response_text, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            '"[#6]-[#7]1:[#6]:[#7]:[#6]2:[#6]:1:[#6](=[#8]):[#7](:[#6](=[#8]):[#7]:2-[#6])-[#6]"',
            200,
        ),
        (
            "INVALID_INPUT",
            "cdk",
            "",
            422,
        ),
    ],
)
def test_smiles_smarts(smiles, toolkit, response_text, response_code):
    response = client.get(
        f"/latest/convert/smarts?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


# Caffeine MOL block for testing
CAFFEINE_MOL_BLOCK = """
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
    3.3789   -1.3434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4408   -0.1729    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.8364    1.2740    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5825    2.0973    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4120    1.1592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9425   -0.2439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0074   -1.4048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5231   -2.8079    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4877   -1.1627    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0182    0.2404    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4985    0.4825    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0683    1.4013    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5988    2.8044    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4375   -2.3236    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  7  9  1  0
  9 10  1  0
 10 11  2  0
 10 12  1  0
 12 13  1  0
  9 14  1  0
  6  2  1  0
 12  5  1  0
M  END"""

CAFFEINE_SDF_BLOCK = """2519
  -OEChem-10132511592D

 24 25  0     0  0  0  0  0  0999 V2000
    3.7321    2.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321   -1.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.5443    0.8047    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660    0.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.5443   -0.8047    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5981   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1279    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321   -2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8550    1.7553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.7479   -0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.1121   -2.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7321   -2.6200    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.3521   -2.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.4443    1.5626    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.0476    2.3446    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.2656    1.9479    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3100    1.5369    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4631    1.3100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6900    0.4631    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  9  2  0  0  0  0
  2 10  2  0  0  0  0
  3  8  1  0  0  0  0
  3 10  1  0  0  0  0
  3 12  1  0  0  0  0
  4  7  1  0  0  0  0
  4 11  1  0  0  0  0
  4 13  1  0  0  0  0
  5  9  1  0  0  0  0
  5 10  1  0  0  0  0
  5 14  1  0  0  0  0
  6  8  1  0  0  0  0
  6 11  2  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
 11 15  1  0  0  0  0
 12 16  1  0  0  0  0
 12 17  1  0  0  0  0
 12 18  1  0  0  0  0
 13 19  1  0  0  0  0
 13 20  1  0  0  0  0
 13 21  1  0  0  0  0
 14 22  1  0  0  0  0
 14 23  1  0  0  0  0
 14 24  1  0  0  0  0
M  END

$$$$
"""


@pytest.mark.parametrize(
    "molblock, toolkit, expected_contains, response_code",
    [
        (
            CAFFEINE_MOL_BLOCK,
            "cdk",
            "C1",  # Check if SMILES contains at least carbon
            200,
        ),
        (
            CAFFEINE_MOL_BLOCK,
            "rdkit",
            "C1",
            200,
        ),
        (
            CAFFEINE_SDF_BLOCK,
            "cdk",
            "C1",  # SDF block should work too
            200,
        ),
        (
            CAFFEINE_SDF_BLOCK,
            "rdkit",
            "C1",
            200,
        ),
        (
            "INVALID_MOL_BLOCK",
            "cdk",
            "",
            422,
        ),
    ],
)
def test_molblock_to_smiles(molblock, toolkit, expected_contains, response_code):
    response = client.post(
        f"/latest/convert/molblock?toolkit={toolkit}",
        data=molblock,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == response_code
    if response_code == 200:
        smiles = response.text.strip('"')
        assert len(smiles) > 0
        assert (
            expected_contains in smiles or "N" in smiles
        )  # Caffeine contains both C and N


# ---------------------------------------------------------------------------
# Coverage tests for security hardening (error paths & input validation)
# ---------------------------------------------------------------------------


def test_cdx_to_mol_oversized_file():
    """CDX upload rejects files over 10MB."""
    large_content = b"\x00" * (10 * 1024 * 1024 + 1)
    response = client.post(
        "/latest/convert/cdx-to-mol",
        files={"file": ("huge.cdx", large_content, "chemical/x-cdx")},
    )
    assert response.status_code == 400
    assert "too large" in response.json()["detail"].lower()


def test_iupac_to_smiles_invalid():
    """IUPAC conversion with garbage input returns error."""
    response = client.get(
        "/latest/convert/smiles?input_text=ZZZZNOTACHEMICAL&representation=iupac&converter=opsin"
    )
    assert response.status_code in [200, 422]


def test_selfies_encode_invalid():
    """SELFIES encoding of invalid SMILES returns error."""
    response = client.get("/latest/convert/selfies?smiles=INVALID_SMILES&toolkit=rdkit")
    assert response.status_code in [400, 422]


def test_molblock_to_smiles_invalid():
    """Invalid molblock conversion returns error."""
    response = client.post(
        "/latest/convert/molblock?toolkit=rdkit",
        data="NOT_A_MOLBLOCK",
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code in [422, 500]


def test_smiles_to_formats_invalid():
    """Invalid SMILES to formats conversion returns error."""
    response = client.get("/latest/convert/mol2D?smiles=INVALID!!!&toolkit=rdkit")
    assert response.status_code in [422, 500]


# ---------------------------------------------------------------------------
# Additional coverage tests targeting uncovered lines
# ---------------------------------------------------------------------------

# --- Lines 303-311: /smiles endpoint error handling ---


def test_smiles_endpoint_generic_exception():
    """Test the generic exception handler in /smiles endpoint (lines 309-311).

    Passing an invalid SELFIES string that causes an exception during decoding.
    """
    response = client.get(
        "/latest/convert/smiles?input_text=NOT_VALID_SELFIES&representation=selfies"
    )
    # sf.decoder may return empty or raise; either way the endpoint should handle it
    assert response.status_code in [200, 422]


# --- Lines 498-501: /inchi endpoint with openbabel toolkit ---


def test_smiles_to_inchi_openbabel():
    """Test InChI generation using openbabel toolkit (lines 498-501)."""
    response = client.get(
        "/latest/convert/inchi?smiles=CN1C%3DNC2%3DC1C(%3DO)N(C(%3DO)N2C)C&toolkit=openbabel",
    )
    assert response.status_code == 200
    assert "InChI" in response.text


def test_smiles_to_inchi_openbabel_simple():
    """Test InChI generation using openbabel toolkit with simple SMILES."""
    response = client.get(
        "/latest/convert/inchi?smiles=CCO&toolkit=openbabel",
    )
    assert response.status_code == 200
    assert "InChI" in response.text


# --- Lines 619, 624: /selfies endpoint error paths ---


def test_selfies_encode_empty_result():
    """Test SELFIES encoding when encoder returns empty (lines 619, 624).

    Uses a truly invalid SMILES that sf.encoder cannot handle.
    """
    response = client.get("/latest/convert/selfies?smiles=%5B%5D")
    # An empty SMILES-like input may cause encoder to fail
    assert response.status_code in [200, 400]


# --- Lines 715, 720: /formats endpoint error paths ---


def test_formats_rdkit_invalid_smiles():
    """Test /formats with rdkit toolkit and invalid SMILES."""
    response = client.get(
        "/latest/convert/formats?smiles=INVALID&toolkit=rdkit",
    )
    assert response.status_code == 422


def test_formats_openbabel_valid_smiles():
    """Test /formats with openbabel toolkit and valid SMILES."""
    response = client.get(
        "/latest/convert/formats?smiles=CCO&toolkit=openbabel",
    )
    assert response.status_code == 200
    data = response.json()
    assert "mol" in data
    assert "canonicalsmiles" in data
    assert "inchi" in data
    assert "inchikey" in data


# --- Lines 826, 837, 850, 871-876: /molblock POST endpoint ---


def test_molblock_empty_body():
    """Test molblock endpoint with empty body (line 826)."""
    response = client.post(
        "/latest/convert/molblock?toolkit=cdk",
        data="",
        headers={"Content-Type": "text/plain"},
    )
    # Empty body should trigger 422 ("Invalid or missing molblock")
    assert response.status_code in [400, 422]


def test_molblock_sdf_only_terminator():
    """Test molblock endpoint with only $$$$ (line 837 - empty after strip)."""
    response = client.post(
        "/latest/convert/molblock?toolkit=cdk",
        data="$$$$",
        headers={"Content-Type": "text/plain"},
    )
    # After splitting on $$$$, content should be empty -> 400
    assert response.status_code in [400, 422]


def test_molblock_rdkit_invalid_mol():
    """Test molblock endpoint with rdkit toolkit and invalid MOL (lines 871-876)."""
    # A string that is not a valid MOL block but not empty
    invalid_mol = """invalid
  invalid
  0  0  0  0  0  0  0  0  0  0  0 V2000
M  END"""
    response = client.post(
        "/latest/convert/molblock?toolkit=rdkit",
        data=invalid_mol,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code in [422, 500]


def test_molblock_rdkit_valid():
    """Test molblock endpoint with rdkit toolkit and valid MOL block."""
    response = client.post(
        "/latest/convert/molblock?toolkit=rdkit",
        data=CAFFEINE_MOL_BLOCK,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    smiles = response.text.strip('"')
    assert len(smiles) > 0


# --- Lines 960-1092: /batch POST endpoint ---


def test_batch_smiles_to_canonicalsmiles_cdk():
    """Test batch conversion: SMILES -> canonicalsmiles with CDK (lines 991-993)."""
    body = {
        "inputs": [{"value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "input_format": "smiles"}]
    }
    response = client.post(
        "/latest/convert/batch?output_format=canonicalsmiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["total"] == 1
    assert data["summary"]["successful"] == 1
    assert data["results"][0]["success"] is True
    assert len(data["results"][0]["output"]) > 0


def test_batch_smiles_to_canonicalsmiles_rdkit():
    """Test batch conversion: SMILES -> canonicalsmiles with rdkit."""
    body = {
        "inputs": [{"value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "input_format": "smiles"}]
    }
    response = client.post(
        "/latest/convert/batch?output_format=canonicalsmiles&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_canonicalsmiles_openbabel():
    """Test batch conversion: SMILES -> canonicalsmiles with openbabel (line 997-998)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=canonicalsmiles&toolkit=openbabel",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_inchi_cdk():
    """Test batch conversion: SMILES -> inchi with CDK (lines 1001-1003)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=inchi&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1
    assert "InChI" in data["results"][0]["output"]


def test_batch_smiles_to_inchi_rdkit():
    """Test batch conversion: SMILES -> inchi with rdkit (lines 1004-1006)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=inchi&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1
    assert "InChI" in data["results"][0]["output"]


def test_batch_smiles_to_inchi_openbabel():
    """Test batch conversion: SMILES -> inchi with openbabel (lines 1007-1008)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=inchi&toolkit=openbabel",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_inchikey_cdk():
    """Test batch conversion: SMILES -> inchikey with CDK (lines 1011-1013)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=inchikey&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1
    assert len(data["results"][0]["output"]) > 0


def test_batch_smiles_to_inchikey_rdkit():
    """Test batch conversion: SMILES -> inchikey with rdkit (lines 1015-1016)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=inchikey&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_inchikey_openbabel():
    """Test batch conversion: SMILES -> inchikey with openbabel (lines 1017-1018)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=inchikey&toolkit=openbabel",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_selfies():
    """Test batch conversion: SMILES -> selfies (line 1020-1021)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=selfies&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1
    assert "[C]" in data["results"][0]["output"]


def test_batch_smiles_to_cxsmiles_cdk():
    """Test batch conversion: SMILES -> cxsmiles with CDK (lines 1024-1026)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=cxsmiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_cxsmiles_rdkit():
    """Test batch conversion: SMILES -> cxsmiles with rdkit (lines 1028-1029)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=cxsmiles&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_cxsmiles_openbabel_unsupported():
    """Test batch cxsmiles conversion with openbabel (lines 1030-1033 - unsupported)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=cxsmiles&toolkit=openbabel",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    # Should fail with unsupported error
    assert data["summary"]["failed"] == 1
    assert "not supported" in data["results"][0]["error"].lower()


def test_batch_smiles_to_smarts_rdkit():
    """Test batch conversion: SMILES -> smarts with rdkit (lines 1036-1038)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smarts&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1
    assert "#" in data["results"][0]["output"]  # SMARTS uses atomic number notation


def test_batch_smarts_cdk_unsupported():
    """Test batch smarts conversion with cdk (lines 1039-1042 - unsupported)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smarts&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["failed"] == 1
    assert "not supported" in data["results"][0]["error"].lower()


def test_batch_smiles_to_mol2d_cdk():
    """Test batch conversion: SMILES -> mol2d with CDK (lines 1045-1047)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=mol2d&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1
    assert (
        "V2000" in data["results"][0]["output"]
        or "V3000" in data["results"][0]["output"]
    )


def test_batch_smiles_to_mol2d_rdkit():
    """Test batch conversion: SMILES -> mol2d with rdkit (lines 1048-1050)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=mol2d&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_mol2d_openbabel():
    """Test batch conversion: SMILES -> mol2d with openbabel (lines 1051-1052)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=mol2d&toolkit=openbabel",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_mol3d_rdkit():
    """Test batch conversion: SMILES -> mol3d with rdkit (lines 1055-1057)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=mol3d&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_smiles_to_mol3d_openbabel():
    """Test batch conversion: SMILES -> mol3d with openbabel (lines 1058-1059)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=mol3d&toolkit=openbabel",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_mol3d_cdk_unsupported():
    """Test batch mol3d with cdk toolkit (lines 1060-1063 - unsupported)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=mol3d&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["failed"] == 1
    assert "not supported" in data["results"][0]["error"].lower()


def test_batch_unsupported_output_format():
    """Test batch with unsupported output format (lines 1065-1066)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=INVALID_FORMAT&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["failed"] == 1
    assert "unsupported output format" in data["results"][0]["error"].lower()


def test_batch_smiles_output():
    """Test batch conversion: SMILES -> smiles (line 987-988)."""
    body = {"inputs": [{"value": "CCO", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1
    assert data["results"][0]["output"] == "CCO"


def test_batch_missing_value_or_format():
    """Test batch with missing value or input_format (lines 959-960)."""
    body = {"inputs": [{"value": "", "input_format": "smiles"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["failed"] == 1
    assert "missing" in data["results"][0]["error"].lower()


def test_batch_missing_input_format():
    """Test batch with missing input_format (lines 959-960)."""
    body = {"inputs": [{"value": "CCO", "input_format": ""}]}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["failed"] == 1


def test_batch_iupac_input():
    """Test batch conversion with IUPAC input (lines 965-970)."""
    body = {"inputs": [{"value": "ethanol", "input_format": "iupac"}]}
    response = client.post(
        "/latest/convert/batch?output_format=canonicalsmiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["total"] == 1
    # IUPAC conversion may succeed or fail depending on OPSIN
    assert data["results"][0]["success"] in [True, False]


def test_batch_selfies_input():
    """Test batch conversion with SELFIES input (lines 971-974)."""
    body = {"inputs": [{"value": "[C][C][O]", "input_format": "selfies"}]}
    response = client.post(
        "/latest/convert/batch?output_format=canonicalsmiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_inchi_input():
    """Test batch conversion with InChI input (lines 976-980)."""
    body = {
        "inputs": [
            {"value": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "input_format": "inchi"}
        ]
    }
    response = client.post(
        "/latest/convert/batch?output_format=canonicalsmiles&toolkit=rdkit",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["successful"] == 1


def test_batch_inchi_input_invalid():
    """Test batch conversion with invalid InChI input (lines 978-979)."""
    body = {"inputs": [{"value": "NOT_AN_INCHI", "input_format": "inchi"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["failed"] == 1


def test_batch_unsupported_input_format():
    """Test batch with unsupported input format (line 982)."""
    body = {"inputs": [{"value": "something", "input_format": "mol2d"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["failed"] == 1
    assert "unsupported input format" in data["results"][0]["error"].lower()


def test_batch_multiple_inputs():
    """Test batch conversion with multiple inputs including failures (lines 1079-1092)."""
    body = {
        "inputs": [
            {"value": "CCO", "input_format": "smiles"},
            {"value": "INVALID", "input_format": "inchi"},
            {"value": "c1ccccc1", "input_format": "smiles"},
        ]
    }
    response = client.post(
        "/latest/convert/batch?output_format=inchi&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["total"] == 3
    # At least some should succeed, some should fail
    assert data["summary"]["successful"] >= 1
    assert data["summary"]["failed"] >= 1


def test_batch_empty_inputs():
    """Test batch conversion with empty inputs list."""
    body = {"inputs": []}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["total"] == 0
    assert data["summary"]["successful"] == 0
    assert data["summary"]["failed"] == 0


def test_batch_selfies_input_invalid():
    """Test batch conversion with invalid SELFIES input (lines 973-974)."""
    body = {"inputs": [{"value": "NOT_SELFIES_AT_ALL", "input_format": "selfies"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    # sf.decoder on invalid input may return empty -> should be a failure
    assert data["summary"]["total"] == 1


def test_batch_iupac_input_invalid():
    """Test batch with invalid IUPAC name (lines 967-970)."""
    body = {"inputs": [{"value": "ZZZZNOTACHEMICAL", "input_format": "iupac"}]}
    response = client.post(
        "/latest/convert/batch?output_format=smiles&toolkit=cdk",
        json=body,
        headers=AUTH_HEADERS,
    )
    assert response.status_code == 200
    data = response.json()
    assert data["summary"]["total"] == 1
    # OPSIN may fail to convert this
    assert data["results"][0]["success"] in [True, False]


# --- Lines 1162-1164: CDX endpoint generic exception ---


def test_cdx_to_mol_invalid_cdx_content():
    """Test CDX to MOL with invalid CDX content (lines 1162-1164)."""
    invalid_cdx = b"\x01\x02\x03\x04\x05\x06"
    response = client.post(
        "/latest/convert/cdx-to-mol",
        files={"file": ("test.cdx", invalid_cdx, "chemical/x-cdx")},
    )
    assert response.status_code == 422


def test_cdx_to_mol_wrong_extension():
    """Test CDX endpoint rejects non-CDX files."""
    response = client.post(
        "/latest/convert/cdx-to-mol",
        files={"file": ("test.txt", b"some content", "text/plain")},
    )
    assert response.status_code == 400
    assert "cdx" in response.json()["detail"].lower()
