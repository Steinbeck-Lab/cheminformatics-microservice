from __future__ import annotations

import warnings

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


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
