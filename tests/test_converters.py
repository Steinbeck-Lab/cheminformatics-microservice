import pytest
import warnings
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
    response = client.get(f"/latest/convert/mol2D?smiles={smiles}&toolkit={toolkit}")
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
            "cdk",
            422,
        ),
    ],
)
def test_create_3D_coordinates(smiles, toolkit, response_code):
    response = client.get(f"/latest/convert/mol3D?smiles={smiles}&toolkit={toolkit}")
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
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "openbabel",
            '"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"',
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
def test_smiles_to_inchi(smiles, toolkit, response_text, response_code):
    response = client.get(f"/latest/convert/inchi?smiles={smiles}&toolkit={toolkit}")
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
    response = client.get(f"/latest/convert/inchikey?smiles={smiles}&toolkit={toolkit}")
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
    response = client.get(f"/latest/convert/cxsmiles?smiles={smiles}&toolkit={toolkit}")
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
        f"/latest/convert/canonicalsmiles?smiles={smiles}&toolkit={toolkit}"
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "smiles, response_text, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            '"1,3,7-trimethylpurine-2,6-dione"',
            200,
        ),
        (
            "INVALID_INPUT",
            "",
            422,
        ),
    ],
)
def test_smiles_to_iupac(smiles, response_text, response_code):
    response = client.get(f"/latest/convert/iupac?smiles={smiles}")
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
            '"CN1C=NC2=C1C(=O)N(C)C(=O)N2C"',
            200,
        ),
        (
            "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]",
            "selfies",
            '"CN1C=NC2=C1C(=O)N(C)C(=O)N2C"',
            200,
        ),
    ],
)
def test_iupac_or_selfies_to_smiles(
    input, representation, response_text, response_code
):
    response = client.get(
        f"/latest/convert/smiles?input_text={input}&representation={representation}"
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
    response = client.get(f"/latest/convert/formats?smiles={smiles}&toolkit={toolkit}")
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
