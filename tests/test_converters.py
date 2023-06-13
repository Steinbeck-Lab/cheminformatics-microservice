import pytest
import warnings
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)


def test_converters_index():
    response = client.get("/v1/convert/")
    assert response.status_code == 200
    assert response.json() == {
        "module": "converters",
        "message": "Successful",
        "status": 200,
    }


def test_SMILES_Mol():
    response = client.get(
        "/v1/convert/mol2D?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&generator=cdk"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/plain; charset=utf-8"


def test_SMILES_Generate3DConformer():
    response = client.get(
        "/v1/convert/mol3D?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&generator=rdkit"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/plain; charset=utf-8"


def test_SMILES_Canonicalise():
    response = client.get(
        "/v1/convert/canonicalsmiles?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '"CN1C=NC2=C1C(=O)N(C)C(=O)N2C"'


def test_SMILES_to_InChI():
    response = client.get("/v1/convert/inchi?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert (
        response.text
        == '"InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"'
    )


def test_SMILES_to_InChIKey():
    response = client.get("/v1/convert/inchikey?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '"RYYVLZVUVIJVGH-UHFFFAOYSA-N"'


def test_SMILES_to_CXSMILES():
    response = client.get("/v1/convert/cxsmiles?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert (
        response.text
        == '"CN1C=NC2=C1C(=O)N(C)C(=O)N2C |(2.68,2.45,;2.22,1.02,;3.1,-0.19,;2.22,-1.4,;0.8,-0.94,;0.8,0.56,;-0.5,1.31,;-0.5,2.81,;-1.8,0.56,;-3.1,1.31,;-1.8,-0.94,;-3.1,-1.69,;-0.5,-1.69,;-0.5,-3.19,)|"'
    )


def test_SMILES_convert_to_Formats():
    response = client.get(
        "/v1/convert/formats?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&generator=rdkit"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.json() == {
        "mol": "\n     RDKit          2D\n\n 14 15  0  0  0  0  0  0  0  0999 V2000\n    2.7760    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2760    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3943    1.2135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0323    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0323   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3943   -1.2135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7062   -2.6807    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1328   -3.1443    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4086   -3.6844    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.8351   -3.2209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9499   -4.2246    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1470   -1.7537    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.5736   -1.2902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0967   -5.1517    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  2  0\n  6  7  1  0\n  7  8  2  0\n  7  9  1  0\n  9 10  1  0\n 10 11  2  0\n 10 12  1  0\n 12 13  1  0\n  9 14  1  0\n  6  2  1  0\n 12  5  1  0\nM  END\n",
        "canonicalsmiles": "CN1C(=O)C2=C(N=CN2C)N(C)C1=O",
        "inchi": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
        "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
    }


def test_SMILES_to_IUPACname():
    response = client.get("/v1/convert/iupac?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '"1,3,7-trimethylpurine-2,6-dione"'


def test_IUPACname_or_SELFIES_to_SMILES():
    response = client.get(
        "/v1/convert/smiles?input_text=1,3,7-trimethylpurine-2,6-dione&representation=iupac"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '"CN1C=NC2=C1C(=O)N(C)C(=O)N2C"'

    response = client.get(
        "/v1/convert/smiles?input_text=[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]&representation=selfies"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '"CN1C=NC2=C1C(=O)N(C(=O)N2C)C"'


def test_encode_SELFIES():
    response = client.get("/v1/convert/selfies?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert (
        response.text
        == '"[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]"'
    )


# Filter out DeprecationWarning messages
@pytest.fixture(autouse=True)
def ignore_deprecation_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
