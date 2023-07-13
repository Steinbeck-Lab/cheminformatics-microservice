import pytest
from fastapi.testclient import TestClient
from app.main import app


client = TestClient(app)


def test_chem_index():
    response = client.get("/v1/depict/")
    assert response.status_code == 200
    assert response.json() == {
        "module": "depict",
        "message": "Successful",
        "status": 200,
    }


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.mark.parametrize(
    "smiles, generator, width, height, rotate, CIP, unicolor",
    [
        ("CCO", "cdk", 512, 512, 0, False, False),
        ("CCO", "rdkit", 512, 512, 0, False, False),
    ],
)
def test_depict2D_molecule(smiles, generator, width, height, rotate, CIP, unicolor):
    response = client.get(
        f"/v1/depict/2D?smiles={smiles}&generator={generator}&width={width}&height={height}&rotate={rotate}&CIP={CIP}&unicolor={unicolor}"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "image/svg+xml"


def test_depict3D_molecule(test_smiles):
    response = client.get(f"/v1/depict/3D?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/html; charset=utf-8"


# Run the tests
if __name__ == "__main__":
    pytest.main()
