import pytest
from fastapi.testclient import TestClient
from app.main import app


client = TestClient(app)


def test_chem_index():
    response = client.get("/latest/depict/")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}


@pytest.mark.parametrize(
    "smiles, generator, width, height, rotate, CIP, unicolor,response_code",
    [
        ("CCO", "cdk", 512, 512, 0, False, False, 200),
        ("CCO", "rdkit", 512, 512, 0, False, False, 200),
        ("INVALID_INPUT", "cdk", 512, 512, 0, False, False, 422),
    ],
)
def test_depict2D_molecule(
    smiles, generator, width, height, rotate, CIP, unicolor, response_code
):
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={generator}&width={width}&height={height}&rotate={rotate}&CIP={CIP}&unicolor={unicolor}"
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"


@pytest.mark.parametrize(
    "smiles, toolkit,response_code",
    [
        ("CCO", "openbabel", 200),
        ("CCO", "rdkit", 200),
        ("INVALID_INPUT", "openbabel", 422),
    ],
)
def test_depict3D_molecule(smiles, toolkit, response_code):
    response = client.get(f"/latest/depict/3D?smiles={smiles}&toolkit={toolkit}")
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json; charset=utf-8"


# Run the tests
if __name__ == "__main__":
    pytest.main()
