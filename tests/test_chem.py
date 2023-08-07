import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


def test_chem_index():
    response = client.get("/latest/chem/")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}


def test_smiles_to_stereo_isomers(test_smiles):
    response = client.get(f"/latest/chem/stereoisomers?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '["Cn1c(=O)c2c(ncn2C)n(C)c1=O"]'


def test_SMILES_Descriptors_returns_descriptors(test_smiles):
    response = client.get(f"/latest/chem/descriptors?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert (
        response.text
        == '{"atom_count":24,"heavy_atom_count":14,"molecular_weight":194.19,"exactmolecular_weight":194.08038,"alogp":-1.03,"rotatable_bond_count":0,"topological_polar_surface_area":61.82,"hydrogen_bond_acceptors":6,"hydrogen_bond_donors":0,"hydrogen_bond_acceptors_lipinski":6,"hydrogen_bond_donors_lipinski":0,"lipinski_rule_of_five_violations":0,"aromatic_rings_count":2,"qed_drug_likeliness":0.54,"formal_charge":0,"fractioncsp3":0.375,"number_of_minimal_rings":2,"van_der_walls_volume":"None","linear_sugars":false,"circular_sugars":false,"murko_framework":"N1=C[N]C2=C1NCNC2","nplikeness":-1.09}'
    )

    response = client.get(f"/latest/chem/descriptors?smiles={test_smiles}&format=html")
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/html; charset=utf-8"


@pytest.mark.parametrize(
    "smiles, expected_score",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", -1.09),
        ("CC=O", 0.96),
    ],
)
def test_NPlikeliness_Score(smiles, expected_score):
    response = client.get(f"/latest/chem/nplikeness/score?smiles={smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.json() == expected_score


"""
def test_ClassyFire_Classify(test_smiles):
    response = client.get(f"/latest/chem/classyfire/classify?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
"""


@pytest.mark.parametrize(
    "smiles, toolkit, expected",
    [
        ("CC,CC", "cdk", '"1.00000"'),
        ("CC,CC", "rdkit", "1.0"),
        (
            "CC,CC,CC",
            "cdk",
            "<table><tr><th></th><th>0</th><th>1</th><th>2</th></tr><tr><td>0</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr><tr><td>1</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr><tr><td>2</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr></table>",
        ),
    ],
)
def test_tanimoto_similarity(smiles, toolkit, expected):
    response = client.get(f"/latest/chem/tanimoto?smiles={smiles}&toolkit={toolkit}")
    assert response.status_code == 200
    assert response.text == expected


@pytest.mark.parametrize(
    "smiles, fix, expected",
    [
        ("CCO", False, '{"smi":"CCO","messages":["No Errors Found"]}'),
        ("CCO", True, '{"smi":"CCO","messages":["No Errors Found"]}'),
    ],
)
def test_check_errors(smiles, fix, expected):
    response = client.get(f"/latest/chem/errors?smiles={smiles}&fix={fix}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == expected


def test_hose_codes(test_smiles):
    response = client.get(
        f"/latest/chem/HOSEcode?smiles={test_smiles}&spheres=0&toolkit=cdk&ringsize=false"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert (
        response.text
        == '["C-4;(//)","N-3;(//)","C-3;(//)","N-2;(//)","C-3;(//)","C-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-4;(//)","C-4;(//)"]'
    )
    # Add more assertions for the expected response data


def test_coconut_preprocessing(test_smiles):
    response = client.get(f"/latest/chem/coconut/pre-processing?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


# Run the tests
if __name__ == "__main__":
    pytest.main()
