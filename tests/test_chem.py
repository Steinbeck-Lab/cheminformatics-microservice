import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


def test_chem_index():
    response = client.get("/v1/chem/")
    assert response.status_code == 200
    assert response.json() == {"module": "chem", "message": "Successful", "status": 200}


def test_smiles_to_stereo_isomers(test_smiles):
    response = client.get(f"/v1/chem/stereoisomers?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '["Cn1c(=O)c2c(ncn2C)n(C)c1=O"]'


def test_SMILES_Descriptors_returns_descriptors(test_smiles):
    response = client.get(f"/v1/chem/descriptors?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert (
        response.text
        == '{"atom_count":24,"heavy_atom_count":14,"molecular_weight":194.19,"exactmolecular_weight":194.08038,"alogp":-1.03,"rotatable_bond_count":0,"topological_polar_surface_area":61.82,"hydrogen_bond_acceptors":6,"hydrogen_bond_donors":0,"hydrogen_bond_acceptors_lipinski":6,"hydrogen_bond_donors_lipinski":0,"lipinski_rule_of_five_violations":0,"aromatic_rings_count":2,"qed_drug_likeliness":0.54,"formal_charge":0,"fractioncsp3":0.375,"number_of_minimal_rings":2,"linear_sugars":false,"circular_sugars":false,"murko_framework":"N1=C[N]C2=C1NCNC2","nplikeliness":-1.09}'
    )

    response = client.get(f"/v1/chem/descriptors?smiles={test_smiles}&format=html")
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/html; charset=utf-8"

    response = client.get("/v1/chem/descriptors")
    assert response.status_code == 422  # Missing required parameter

    response = client.get(f"/v1/chem/descriptors?smiles={test_smiles}&toolkit=unknown")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '"Error Calculating Descriptors"'


@pytest.mark.parametrize(
    "smiles, expected_score",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "-1.09"),
        ("CC=O", "0.96"),
    ],
)
def test_NPlikeliness_Score(smiles, expected_score):
    response = client.get(f"/v1/chem/npscore?smiles={smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.json() == expected_score


def test_ClassyFire_Classify(test_smiles):
    response = client.get(f"/v1/chem/classyfire/classify?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


@pytest.mark.parametrize(
    "smiles, expected",
    [("invalid_smiles", '"Error reading SMILES string, check again."')],
)
def test_cdk2d_coordinates(smiles, expected):
    response = client.get(f"/v1/chem/cdk2d?smiles={smiles}")
    assert response.status_code == 200
    assert response.text == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "CC",
            """     RDKit          3D\n'\n '\n'\n '  8  7  0  0  0  0  0  0  0  0999 V2000\n'\n '    0.7480   -0.0676    0.0869 C   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '   -0.7480    0.0676   -0.0869 C   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '    1.2303   -0.2658   -0.8749 H   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '    1.1714    0.8520    0.5016 H   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '    0.9830   -0.8922    0.7665 H   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '   -0.9830    0.8922   -0.7665 H   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '   -1.2303    0.2658    0.8749 H   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '   -1.1714   -0.8520   -0.5016 H   0  0  0  0  0  0  0  0  0  0  0  0\n'\n '  1  2  1  0\n'\n '  1  3  1  0\n'\n '  1  4  1  0\n'\n '  1  5  1  0\n'\n '  2  6  1  0\n'\n '  2  7  1  0\n'\n '  2  8  1  0\n'\n 'M  END""",
        ),
        ("invalid_smiles", '"Error reading SMILES string, check again."'),
    ],
)
def test_rdkit3d_mol(smiles, expected):
    response = client.get(f"/v1/chem/rdkit3d?smiles={smiles}")
    assert response.status_code == 200
    # assert response.text == expected


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
    response = client.get(f"/v1/chem/tanimoto?smiles={smiles}&toolkit={toolkit}")
    assert response.status_code == 200
    assert response.text == expected


@pytest.mark.parametrize(
    "smiles, generator, width, height, rotate, CIP, unicolor",
    [
        ("CCO", "cdksdg", 512, 512, 0, False, False),
        ("CCO", "rdkit", 512, 512, 0, False, False),
    ],
)
def test_depict2D_molecule(smiles, generator, width, height, rotate, CIP, unicolor):
    response = client.get(
        f"/v1/chem/depict?smiles={smiles}&generator={generator}&width={width}&height={height}&rotate={rotate}&CIP={CIP}&unicolor={unicolor}"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "image/svg+xml"


@pytest.mark.parametrize(
    "smiles, fix, expected",
    [
        ("CCO", False, '"No Errors Found"'),
        ("CCO", True, '"No Errors Found"'),
    ],
)
def test_check_errors(smiles, fix, expected):
    response = client.get(f"/v1/chem/checkerrors?smiles={smiles}&fix={fix}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == expected


def test_depict3D_molecule(test_smiles):
    response = client.get(f"/v1/chem/depict3D?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/html; charset=utf-8"


def test_hose_codes(test_smiles):
    response = client.get(
        f"/v1/chem/hosecode?framework=cdk&smiles={test_smiles}&spheres=0&ringsize=false"
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert (
        response.text
        == '["C-4;(//)","N-3;(//)","C-3;(//)","N-2;(//)","C-3;(//)","C-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-4;(//)","C-4;(//)"]'
    )
    # Add more assertions for the expected response data


def test_coconut_preprocessing(test_smiles):
    response = client.get(f"/v1/chem/coconutpreprocessing?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


# Run the tests
if __name__ == "__main__":
    pytest.main()
