import pytest
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def molfile():
    return """
  CDK     08302311362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.9743    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    1.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    2.8125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6238    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  END
"""


def test_chem_index():
    response = client.get("/latest/chem/")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}


def test_successful_smiles_to_stereo_isomers(test_smiles):
    response = client.get(f"/latest/chem/stereoisomers?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == '["Cn1c(=O)c2c(ncn2C)n(C)c1=O"]'


def test_exception_smiles_to_stereo_isomers():
    invalid_smiles = "INVALID_SMILES"
    response = client.get(f"/latest/chem/stereoisomers?smiles={invalid_smiles}")
    assert response.status_code == 400
    assert response.headers["content-type"] == "application/json"


def test_smiles_descriptors(test_smiles):
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
    "multiple_smiles, toolkit, success_response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C, CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            "cdk",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C, CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            "rdkit",
            200,
        ),
    ],
)
def test_successful_smiles_descriptors_multiple(
    multiple_smiles, toolkit, success_response_code
):
    response = client.get(
        f"/latest/chem/descriptors/multiple?smiles={multiple_smiles}&toolkit={toolkit}"
    )
    assert response.status_code == success_response_code
    assert response.headers["content-type"] == "application/json"


@pytest.mark.parametrize(
    "invalid_input, toolkit, exception_response_code",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "cdk", 400),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "rdkit", 400),
    ],
)
def test_exception_smiles_descriptors_multiple(
    invalid_input, toolkit, exception_response_code
):
    response = client.get(
        f"/latest/chem/descriptors/multiple?smiles={invalid_input}&toolkit={toolkit}"
    )
    assert response.status_code == exception_response_code


@pytest.mark.parametrize(
    "smiles, expected_score", [("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", -1.09), ("CC=O", 0.96)]
)
def test_successful_NPlikeliness_Score(smiles, expected_score):
    response = client.get(f"/latest/chem/nplikeness/score?smiles={smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.json() == expected_score


@pytest.mark.parametrize(
    "invalid_smiles, exception_response_code", [("INVALID_SMILES", 400)]
)
def test_exception_NPlikeliness_Score(invalid_smiles, exception_response_code):
    response = client.get(f"/latest/chem/nplikeness/score?smiles={invalid_smiles}")
    assert response.status_code == exception_response_code


def test_successful_classyFire_classify(test_smiles):
    response = client.get(f"/latest/chem/classyfire/classify?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


def test_successful_classyFire_result():
    job_id = 11212508
    response = client.get(f"/latest/chem/classyfire/{job_id}/result")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


# def test_exception_classyFire_classify(test_smiles):
#     response = client.get(f"/latest/chem/classyfire/classify?smiles={test_smiles}")
#     assert response.status_code == 200
#     assert response.headers["content-type"] == "application/json"


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
def test_successful_tanimoto_similarity(smiles, toolkit, expected):
    response = client.get(f"/latest/chem/tanimoto?smiles={smiles}&toolkit={toolkit}")
    assert response.status_code == 200
    assert response.text == expected


@pytest.mark.parametrize(
    "invalid_smiles, toolkit, exception_response_code", [("INVALID_SMILES", "cdk", 400)]
)
def test_exception_tanimoto_similarity(
    invalid_smiles, toolkit, exception_response_code
):
    response = client.get(f"/latest/chem/tanimoto?smiles={invalid_smiles}&toolkit={toolkit}")
    assert response.status_code == exception_response_code


@pytest.mark.parametrize(
    "smiles, fix, expected",
    [
        ("CCO", False, '{"smi":"CCO","messages":["No Errors Found"]}'),
        ("CCO", True, '{"smi":"CCO","messages":["No Errors Found"]}'),
    ],
)
def test_successful_check_errors(smiles, fix, expected):
    response = client.get(f"/latest/chem/errors?smiles={smiles}&fix={fix}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"
    assert response.text == expected


@pytest.mark.parametrize(
    "invalid_smiles, fix, exception_response_code", [("INVALID_SMILES", False, 400)]
)
def test_exception_check_errors(invalid_smiles, fix, exception_response_code):
    response = client.get(f"/latest/chem/errors?smiles={invalid_smiles}&fix={fix}")
    assert response.status_code == exception_response_code


@pytest.mark.parametrize(
    "smiles, toolkit, ringsize, success_response_code, response_text",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            False,
            200,
            '["C-4;(//)","N-3;(//)","C-3;(//)","N-2;(//)","C-3;(//)","C-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-4;(//)","C-4;(//)"]',
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            True,
            200,
            '["C-4;(//)","N-3;(//)","C-3;(//)","N-2;(//)","C-3;(//)","C-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-3;(//)","O-1;(//)","N-3;(//)","C-4;(//)"]',
        ),
    ],
)
def test_successful_hose_codes(
    smiles, toolkit, ringsize, success_response_code, response_text
):
    response = client.get(
        f"/latest/chem/HOSEcode?smiles={smiles}&spheres=0&toolkit={toolkit}&ringsize={ringsize}"
    )
    assert response.status_code == {success_response_code}
    assert response.headers["content-type"] == "application/json"
    assert response.text == {response_text}


@pytest.mark.parametrize(
    "invalid_smiles, toolkit, ringsize, exception_response_code",
    [("INVALID_SMILES", "cdk", False, 500)],
)
def test_exception_hose_codes(
    invalid_smiles, toolkit, ringsize, exception_response_code
):
    response = client.get(
        f"/latest/chem/HOSEcode?smiles={invalid_smiles}&spheres=0&toolkit={toolkit}&ringsize={ringsize}"
    )
    assert response.status_code == exception_response_code


def test_success_standardize_mol(molfile):
    response = client.post(
        "/latest/chem/standardize", data=molfile, headers={"Content-Type": "text/plain"}
    )
    assert response.status_code == 200
    assert "standardized_mol" in response.json()
    assert "canonical_smiles" in response.json()
    assert "inchi" in response.json()
    assert "inchikey" in response.json()


@pytest.mark.parametrize(
    "invalid_molfile, exception_response_code", [(" CDK     08302311362D", 400)]
)
def test_exception_standardize_mol(invalid_molfile, exception_response_code):
    response = client.post(
        "/latest/chem/standardize",
        data=invalid_molfile,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == exception_response_code


def test_successful_coconut_preprocessing(test_smiles):
    response = client.get(f"/latest/chem/coconut/pre-processing?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


def test_exception_coconut_preprocessing():
    invalid_smiles = "INVALID_SMILES"
    response = client.get(
        f"/latest/chem/coconut/pre-processing?smiles={invalid_smiles}"
    )
    assert response.status_code == 400


# Run the tests
if __name__ == "__main__":
    pytest.main()
