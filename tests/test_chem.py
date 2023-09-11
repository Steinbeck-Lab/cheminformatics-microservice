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


@pytest.mark.parametrize(
    "smiles, response_text, response_code",
    [
        (
            "CC",
            '["CC"]',
            200,
        ),
        (
            "INVALID_INPUT",
            "",
            422,
        ),
    ],
)
def test_smiles_to_stereo_isomers(smiles, response_text, response_code):
    response = client.get(f"/latest/chem/stereoisomers?smiles={smiles}")
    assert response.status_code == response_code
    if smiles != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "smiles, format, response_code",
    [
        (
            "CC",
            "html",
            200,
        ),
        (
            "CC",
            "json",
            200,
        ),
        (
            "INVALID_INPUT",
            "json",
            422,
        ),
    ],
)
def test_smiles_descriptors(smiles, format, response_code):
    response = client.get(f"/latest/chem/descriptors?smiles={smiles}&format={format}")
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "multiple_smiles, toolkit, response_code",
    [
        (
            "CC,CCO",
            "cdk",
            200,
        ),
        (
            "CC,CCO",
            "rdkit",
            200,
        ),
        ("CC", "cdk", 422),
    ],
)
def test_smiles_descriptors_multiple(multiple_smiles, toolkit, response_code):
    response = client.get(
        f"/latest/chem/descriptors/multiple?smiles={multiple_smiles}&toolkit={toolkit}"
    )
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "smiles, expected_score, response_code",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", -1.09, 200),
        ("CC=O", 0.96, 200),
        ("INVALID_INPUT", 0, 422),
    ],
)
def test_NPlikeliness_Score(smiles, expected_score, response_code):
    response = client.get(f"/latest/chem/nplikeness/score?smiles={smiles}")
    assert response.status_code == response_code
    if smiles != "INVALID_INPUT":
        assert response.json() == expected_score


def test_successful_classyFire_classify(test_smiles):
    response = client.get(f"/latest/chem/classyfire/classify?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


def test_successful_classyFire_result():
    job_id = 11212508
    response = client.get(f"/latest/chem/classyfire/{job_id}/result")
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, expected, response_code",
    [
        ("CC,CC", "cdk", '"1.00000"', 200),
        ("CC,CC", "rdkit", "1.0", 200),
        (
            "CC,CC,CC",
            "cdk",
            "<table><tr><th></th><th>0</th><th>1</th><th>2</th></tr><tr><td>0</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr><tr><td>1</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr><tr><td>2</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr></table>",
            200,
        ),
        ("INVALID_INPUT", "cdk", "", 422),
    ],
)
def test_successful_tanimoto_similarity(smiles, toolkit, expected, response_code):
    response = client.get(f"/latest/chem/tanimoto?smiles={smiles}&toolkit={toolkit}")
    assert response.status_code == response_code
    if smiles != "INVALID_INPUT":
        assert response.text == expected


@pytest.mark.parametrize(
    "smiles, fix, expected, response_code",
    [
        ("CCO", False, '{"smi":"CCO","messages":["No Errors Found"]}', 200),
        ("CCO", True, '{"smi":"CCO","messages":["No Errors Found"]}', 200),
        ("INVALID_INPUT", False, "", 422),
    ],
)
def test_successful_check_errors(smiles, fix, expected, response_code):
    response = client.get(f"/latest/chem/errors?smiles={smiles}&fix={fix}")
    assert response.status_code == response_code
    if smiles != "INVALID_INPUT":
        assert response.text == expected


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
    assert response.status_code == success_response_code
    assert response.text == response_text


# @pytest.mark.parametrize(
#     "invalid_smiles, toolkit, ringsize, exception_response_code",
#     [("test", "cdk", False, 500)],
# )
# def test_exception_hose_codes(
#     invalid_smiles, toolkit, ringsize, exception_response_code
# ):
#     response = client.get(
#         f"/latest/chem/HOSEcode?smiles={invalid_smiles}&spheres=0&toolkit={toolkit}&ringsize={ringsize}"
#     )
#     assert response.status_code == exception_response_code


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
    "invalid_molfile, exception_response_code", [(" CDK     08302311362D", 422)]
)
def test_exception_standardize_mol(invalid_molfile, exception_response_code):
    response = client.post(
        "/latest/chem/standardize",
        data=invalid_molfile,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == exception_response_code


@pytest.mark.parametrize(
    "smiles, response_code",
    [("CCO", 200), ("INVALID_INPUT", 422)],
)
def test_successful_coconut_preprocessing(smiles, response_code):
    response = client.get(f"/latest/chem/coconut/pre-processing?smiles={smiles}")
    assert response.status_code == response_code


# Run the tests
if __name__ == "__main__":
    pytest.main()
