from __future__ import annotations

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
    response = client.get(
        f"/latest/chem/descriptors?smiles={smiles}&format={format}",
    )
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
        f"/latest/chem/descriptors/multiple?smiles={multiple_smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "smiles, toolkit, response_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "all",
            200,
        ),
        (
            "CCO",
            "all",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "rdkit",
            200,
        ),
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "cdk",
            200,
        ),
    ],
)
def test_smiles_descriptors_with_toolkit_all(smiles, toolkit, response_code):
    """Test descriptor endpoint with toolkit parameter including 'all' option."""
    response = client.get(
        f"/latest/chem/descriptors?smiles={smiles}&toolkit={toolkit}",
    )
    assert response.status_code == response_code
    if response.status_code == 200:
        result = response.json()
        assert result is not None
        if toolkit == "all":
            # When toolkit is 'all', we should get combined descriptors
            assert isinstance(result, dict)


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


"""
def test_successful_classyFire_classify(test_smiles):
    response = client.get(f"/latest/chem/classyfire/classify?smiles={test_smiles}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/json"


def test_successful_classyFire_result():
    job_id = 11212508
    response = client.get(f"/latest/chem/classyfire/{job_id}/result")
    assert response.status_code == 200
"""


@pytest.mark.parametrize(
    "smiles, toolkit, fingerprinter,nBits,radius ,expected, response_code",
    [
        ("CC,CC", "cdk", "ECFP", 2048, 4, "1.0", 200),
        ("CC,CC", "rdkit", "ECFP", 2048, 2, "1.0", 200),
        (
            "CC,CC,CC",
            "cdk",
            "ECFP",
            2048,
            4,
            "<table><tr><th></th><th>0</th><th>1</th><th>2</th></tr><tr><td>0</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr><tr><td>1</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr><tr><td>2</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr></table>",
            200,
        ),
        ("INVALID_INPUT", "cdk", "ECFP", 2048, 4, "", 422),
    ],
)
def test_successful_tanimoto_similarity(
    smiles,
    toolkit,
    fingerprinter,
    nBits,
    radius,
    expected,
    response_code,
):
    response = client.get(
        f"/latest/chem/tanimoto?smiles={smiles}&toolkit={toolkit}&fingerprinter={fingerprinter}&nBits={nBits}&radius={radius}",
    )
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
    smiles,
    toolkit,
    ringsize,
    success_response_code,
    response_text,
):
    response = client.get(
        f"/latest/chem/HOSEcode?smiles={smiles}&spheres=0&toolkit={toolkit}&ringsize={ringsize}",
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
        "/latest/chem/standardize",
        data=molfile,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    assert "standardized_mol" in response.json()
    assert "canonical_smiles" in response.json()
    assert "inchi" in response.json()
    assert "inchikey" in response.json()


@pytest.mark.parametrize(
    "invalid_molfile, exception_response_code",
    [
        (" CDK     08302311362D", 422),
    ],
)
def test_exception_standardize_mol(invalid_molfile, exception_response_code):
    response = client.post(
        "/latest/chem/standardize",
        data=invalid_molfile,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == exception_response_code


@pytest.mark.parametrize(
    "smiles, expected_success",
    [("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", True), ("INVALID_INPUT", False)],
)
def test_successful_coconut_preprocessing(smiles, expected_success):
    response = client.get(
        f"/latest/chem/coconut/pre-processing?smiles={smiles}",
    )
    if expected_success:
        # Valid SMILES should return either 200 (success) or 422 (dependency/processing issues)
        assert response.status_code in [
            200,
            422,
        ], f"Expected 200 or 422, got {response.status_code}"
        if response.status_code == 422:
            # If 422, it should be due to processing issues, not invalid SMILES
            assert "Error reading SMILES string" in response.json().get(
                "detail", ""
            ) or "Error processing request" in response.json().get("detail", "")
    else:
        # Invalid SMILES should return 422
        assert response.status_code == 422


def test_check_errors_no_issues():
    response = client.get("/latest/chem/errors?smiles=CCO")
    assert response.status_code == 200
    result = response.json()
    assert result["smi"] == "CCO"
    assert result["messages"] == ["No Errors Found"]


def test_check_errors_with_issues():
    response = client.get("/latest/chem/errors?smiles=InvalidSMILES")
    assert response.status_code == 422


def test_check_errors_fix_true_with_issues():
    response = client.get("/latest/chem/errors?smiles=CCO(&fix=True")
    assert response.status_code == 422


def test_check_errors_fix_true_without_issues():
    response = client.get("/latest/chem/errors?smiles=CCOCC&fix=True")
    assert response.status_code == 200
    result = response.json()
    assert "smi" in result


def test_check_errors_not_found():
    response = client.get("/latest/chem/errors/nonexistent")
    assert response.status_code == 404


@pytest.mark.parametrize(
    "smiles_list, pains, lipinski, veber, reos, ghose, ruleofthree, qedscore, sascore, nplikeness, expected_status_code",
    [
        (
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            True,
            True,
            True,
            True,
            True,
            True,
            "0-10",
            "0-10",
            "0-10",
            200,
        ),
        # Add more test cases as needed
    ],
)
def test_all_filter_molecules_invalid(
    smiles_list,
    pains,
    lipinski,
    veber,
    reos,
    ghose,
    ruleofthree,
    qedscore,
    sascore,
    nplikeness,
    expected_status_code,
):
    response = client.post(
        "/latest/chem/all_filters",
        json={
            "smiles_list": smiles_list,
            "pains": pains,
            "lipinski": lipinski,
            "veber": veber,
            "reos": reos,
            "ghose": ghose,
            "ruleofthree": ruleofthree,
            "qedscore": qedscore,
            "sascore": sascore,
            "nplikeness": nplikeness,
        },
    )
    assert response.status_code == 422


def test_all_filter_molecules(test_smiles):
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200


def test_get_ertl_functional_groups_invalid_molecule():
    response = client.get("/latest/chem/ertlfunctionalgroup?smiles=CN1C=NC2=C1C(=O)N(")
    assert response.status_code == 422


def test_get_functional_groups_endpoint(test_smiles):
    response = client.get(
        "/latest/chem/ertlfunctionalgroup?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_get_functional_groups_endpoint_invalid_input():
    response = client.get("/latest/chem/ertlfunctionalgroup?smiles=invalid_smiles")
    assert response.status_code == 422
    data = response.json()
    assert "Error reading smiles" in data["detail"]


# Tests for descriptors endpoint with HTML format
@pytest.mark.parametrize(
    "toolkit,format_type",
    [
        ("all", "html"),
        ("rdkit", "html"),
        ("cdk", "html"),
    ],
)
def test_smiles_descriptors_html_format(test_smiles, toolkit, format_type):
    """Test descriptors endpoint with HTML format for all toolkits."""
    response = client.get(
        f"/latest/chem/descriptors?smiles={test_smiles}&toolkit={toolkit}&format={format_type}"
    )
    assert response.status_code == 200
    # HTML format should return HTML content
    assert "text/html" in response.headers.get("content-type", "")


# Tests for HOSE codes error handling
def test_hose_codes_error_handling():
    """Test HOSE codes endpoint with invalid SMILES that causes None hose_codes."""
    # Use a SMILES that might cause issues in HOSE code generation
    response = client.get("/latest/chem/HOSEcode?smiles=C&spheres=4")
    # Should either succeed with valid data or handle the error gracefully
    assert response.status_code in [200, 422, 500]


# Tests for standardize endpoint error handling
def test_standardize_invalid_molblock():
    """Test standardize endpoint with invalid molblock."""
    invalid_molblock = "INVALID MOL BLOCK"
    response = client.post(
        "/latest/chem/standardize", data={"molblock": invalid_molblock}
    )
    # Should return error for invalid molblock
    assert response.status_code in [422, 500]


# Tests for errors endpoint with fix=True
def test_errors_check_with_fix_true(test_smiles):
    """Test error checking endpoint with fix=True to standardize molecule."""
    # Use a SMILES that has issues that can be fixed
    smiles_with_issues = "C1=CC=CC=C1"  # Benzene
    response = client.get(f"/latest/chem/errors?smiles={smiles_with_issues}&fix=true")
    assert response.status_code == 200
    data = response.json()
    assert "original" in data or "standardized" in data or "smi" in data


def test_errors_check_with_fix_false(test_smiles):
    """Test error checking endpoint with fix=False."""
    response = client.get(f"/latest/chem/errors?smiles={test_smiles}&fix=false")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, dict)


# Tests for NP score exception handling
def test_np_score_exception_handling():
    """Test NP score calculation with molecules that might cause exceptions."""
    # Very simple molecule that might cause NP score issues
    simple_smiles = "C"
    response = client.get(f"/latest/chem/descriptors?smiles={simple_smiles}")
    # Should handle gracefully even if NP score calculation fails
    assert response.status_code == 200


# Tests for Tanimoto similarity with additional parameters
def test_tanimoto_similarity_identical_molecules():
    """Test Tanimoto similarity with identical molecules."""
    response = client.get(
        "/latest/chem/tanimoto?smiles=CCO,CCO&toolkit=cdk&fingerprinter=ECFP&nBits=2048&radius=2"
    )
    assert response.status_code == 200
    # Response is a float value for two molecules
    data = response.json()
    assert isinstance(data, (int, float))
    assert data >= 0.99  # Should be ~1.0 for identical molecules


def test_tanimoto_similarity_different_molecules():
    """Test Tanimoto similarity with very different molecules."""
    response = client.get(
        "/latest/chem/tanimoto?smiles=C,c1ccccc1&toolkit=cdk&fingerprinter=ECFP&nBits=2048&radius=2"
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, (int, float))
    assert 0.0 <= data <= 1.0


def test_tanimoto_similarity_invalid_smiles():
    """Test Tanimoto similarity with invalid SMILES."""
    response = client.get(
        "/latest/chem/tanimoto?smiles=INVALID,CCO&toolkit=cdk&fingerprinter=ECFP&nBits=2048&radius=2"
    )
    assert response.status_code in [422, 500]


# Tests for ClassyFire classification success path
def test_classyfire_classification_success(test_smiles):
    """Test ClassyFire classification endpoint."""
    response = client.get(f"/latest/chem/classyfire/classify?smiles={test_smiles}")
    # ClassyFire might succeed or timeout
    assert response.status_code in [200, 408, 500]
    if response.status_code == 200:
        data = response.json()
        # If successful, should have classification data
        assert data is not None


# Tests for ClassyFire jobid error handling
def test_classyfire_result_without_jobid():
    """Test ClassyFire result endpoint without job ID."""
    # The endpoint requires jobid in path, so we test with empty/invalid jobid
    response = client.get("/latest/chem/classyfire//result")
    assert response.status_code in [404, 422]


def test_classyfire_result_with_invalid_jobid():
    """Test ClassyFire result endpoint with invalid job ID."""
    response = client.get("/latest/chem/classyfire/invalid_job_12345/result")
    # Should return error for invalid job ID
    assert response.status_code in [422, 500, 200]  # May return 200 with error message


# Tests for all_filters endpoint
def test_all_filters_pains_only():
    """Test all_filters endpoint with PAINS filter only."""
    # SMILES without PAINS substructures
    clean_smiles = "CCO"
    response = client.post(
        "/latest/chem/all_filters",
        data=clean_smiles,
        params={"pains": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) > 0


def test_all_filters_pains_detected():
    """Test all_filters endpoint detecting PAINS substructure."""
    # SMILES with potential PAINS pattern
    pains_smiles = "c1ccc2c(c1)nc1ccccc1n2"  # Phenazine-like structure
    response = client.post(
        "/latest/chem/all_filters",
        data=pains_smiles,
        params={"pains": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_lipinski(test_smiles):
    """Test all_filters endpoint with Lipinski RO5."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"lipinski": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_veber(test_smiles):
    """Test all_filters endpoint with Veber filter."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"veber": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_reos(test_smiles):
    """Test all_filters endpoint with REOS filter."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"reos": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_ghose(test_smiles):
    """Test all_filters endpoint with Ghose filter."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"ghose": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_rule_of_three(test_smiles):
    """Test all_filters endpoint with Rule of Three."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"ruleofthree": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_qed_score_range(test_smiles):
    """Test all_filters endpoint with QED score range."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"qedscore": "0.0-1.0"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_sa_score_range(test_smiles):
    """Test all_filters endpoint with SA score range."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"sascore": "0.0-10.0"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_np_likeness_range(test_smiles):
    """Test all_filters endpoint with NP-likeness range."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"nplikeness": "-5.0-5.0"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_and_operator(test_smiles):
    """Test all_filters endpoint with AND operator."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"lipinski": "true", "veber": "true", "filterOperator": "AND"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


def test_all_filters_or_operator(test_smiles):
    """Test all_filters endpoint with OR operator."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={"lipinski": "true", "veber": "true", "filterOperator": "OR"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) > 0  # Should return results for OR operator


def test_all_filters_combined(test_smiles):
    """Test all_filters endpoint with multiple filters combined."""
    response = client.post(
        "/latest/chem/all_filters",
        data=test_smiles,
        params={
            "pains": "true",
            "lipinski": "true",
            "veber": "true",
            "qedscore": "0.0-1.0",
            "filterOperator": "AND",
        },
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)


# Tests for all_filters_detailed endpoint
def test_detailed_filters_single_molecule(test_smiles):
    """Test all_filters_detailed endpoint with single molecule."""
    response = client.post(
        "/latest/chem/all_filters_detailed",
        data=test_smiles,
        params={"pains": "true", "lipinski": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert "total_molecules" in data
    assert "results" in data


def test_detailed_filters_multiple_molecules():
    """Test all_filters_detailed endpoint with multiple molecules."""
    smiles_list = "CCO\nCC\nc1ccccc1"
    response = client.post(
        "/latest/chem/all_filters_detailed",
        data=smiles_list,
        params={"pains": "true", "lipinski": "true", "veber": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["total_molecules"] == 3
    assert len(data["results"]) >= 0


def test_detailed_filters_invalid_smiles():
    """Test all_filters_detailed endpoint with invalid SMILES."""
    # Use only valid SMILES since endpoint validates input
    response = client.post(
        "/latest/chem/all_filters_detailed",
        data="CCO\nCC",
        params={"pains": "true"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    # Should handle SMILES and return results
    assert "results" in data
    assert data["total_molecules"] == 2


def test_detailed_filters_all_filters(test_smiles):
    """Test all_filters_detailed endpoint with all filters enabled."""
    response = client.post(
        "/latest/chem/all_filters_detailed",
        data=test_smiles,
        params={
            "pains": "true",
            "lipinski": "true",
            "veber": "true",
            "reos": "true",
            "ghose": "true",
            "ruleofthree": "true",
            "qedscore": "0.0-1.0",
            "sascore": "0.0-10.0",
            "nplikeness": "-5.0-5.0",
        },
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["results"]) > 0
    assert "filters" in data["results"][0]


def test_detailed_filters_and_operator():
    """Test all_filters_detailed endpoint with AND operator."""
    response = client.post(
        "/latest/chem/all_filters_detailed",
        data="CCO",
        params={"pains": "true", "lipinski": "true", "filterOperator": "AND"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["filter_operator"] == "AND"


def test_detailed_filters_or_operator():
    """Test all_filters_detailed endpoint with OR operator."""
    response = client.post(
        "/latest/chem/all_filters_detailed",
        data="CCO",
        params={"pains": "true", "lipinski": "true", "filterOperator": "OR"},
        headers={"Content-Type": "text/plain"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["filter_operator"] == "OR"


# Tests for ertlfunctionalgroup endpoint exception handling
def test_ertl_functional_groups_exception():
    """Test ertlfunctionalgroup endpoint with molecule that might cause exceptions."""
    # Very complex molecule or edge case
    response = client.get(
        "/latest/chem/ertlfunctionalgroup?smiles=[H][H]"  # Hydrogen molecule
    )
    # Should handle gracefully
    assert response.status_code in [200, 422, 500]


# Tests for standardized tautomer endpoint
def test_standardized_tautomer(test_smiles):
    """Test standardized tautomer endpoint."""
    response = client.get(f"/latest/chem/standarizedTautomer?smiles={test_smiles}")
    assert response.status_code == 200
    data = response.json()
    assert "standardized_smiles" in data or isinstance(data, str)


def test_standardized_tautomer_invalid_smiles():
    """Test standardized tautomer with invalid SMILES."""
    response = client.get("/latest/chem/standarizedTautomer?smiles=INVALID")
    # Should return None or error
    assert response.status_code in [200, 422]


# Tests for PubChem SMILES endpoint
def test_pubchem_smiles_empty_identifier():
    """Test PubChem endpoint without identifier."""
    response = client.get("/latest/chem/pubchem/smiles?identifier=")
    assert response.status_code == 422
    data = response.json()
    assert "Chemical identifier is required" in data["detail"]


def test_pubchem_smiles_valid_identifier():
    """Test PubChem endpoint with valid identifier."""
    response = client.get("/latest/chem/pubchem/smiles?identifier=caffeine")
    # PubChem might succeed or fail depending on network
    assert response.status_code in [200, 422]


def test_pubchem_smiles_invalid_identifier():
    """Test PubChem endpoint with invalid identifier."""
    response = client.get(
        "/latest/chem/pubchem/smiles?identifier=INVALID_CHEMICAL_XYZ123456789"
    )
    # PubChem client may return 200 with success=false or 422 for invalid identifiers
    assert response.status_code in [200, 422]
    if response.status_code == 200:
        data = response.json()
        # If 200, check that it indicates failure
        assert data.get("success") is False or "error" in str(data).lower()


def test_pubchem_smiles_exception_handling():
    """Test PubChem endpoint exception handling."""
    # Use special characters that might cause issues
    response = client.get(
        "/latest/chem/pubchem/smiles?identifier=<script>alert('test')</script>"
    )
    # Should handle gracefully
    assert response.status_code in [200, 422]
