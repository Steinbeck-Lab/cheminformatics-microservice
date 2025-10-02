from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.main import app


client = TestClient(app)


def test_tool_index():
    response = client.get("/latest/tools/")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}


def test_generate_structures_exception():
    molecular_formula = "invalid_formula"
    response = client.get(
        f"/latest/tools/generate-structures?molecular_formula={molecular_formula}",
    )
    # Should return 500 for invalid formula that causes surge to fail
    assert response.status_code == 500


@pytest.mark.parametrize(
    "input, response_code",
    [
        ("C6H6", 200),
        ("C4H8", 200),
    ],
)
def test_generate_structures(input, response_code):
    response = client.get(
        f"/latest/tools/generate-structures?molecular_formula={input}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"

    # Check the new response format
    response_data = response.json()
    assert "message" in response_data
    assert "output" in response_data
    assert response_data["message"] == "Success"

    # Check output structure
    output = response_data["output"]
    assert "total_count" in output
    assert "generated_count" in output
    assert "structures" in output
    assert "settings" in output
    assert "formula" in output
    assert "limit_applied" in output

    # Check data types
    assert isinstance(output["total_count"], int)
    assert isinstance(output["generated_count"], int)
    assert isinstance(output["structures"], list)
    assert isinstance(output["settings"], dict)
    assert isinstance(output["formula"], str)
    assert isinstance(output["limit_applied"], bool)

    # Check that we have structures
    assert len(output["structures"]) > 0
    assert output["formula"] == input


def test_generate_structures_heavy_atom_limit():
    # Test formula with more than 10 heavy atoms (C15H32 has 15 carbon atoms)
    molecular_formula = "C15H32"
    response = client.get(
        f"/latest/tools/generate-structures?molecular_formula={molecular_formula}",
    )
    assert response.status_code == 400
    assert "heavy atoms" in response.json()["detail"]


@pytest.mark.parametrize(
    "input,response_text, response_code",
    [
        (
            "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
            '"The molecule contains Linear and Circular sugars"',
            200,
        ),
        (
            "C1=CC=C(C(=C1)CO)OC2C(C(C(C(O2)CO)O)O)O",
            '"The molecule contains only Circular sugar"',
            200,
        ),
        (
            "CC1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=CC=C3OCC(C(C(C(C(=O)OCC(C(C(C(C=O)O)O)(O)OCC(C(C(C(C=O)O)O)O)O)O)O)O)O)O",
            '"The molecule contains only Linear sugar"',
            200,
        ),
        (
            "CC",
            '"The molecule contains no sugar"',
            200,
        ),
        ("INVALID_INPUT", "", 422),
    ],
)
def test_sugars_info(input, response_text, response_code):
    response = client.get(f"/latest/tools/sugars-info?smiles={input}")
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if input != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "input,response_text, response_code",
    [
        ("OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O", '"C(C1C(C(C(CO1)O)O)O)O"', 200),
        ("INVALID_INPUT", "", 422),
    ],
)
def test_remove_linear_sugars(input, response_text, response_code):
    response = client.get(f"/latest/tools/remove-linear-sugars?smiles={input}")
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if input != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "input,response_text, response_code",
    [
        (
            "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
            '"C(C(C(C(CO)O)O)O)O"',
            200,
        ),
        ("INVALID_INPUT", "", 422),
    ],
)
def test_remove_circular_sugars(input, response_text, response_code):
    response = client.get(
        f"/latest/tools/remove-circular-sugars?smiles={input}",
    )
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if input != "INVALID_INPUT":
        assert response.text == response_text


@pytest.mark.parametrize(
    "input,response_text, response_code",
    [
        (
            "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
            '"C1=C(C=C(C(=C1)O)O)C2C3=C(C=C(C=O)C2C(=O)OC4CC(=CC(C4O)O)C(=O)O)C=C(C(=C3)O)O"',
            200,
        ),
        ("INVALID_INPUT", "", 422),
    ],
)
def test_remove_sugars(input, response_text, response_code):
    response = client.get(f"/latest/tools/remove-sugars?smiles={input}")
    assert response.status_code == response_code
    assert response.headers["content-type"] == "application/json"
    if input != "INVALID_INPUT":
        assert response.text == response_text
