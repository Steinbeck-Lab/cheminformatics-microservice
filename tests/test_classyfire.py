import pytest
from unittest.mock import patch, MagicMock
from app.modules.classyfire import classify, result


@pytest.fixture
def valid_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def invalid_smiles():
    return "invalid_smiles"


@pytest.mark.asyncio
@patch("app.modules.classyfire.requests.post")
@patch("app.modules.classyfire.requests.get")
async def test_valid_classyfire(mock_get, mock_post, valid_smiles):
    # Mock the initial classification request
    mock_post_response = MagicMock()
    mock_post_response.json.return_value = {
        "id": "12345",
        "query_type": "STRUCTURE",
        "query_input": valid_smiles,
    }
    mock_post_response.raise_for_status.return_value = None
    mock_post.return_value = mock_post_response

    # Mock the result retrieval request
    mock_get_response = MagicMock()
    mock_get_response.json.return_value = {
        "id": "12345",
        "classification_status": "Done",
        "entities": [{"class": {"name": "Imidazopyrimidines"}}],
    }
    mock_get_response.raise_for_status.return_value = None
    mock_get.return_value = mock_get_response

    result_ = await classify(valid_smiles)
    assert result_["query_type"] == "STRUCTURE"
    id_ = result_["id"]

    classified = await result(id_)
    assert classified["classification_status"] == "Done"
    assert classified["entities"][0]["class"]["name"] == "Imidazopyrimidines"


@pytest.mark.asyncio
@patch("app.modules.classyfire.requests.post")
@patch("app.modules.classyfire.requests.get")
async def test_invalid_classyfire(mock_get, mock_post, invalid_smiles):
    # Mock the initial classification request
    mock_post_response = MagicMock()
    mock_post_response.json.return_value = {
        "id": "12346",
        "query_type": "STRUCTURE",
        "query_input": invalid_smiles,
    }
    mock_post_response.raise_for_status.return_value = None
    mock_post.return_value = mock_post_response

    # Mock the result retrieval request
    mock_get_response = MagicMock()
    mock_get_response.json.return_value = {
        "id": "12346",
        "classification_status": "Done",
        "invalid_entities": [
            {"report": ["Cannot process the input SMILES string, please check again"]}
        ],
    }
    mock_get_response.raise_for_status.return_value = None
    mock_get.return_value = mock_get_response

    result_ = await classify(invalid_smiles)
    assert result_["query_input"] == "invalid_smiles"
    id_ = result_["id"]

    classified = await result(id_)
    assert classified["classification_status"] == "Done"
    assert (
        classified["invalid_entities"][0]["report"][0]
        == "Cannot process the input SMILES string, please check again"
    )
