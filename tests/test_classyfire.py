import pytest
from unittest.mock import patch, MagicMock, AsyncMock
import httpx
from app.modules.classyfire import classify, result


@pytest.fixture
def valid_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def invalid_smiles():
    return "invalid_smiles"


@pytest.mark.asyncio
@patch("app.modules.classyfire.httpx.AsyncClient")
async def test_valid_classyfire(mock_async_client, valid_smiles):
    # Create mock client instance and context manager
    mock_client_instance = MagicMock()
    mock_async_client.return_value.__aenter__.return_value = mock_client_instance

    # Mock the POST response for classify
    mock_post_response = MagicMock()
    mock_post_response.json.return_value = {
        "id": "12345",
        "query_type": "STRUCTURE",
        "query_input": valid_smiles,
    }
    mock_post_response.raise_for_status.return_value = None
    mock_client_instance.post = AsyncMock(return_value=mock_post_response)

    # Mock the GET response for result
    mock_get_response = MagicMock()
    mock_get_response.json.return_value = {
        "id": "12345",
        "classification_status": "Done",
        "entities": [{"class": {"name": "Imidazopyrimidines"}}],
    }
    mock_get_response.raise_for_status.return_value = None
    mock_client_instance.get = AsyncMock(return_value=mock_get_response)

    result_ = await classify(valid_smiles)
    assert result_["query_type"] == "STRUCTURE"
    id_ = result_["id"]

    classified = await result(id_)
    assert classified["classification_status"] == "Done"
    assert classified["entities"][0]["class"]["name"] == "Imidazopyrimidines"


@pytest.mark.asyncio
@patch("app.modules.classyfire.httpx.AsyncClient")
async def test_invalid_classyfire(mock_async_client, invalid_smiles):
    # Create mock client instance and context manager
    mock_client_instance = MagicMock()
    mock_async_client.return_value.__aenter__.return_value = mock_client_instance

    # Mock the POST response for classify
    mock_post_response = MagicMock()
    mock_post_response.json.return_value = {
        "id": "12346",
        "query_type": "STRUCTURE",
        "query_input": invalid_smiles,
    }
    mock_post_response.raise_for_status.return_value = None
    mock_client_instance.post = AsyncMock(return_value=mock_post_response)

    # Mock the GET response for result
    mock_get_response = MagicMock()
    mock_get_response.json.return_value = {
        "id": "12346",
        "classification_status": "Done",
        "invalid_entities": [
            {"report": ["Cannot process the input SMILES string, please check again"]}
        ],
    }
    mock_get_response.raise_for_status.return_value = None
    mock_client_instance.get = AsyncMock(return_value=mock_get_response)

    result_ = await classify(invalid_smiles)
    assert result_["query_input"] == "invalid_smiles"
    id_ = result_["id"]

    classified = await result(id_)
    assert classified["classification_status"] == "Done"
    assert (
        classified["invalid_entities"][0]["report"][0]
        == "Cannot process the input SMILES string, please check again"
    )


@pytest.mark.asyncio
@patch("app.modules.classyfire.httpx.AsyncClient")
async def test_classify_http_error(mock_async_client, valid_smiles):
    """Test that classify properly raises HTTPError exceptions"""
    # Create mock client instance and context manager
    mock_client_instance = MagicMock()
    mock_async_client.return_value.__aenter__.return_value = mock_client_instance
    
    # Mock the POST to raise an HTTPError
    mock_client_instance.post = AsyncMock(
        side_effect=httpx.HTTPError("Connection failed")
    )
    
    # Should re-raise the HTTPError
    with pytest.raises(httpx.HTTPError):
        await classify(valid_smiles)


@pytest.mark.asyncio
@patch("app.modules.classyfire.httpx.AsyncClient")
async def test_classify_timeout_error(mock_async_client, valid_smiles):
    """Test that classify properly raises timeout errors"""
    # Create mock client instance and context manager
    mock_client_instance = MagicMock()
    mock_async_client.return_value.__aenter__.return_value = mock_client_instance
    
    # Mock the POST to raise a timeout error
    mock_client_instance.post = AsyncMock(
        side_effect=httpx.TimeoutException("Request timed out")
    )
    
    # Should re-raise the HTTPError (TimeoutException is a subclass of HTTPError)
    with pytest.raises(httpx.HTTPError):
        await classify(valid_smiles)


@pytest.mark.asyncio
@patch("app.modules.classyfire.httpx.AsyncClient")
async def test_result_http_error(mock_async_client):
    """Test that result properly raises HTTPError exceptions"""
    # Create mock client instance and context manager
    mock_client_instance = MagicMock()
    mock_async_client.return_value.__aenter__.return_value = mock_client_instance
    
    # Mock the GET to raise an HTTPError
    mock_client_instance.get = AsyncMock(
        side_effect=httpx.HTTPError("Connection failed")
    )
    
    # Should re-raise the HTTPError
    with pytest.raises(httpx.HTTPError):
        await result("12345")


@pytest.mark.asyncio
@patch("app.modules.classyfire.httpx.AsyncClient")
async def test_result_not_found_error(mock_async_client):
    """Test that result properly raises HTTPStatusError for 404"""
    # Create mock client instance and context manager
    mock_client_instance = MagicMock()
    mock_async_client.return_value.__aenter__.return_value = mock_client_instance
    
    # Create a mock response for 404
    mock_response = MagicMock()
    mock_response.status_code = 404
    
    # Mock the GET to raise an HTTPStatusError
    mock_client_instance.get = AsyncMock(
        side_effect=httpx.HTTPStatusError(
            "Not Found",
            request=MagicMock(),
            response=mock_response
        )
    )
    
    # Should re-raise the HTTPError
    with pytest.raises(httpx.HTTPError):
        await result("invalid_id")

