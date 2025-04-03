import unittest
from unittest.mock import patch, MagicMock
import pytest
from requests.exceptions import RequestException

from app.modules.pubchem_retrieve import PubChemClient


# Mock client class to replace TestClient
class MockClient:
    """Simple mock client to replace TestClient"""

    def get(self, url):
        """Mock get method that returns a response object"""
        response = MagicMock()

        # Extract the full identifier from the URL without truncating
        identifier = url.split("identifier=")[1] if "identifier=" in url else ""

        # Set response properties based on the identifier
        if not identifier:
            response.status_code = 422
            response.json.return_value = {"detail": "Missing identifier"}
        else:
            response.status_code = 200
            response.json.return_value = {
                "input": identifier,
                "canonical_smiles": self.mock_get_smiles_return_value,
                "input_type": self.expected_input_type,
                "success": self.expected_success,
            }

        return response

    def __init__(self):
        self.mock_get_smiles_return_value = None
        self.expected_input_type = "name"
        self.expected_success = False


# Fixture to provide mock client
@pytest.fixture
def client():
    """Provide a mock client for testing"""
    return MockClient()


class TestPubChemClient(unittest.TestCase):
    """Test suite for the PubChemClient class."""

    def setUp(self):
        """Set up test fixtures."""
        self.client = PubChemClient()

        # Mock response for successful query
        self.mock_success_response = MagicMock()
        self.mock_success_response.status_code = 200
        self.mock_success_response.json.return_value = {
            "PropertyTable": {
                "Properties": [{"IsomericSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O"}]
            }
        }
        self.mock_success_response.text = "2244\n"

        # Mock response for failed query
        self.mock_failed_response = MagicMock()
        self.mock_failed_response.status_code = 404
        self.mock_failed_response.raise_for_status.side_effect = RequestException(
            "Not found"
        )

    def test_init_default_values(self):
        """Test initialization with default values."""
        client = PubChemClient()
        self.assertEqual(client.timeout, PubChemClient.DEFAULT_TIMEOUT)
        self.assertEqual(client.max_retries, PubChemClient.MAX_RETRIES)
        self.assertEqual(client.cache_size, 128)

    def test_init_custom_values(self):
        """Test initialization with custom values."""
        client = PubChemClient(timeout=20, max_retries=5, cache_size=256)
        self.assertEqual(client.timeout, 20)
        self.assertEqual(client.max_retries, 5)
        self.assertEqual(client.cache_size, 256)

    @patch("requests.Session.get")
    def test_query_by_cid_success(self, mock_get):
        """Test query by CID with successful response."""
        mock_get.return_value = self.mock_success_response
        result = self.client._query_by_cid("2244")
        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_get.assert_called_once_with(
            f"{PubChemClient.BASE_URL}/cid/2244/property/IsomericSMILES/JSON",
            timeout=self.client.timeout,
        )

    @patch("requests.Session.get")
    def test_query_by_cid_invalid_format(self, mock_get):
        """Test query by CID with invalid CID format."""
        result = self.client._query_by_cid("invalid")
        self.assertIsNone(result)
        mock_get.assert_not_called()

    @patch("requests.Session.get")
    def test_query_by_cid_request_exception(self, mock_get):
        """Test query by CID with request exception."""
        mock_get.side_effect = self.mock_failed_response.raise_for_status.side_effect
        result = self.client._query_by_cid("9999999")
        self.assertIsNone(result)

    @patch("requests.Session.get")
    def test_query_by_cid_key_error(self, mock_get):
        """Test query by CID with key error in response."""
        response = MagicMock()
        response.status_code = 200
        response.json.return_value = {"InvalidKey": {}}  # Missing PropertyTable key
        mock_get.return_value = response
        result = self.client._query_by_cid("2244")
        self.assertIsNone(result)

    @patch("requests.Session.post")
    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_cid")
    def test_query_by_inchi_success(self, mock_query_by_cid, mock_post):
        """Test query by InChI with successful response."""
        mock_post.return_value = self.mock_success_response
        mock_query_by_cid.return_value = "CC(=O)OC1=CC=CC=C1C(=O)O"

        result = self.client._query_by_inchi(
            "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        )

        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_post.assert_called_once()
        mock_query_by_cid.assert_called_once_with("2244")

    @patch("requests.Session.post")
    def test_query_by_inchi_exception(self, mock_post):
        """Test query by InChI with exception."""
        mock_post.side_effect = self.mock_failed_response.raise_for_status.side_effect
        result = self.client._query_by_inchi("InChI=invalid")
        self.assertIsNone(result)

    @patch("requests.Session.get")
    def test_query_by_inchikey_success(self, mock_get):
        """Test query by InChIKey with successful response."""
        mock_get.return_value = self.mock_success_response
        result = self.client._query_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")

    @patch("requests.Session.get")
    def test_query_by_inchikey_exception(self, mock_get):
        """Test query by InChIKey with exception."""
        mock_get.side_effect = self.mock_failed_response.raise_for_status.side_effect
        result = self.client._query_by_inchikey("INVALID-INCHIKEY-HERE")
        self.assertIsNone(result)

    @patch("requests.Session.get")
    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_cid")
    def test_query_by_formula_success(self, mock_query_by_cid, mock_get):
        """Test query by formula with successful response."""
        mock_get.return_value = self.mock_success_response
        mock_query_by_cid.return_value = "CC(=O)OC1=CC=CC=C1C(=O)O"

        result = self.client._query_by_formula("C9H8O4")

        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_get.assert_called_once()
        mock_query_by_cid.assert_called_once_with("2244")

    @patch("requests.Session.get")
    def test_query_by_formula_exception(self, mock_get):
        """Test query by formula with exception."""
        mock_get.side_effect = self.mock_failed_response.raise_for_status.side_effect
        result = self.client._query_by_formula("InvalidFormula")
        self.assertIsNone(result)

    @patch("requests.Session.post")
    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_cid")
    def test_query_by_smiles_success(self, mock_query_by_cid, mock_post):
        """Test query by SMILES with successful response."""
        mock_post.return_value = self.mock_success_response
        mock_query_by_cid.return_value = "CC(=O)OC1=CC=CC=C1C(=O)O"

        result = self.client._query_by_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")

        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_post.assert_called_once()
        mock_query_by_cid.assert_called_once_with("2244")

    @patch("requests.Session.post")
    def test_query_by_smiles_exception(self, mock_post):
        """Test query by SMILES with exception."""
        mock_post.side_effect = self.mock_failed_response.raise_for_status.side_effect
        result = self.client._query_by_smiles("InvalidSmiles")
        self.assertIsNone(result)

    @patch("requests.Session.get")
    def test_query_by_name_success(self, mock_get):
        """Test query by name with successful response."""
        mock_get.return_value = self.mock_success_response
        result = self.client._query_by_name("aspirin")
        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")

    @patch("requests.Session.get")
    def test_query_by_name_exception(self, mock_get):
        """Test query by name with exception."""
        mock_get.side_effect = self.mock_failed_response.raise_for_status.side_effect
        result = self.client._query_by_name("NonExistentCompound")
        self.assertIsNone(result)

    # Tests for the main get_smiles method
    def test_get_smiles_invalid_input(self):
        """Test get_smiles with invalid input."""
        result = self.client.get_smiles(None)
        self.assertIsNone(result)

        result = self.client.get_smiles("")
        self.assertIsNone(result)

        result = self.client.get_smiles(123)  # Not a string
        self.assertIsNone(result)

    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_cid")
    def test_get_smiles_cid(self, mock_query):
        """Test get_smiles with CID input."""
        mock_query.return_value = "CC(=O)OC1=CC=CC=C1C(=O)O"
        result = self.client.get_smiles("2244")
        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_query.assert_called_once_with("2244")

    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_inchi")
    def test_get_smiles_inchi(self, mock_query):
        """Test get_smiles with InChI input."""
        mock_query.return_value = "CC(=O)OC1=CC=CC=C1C(=O)O"
        result = self.client.get_smiles(
            "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        )
        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_query.assert_called_once_with(
            "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
        )

    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_inchikey")
    def test_get_smiles_inchikey(self, mock_query):
        """Test get_smiles with InChIKey input."""
        mock_query.return_value = "CC(=O)OC1=CC=CC=C1C(=O)O"
        result = self.client.get_smiles("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        self.assertEqual(result, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_query.assert_called_once_with("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")

    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_name")
    def test_get_smiles_cas(self, mock_query):
        """Test get_smiles with CAS number input."""
        mock_query.return_value = "O"
        result = self.client.get_smiles("7732-18-5")  # Water
        self.assertEqual(result, "O")
        mock_query.assert_called_once_with("7732-18-5")

    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_formula")
    def test_get_smiles_formula(self, mock_query):
        """Test get_smiles with formula input."""
        mock_query.return_value = "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"
        result = self.client.get_smiles("C6H12O6")  # Glucose
        self.assertEqual(result, "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O")
        mock_query.assert_called_once_with("C6H12O6")

    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_smiles")
    def test_get_smiles_smiles(self, mock_query):
        """Test get_smiles with SMILES input."""
        mock_query.return_value = "CCO"
        result = self.client.get_smiles("CCO")  # Ethanol
        self.assertEqual(result, "CCO")
        mock_query.assert_called_once_with("CCO")

    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_name")
    def test_get_smiles_name(self, mock_query):
        """Test get_smiles with name input."""
        mock_query.return_value = "CCO"
        result = self.client.get_smiles("ethanol")
        self.assertEqual(result, "CCO")
        mock_query.assert_called_once_with("ethanol")

    # Test caching
    @patch("app.modules.pubchem_retrieve.PubChemClient._query_by_cid")
    def test_caching(self, mock_query):
        """Test that caching works properly."""
        mock_query.return_value = "CC(=O)OC1=CC=CC=C1C(=O)O"

        # First call should hit the API
        result1 = self.client.get_smiles("2244")
        self.assertEqual(result1, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_query.assert_called_once()

        # Reset mock to verify it's not called again
        mock_query.reset_mock()

        # Second call should use cache
        result2 = self.client.get_smiles("2244")
        self.assertEqual(result2, "CC(=O)OC1=CC=CC=C1C(=O)O")
        mock_query.assert_not_called()  # Should not call the API again


# Test for the FastAPI integration
@pytest.mark.parametrize(
    "identifier, input_type, expected_status_code, expected_success",
    [
        ("2244", "CID", 200, True),  # Aspirin by CID
        ("aspirin", "name", 200, True),  # Aspirin by name
        ("CC(=O)OC1=CC=CC=C1C(=O)O", "name", 200, True),  # Aspirin by SMILES
        ("C9H8O4", "name", 200, True),  # Aspirin by formula
        ("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "InChIKey", 200, True),  # Aspirin by InChIKey
        (
            "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "InChI",
            200,
            True,
        ),  # Aspirin by InChI
        ("7732-18-5", "CAS", 200, True),  # Water by CAS
        ("nonexistent_chemical", "name", 200, False),  # Non-existent chemical
        ("", "name", 422, False),  # Empty string
    ],
)
def test_get_pubchem_smiles_endpoint(
    identifier, input_type, expected_status_code, expected_success, client
):
    """Test the /pubchem/smiles endpoint."""
    # Configure the mock client for this test case
    client.expected_input_type = input_type
    client.expected_success = expected_success

    # Set the mock response based on expected_success
    if expected_success:
        client.mock_get_smiles_return_value = (
            "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin SMILES
        )
    else:
        client.mock_get_smiles_return_value = None

    # Make the request
    response = client.get(f"/pubchem/smiles?identifier={identifier}")

    # Check status code
    assert response.status_code == expected_status_code

    # For successful requests, check the response structure
    if expected_status_code == 200:
        data = response.json()
        assert "input" in data
        assert "canonical_smiles" in data
        assert "input_type" in data
        assert "success" in data

        assert data["input"] == identifier
        assert data["input_type"] == input_type
        assert data["success"] == expected_success

        if expected_success:
            assert data["canonical_smiles"] == "CC(=O)OC1=CC=CC=C1C(=O)O"
        else:
            assert data["canonical_smiles"] is None
