from __future__ import annotations

import time

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


@pytest.fixture
def caffeine_smiles():
    """Standard test SMILES string for caffeine."""
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def set_auth_token(monkeypatch):
    """Set a test auth token in environment."""
    test_token = "test_internal_token_12345"
    monkeypatch.setenv("CMS_INTERNAL_AUTH_TOKEN", test_token)
    # Need to reload the limiter module to pick up the new env var
    from app import limiter as limiter_module
    import importlib

    importlib.reload(limiter_module)
    yield test_token
    # Reload again to reset
    importlib.reload(limiter_module)


@pytest.fixture(autouse=True)
def reset_rate_limits():
    """Reset rate limits before each test to ensure test isolation."""
    from app.limiter import limiter

    # Reset the rate limiter state between tests
    # SlowAPI/limits stores state in the storage backend
    if hasattr(limiter._limiter, "reset"):
        limiter._limiter.reset()
    # Alternative: clear the storage backend directly
    if hasattr(limiter._limiter, "storage"):
        limiter._limiter.storage.reset()
    # Small delay to ensure storage is cleared
    time.sleep(0.1)
    yield
    # Clean up after test
    if hasattr(limiter._limiter, "reset"):
        limiter._limiter.reset()
    if hasattr(limiter._limiter, "storage"):
        limiter._limiter.storage.reset()
    time.sleep(0.1)


class TestRateLimitBasics:
    """Test basic rate limiting functionality."""

    def test_rate_limit_applied_to_mol2d_endpoint(self, caffeine_smiles):
        """Test that rate limit is applied to /convert/mol2D endpoint."""
        # Make multiple requests to hit the rate limit (20/minute)
        responses = []
        for i in range(22):  # Exceeds the 20/minute limit
            response = client.get(
                f"/latest/convert/mol2D?smiles={caffeine_smiles}&toolkit=cdk"
            )
            responses.append(response)
            # Small delay to avoid overwhelming the test
            if response.status_code == 429:
                break

        # At least one request should be rate limited
        status_codes = [r.status_code for r in responses]
        assert (
            429 in status_codes
        ), "Expected at least one 429 Too Many Requests response"

        # The last response should be 429
        assert responses[-1].status_code == 429

    def test_rate_limit_applied_to_mol3d_endpoint(self, caffeine_smiles):
        """Test that rate limit is applied to /convert/mol3D endpoint."""
        # Make multiple requests to hit the rate limit (20/minute)
        responses = []
        for i in range(22):  # Exceeds the 20/minute limit
            response = client.get(
                f"/latest/convert/mol3D?smiles={caffeine_smiles}&toolkit=rdkit"
            )
            responses.append(response)
            if response.status_code == 429:
                break

        # At least one request should be rate limited
        status_codes = [r.status_code for r in responses]
        assert (
            429 in status_codes
        ), "Expected at least one 429 Too Many Requests response"

    def test_rate_limit_applied_to_iupac_smiles_endpoint(self):
        """Test that rate limit is applied to /convert/smiles endpoint."""
        # Make multiple requests to hit the rate limit (10/minute)
        responses = []
        for i in range(12):  # Exceeds the 10/minute limit
            response = client.get(
                "/latest/convert/smiles?input_text=benzene&representation=iupac&converter=opsin"
            )
            responses.append(response)
            if response.status_code == 429:
                break

        # At least one request should be rate limited
        status_codes = [r.status_code for r in responses]
        assert (
            429 in status_codes
        ), "Expected at least one 429 Too Many Requests response"

    def test_rate_limit_applied_to_depict_3d_endpoint(self, caffeine_smiles):
        """Test that rate limit is applied to /depict/3D endpoint (25/minute)."""
        # Make multiple requests to hit the rate limit
        responses = []
        for i in range(27):  # Exceeds the 25/minute limit
            response = client.get(
                f"/latest/depict/3D?smiles={caffeine_smiles}&toolkit=rdkit"
            )
            responses.append(response)
            if response.status_code == 429:
                break

        # At least one request should be rate limited
        status_codes = [r.status_code for r in responses]
        assert (
            429 in status_codes
        ), "Expected at least one 429 Too Many Requests response"

    def test_rate_limit_applied_to_batch_convert_endpoint(self, caffeine_smiles):
        """Test that rate limit is applied to /convert/batch endpoint (10/minute)."""
        # Make multiple requests to hit the rate limit
        batch_data = {"inputs": [{"value": caffeine_smiles, "input_format": "smiles"}]}
        responses = []
        for i in range(12):  # Exceeds the 10/minute limit
            response = client.post(
                "/latest/convert/batch?output_format=canonicalsmiles&toolkit=rdkit",
                json=batch_data,
            )
            responses.append(response)
            if response.status_code == 429:
                break

        # At least one request should be rate limited
        status_codes = [r.status_code for r in responses]
        assert (
            429 in status_codes
        ), "Expected at least one 429 Too Many Requests response"


@pytest.mark.parametrize(
    "endpoint,test_input,expected_limit,description,method",
    [
        ("/latest/convert/mol2D?smiles={input}&toolkit=cdk", "CCO", 20, "mol2D", "GET"),
        (
            "/latest/convert/mol3D?smiles={input}&toolkit=rdkit",
            "CCO",
            20,
            "mol3D",
            "GET",
        ),
        (
            "/latest/depict/3D?smiles={input}&toolkit=rdkit",
            "CCO",
            25,
            "depict 3D",
            "GET",
        ),
        (
            "/latest/convert/smiles?input_text={input}&representation=iupac&converter=opsin",
            "benzene",
            10,
            "IUPAC conversion",
            "GET",
        ),
        (
            "/latest/convert/batch?output_format=canonicalsmiles&toolkit=rdkit",
            "CCO",
            10,
            "batch convert",
            "POST",
        ),
    ],
)
def test_rate_limits_parametrized(
    endpoint, test_input, expected_limit, description, method
):
    """Parametrized test for verifying rate limit configuration on different endpoints.

    This test verifies that each endpoint enforces its configured rate limit correctly:
    - mol2D and mol3D: 20 requests per minute
    - depict 3D: 25 requests per minute
    - IUPAC conversion and batch convert: 10 requests per minute
    """
    formatted_endpoint = endpoint.format(input=test_input)
    batch_data = {"inputs": [{"value": test_input, "input_format": "smiles"}]}
    # Make requests up to the limit - should all succeed
    for i in range(expected_limit):
        if method == "POST":
            response = client.post(formatted_endpoint, json=batch_data)
        else:
            response = client.get(formatted_endpoint)
        assert (
            response.status_code == 200
        ), f"{description}: Request {i+1}/{expected_limit} failed"

    # Next request should be rate limited
    if method == "POST":
        response = client.post(formatted_endpoint, json=batch_data)
    else:
        response = client.get(formatted_endpoint)
    assert (
        response.status_code == 429
    ), f"{description}: Expected rate limit after {expected_limit} requests"


class TestRateLimitNotApplied:
    """Test that endpoints without rate limiting work normally."""

    def test_no_rate_limit_on_health_endpoint(self):
        """Test that health endpoint has no rate limit."""
        # Make many requests - should all succeed
        for i in range(50):
            response = client.get("/health")
            assert response.status_code == 200
            assert response.json() == {"status": "OK"}

    def test_no_rate_limit_on_chem_descriptors_endpoint(self, caffeine_smiles):
        """Test that /chem/descriptors endpoint has no rate limit."""
        # Make many requests - should all succeed
        for i in range(30):
            response = client.get(
                f"/latest/chem/descriptors?smiles={caffeine_smiles}&toolkit=rdkit"
            )
            assert response.status_code == 200

    def test_no_rate_limit_on_chem_stereoisomers_endpoint(self, caffeine_smiles):
        """Test that /chem/stereoisomers endpoint has no rate limit."""
        # Make many requests - should all succeed
        for i in range(30):
            response = client.get(
                f"/latest/chem/stereoisomers?smiles={caffeine_smiles}"
            )
            assert response.status_code == 200

    def test_no_rate_limit_on_convert_inchi_endpoint(self, caffeine_smiles):
        """Test that /convert/inchi endpoint has no rate limit."""
        # Make many requests - should all succeed
        for i in range(30):
            response = client.get(
                f"/latest/convert/inchi?smiles={caffeine_smiles}&toolkit=rdkit"
            )
            assert response.status_code == 200


class TestAuthTokenBypass:
    """Test authentication token bypass functionality."""

    def test_auth_token_bypasses_rate_limit_mol2d(
        self, caffeine_smiles, set_auth_token
    ):
        """Test that valid auth token bypasses rate limit on mol2D endpoint."""
        # Make more requests than the rate limit allows
        headers = {"x-internal-auth": set_auth_token}

        for i in range(25):  # More than the 20/minute limit
            response = client.get(
                f"/latest/convert/mol2D?smiles={caffeine_smiles}&toolkit=cdk",
                headers=headers,
            )
            # All requests should succeed with valid auth token
            assert (
                response.status_code == 200
            ), f"Request {i+1} failed with status {response.status_code}"

    def test_auth_token_bypasses_rate_limit_mol3d(
        self, caffeine_smiles, set_auth_token
    ):
        """Test that valid auth token bypasses rate limit on mol3D endpoint."""
        headers = {"x-internal-auth": set_auth_token}

        for i in range(25):  # More than the 20/minute limit
            response = client.get(
                f"/latest/convert/mol3D?smiles={caffeine_smiles}&toolkit=rdkit",
                headers=headers,
            )
            # All requests should succeed with valid auth token
            assert (
                response.status_code == 200
            ), f"Request {i+1} failed with status {response.status_code}"

    def test_auth_token_bypasses_rate_limit_depict_3d(
        self, caffeine_smiles, set_auth_token
    ):
        """Test that valid auth token bypasses rate limit on depict/3D endpoint."""
        headers = {"x-internal-auth": set_auth_token}

        for i in range(30):  # More than the 25/minute limit
            response = client.get(
                f"/latest/depict/3D?smiles={caffeine_smiles}&toolkit=rdkit",
                headers=headers,
            )
            # All requests should succeed with valid auth token
            assert (
                response.status_code == 200
            ), f"Request {i+1} failed with status {response.status_code}"

    def test_auth_token_bypasses_rate_limit_batch_convert(
        self, caffeine_smiles, set_auth_token
    ):
        """Test that valid auth token bypasses rate limit on batch convert endpoint."""
        headers = {"x-internal-auth": set_auth_token}
        batch_data = {"inputs": [{"value": caffeine_smiles, "input_format": "smiles"}]}

        for i in range(15):  # More than the 10/minute limit
            response = client.post(
                "/latest/convert/batch?output_format=canonicalsmiles&toolkit=rdkit",
                json=batch_data,
                headers=headers,
            )
            # All requests should succeed with valid auth token
            assert (
                response.status_code == 200
            ), f"Request {i+1} failed with status {response.status_code}"

    def test_invalid_auth_token_does_not_bypass_rate_limit(
        self, caffeine_smiles, set_auth_token
    ):
        """Test that invalid auth token does not bypass rate limit."""
        # Use wrong token
        headers = {"x-internal-auth": "wrong_token_12345"}

        responses = []
        for i in range(22):  # More than the 20/minute limit
            response = client.get(
                f"/latest/convert/mol2D?smiles={caffeine_smiles}&toolkit=cdk",
                headers=headers,
            )
            responses.append(response)
            if response.status_code == 429:
                break

        # Should still be rate limited with invalid token
        status_codes = [r.status_code for r in responses]
        assert 429 in status_codes, "Invalid token should not bypass rate limit"

    def test_missing_auth_token_applies_rate_limit(self, caffeine_smiles):
        """Test that missing auth token applies normal rate limiting."""
        responses = []
        for i in range(22):  # More than the 20/minute limit
            response = client.get(
                f"/latest/convert/mol2D?smiles={caffeine_smiles}&toolkit=cdk"
                # No headers
            )
            responses.append(response)
            if response.status_code == 429:
                break

        # Should be rate limited without token
        status_codes = [r.status_code for r in responses]
        assert 429 in status_codes, "Missing token should apply rate limit"


class TestRateLimitKeyFunction:
    """Test the custom rate limit key function logic."""

    def test_rate_limit_key_function_with_valid_token(self, set_auth_token):
        """Test that custom key function returns unique keys for authenticated requests."""
        from app.limiter import custom_rate_limit_key_func
        from starlette.requests import Request

        # Mock request with valid auth token
        scope = {
            "type": "http",
            "method": "GET",
            "headers": [(b"x-internal-auth", set_auth_token.encode())],
        }
        request = Request(scope)

        key1 = custom_rate_limit_key_func(request)

        # Key should start with "auth_bypass_"
        assert key1.startswith(
            "auth_bypass_"
        ), f"Expected key to start with 'auth_bypass_', got {key1}"

    def test_rate_limit_key_function_without_token(self):
        """Test that custom key function returns IP address for unauthenticated requests."""
        from app.limiter import custom_rate_limit_key_func
        from starlette.requests import Request

        # Mock request without auth token
        scope = {
            "type": "http",
            "method": "GET",
            "headers": [],
            "client": ("127.0.0.1", 12345),
        }
        request = Request(scope)

        key = custom_rate_limit_key_func(request)

        # Key should be the IP address
        assert key == "127.0.0.1", f"Expected IP address as key, got {key}"

    def test_rate_limit_key_function_with_invalid_token(self, set_auth_token):
        """Test that custom key function returns IP for invalid token."""
        from app.limiter import custom_rate_limit_key_func
        from starlette.requests import Request

        # Mock request with invalid auth token
        scope = {
            "type": "http",
            "method": "GET",
            "headers": [(b"x-internal-auth", b"wrong_token")],
            "client": ("127.0.0.1", 12345),
        }
        request = Request(scope)

        key = custom_rate_limit_key_func(request)

        # Key should be the IP address (not bypassed)
        assert key == "127.0.0.1", f"Expected IP address for invalid token, got {key}"


class TestRateLimitResponse:
    """Test rate limit response format and headers."""

    def test_rate_limit_response_status_code(self, caffeine_smiles):
        """Test that rate limited requests return 429 status code."""
        # Exhaust rate limit
        for i in range(21):
            response = client.get(
                f"/latest/convert/mol2D?smiles={caffeine_smiles}&toolkit=cdk"
            )
            if response.status_code == 429:
                # Check response format
                assert response.status_code == 429
                # SlowAPI returns a detail message with rate limit info
                assert (
                    "detail" in response.text.lower() or "per" in response.text.lower()
                )
                break
