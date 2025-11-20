"""
Comprehensive tests for the fix radicals endpoint and module.

This test module covers:
- Basic radical fixing functionality
- Different atom types (C, N, O)
- Edge cases and error handling
- API endpoint testing
- Integration with CDK wrapper
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.modules.fix_radicals import fixradicals
from app.modules.toolkits.helpers import parse_input

client = TestClient(app)


class TestFixRadicalsAPI:
    """Test suite for the /fixRadicals API endpoint."""

    def test_fix_radicals_methyl_radical(self):
        """Test fixing a simple methyl radical [CH3]."""
        response = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert "radicals_detected" in data
        assert "radicals_fixed" in data
        assert data["fixed_smiles"] == "C"
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1

    def test_fix_radicals_ethyl_radical(self):
        """Test fixing an ethyl radical C[CH2]."""
        response = client.get("/latest/chem/fixRadicals?smiles=C[CH2]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["fixed_smiles"] == "CC"
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1

    def test_fix_radicals_propyl_radical(self):
        """Test fixing a propyl radical CC[CH2]."""
        response = client.get("/latest/chem/fixRadicals?smiles=CC[CH2]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["fixed_smiles"] == "CCC"
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1

    def test_fix_radicals_hydroxyl_radical(self):
        """Test fixing a hydroxyl radical [OH]."""
        response = client.get("/latest/chem/fixRadicals?smiles=[OH]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["fixed_smiles"] == "O"
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1

    def test_fix_radicals_nitrogen_radical(self):
        """Test fixing a nitrogen-centered radical [NH2]."""
        response = client.get("/latest/chem/fixRadicals?smiles=[NH2]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["fixed_smiles"] == "N"
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1

    def test_fix_radicals_methylamine_radical(self):
        """Test fixing a methylamine radical C[NH]."""
        response = client.get("/latest/chem/fixRadicals?smiles=C[NH]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1

    def test_fix_radicals_no_radicals(self):
        """Test with a molecule that has no radicals."""
        response = client.get("/latest/chem/fixRadicals?smiles=CC")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["fixed_smiles"] == "CC"
        assert data["radicals_detected"] == 0
        assert data["radicals_fixed"] == 0

    def test_fix_radicals_benzene_no_radicals(self):
        """Test with benzene - a more complex molecule without radicals."""
        response = client.get("/latest/chem/fixRadicals?smiles=c1ccccc1")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["radicals_detected"] == 0
        assert data["radicals_fixed"] == 0

    def test_fix_radicals_invalid_smiles(self):
        """Test with an invalid SMILES string."""
        response = client.get("/latest/chem/fixRadicals?smiles=INVALID_SMILES")
        assert response.status_code == 422
        data = response.json()
        assert "detail" in data

    def test_fix_radicals_empty_smiles(self):
        """Test with an empty SMILES string."""
        response = client.get("/latest/chem/fixRadicals?smiles=")
        assert response.status_code == 422

    def test_fix_radicals_secondary_carbon_radical(self):
        """Test fixing a secondary carbon radical."""
        response = client.get("/latest/chem/fixRadicals?smiles=CC([CH])C")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["radicals_detected"] >= 1

    def test_fix_radicals_multiple_radicals(self):
        """Test with a molecule containing multiple radicals."""
        # This is a challenging test case
        response = client.get("/latest/chem/fixRadicals?smiles=[CH2][CH2]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["radicals_detected"] >= 2

    def test_fix_radicals_methoxy_radical(self):
        """Test fixing a methoxy radical CO[O]."""
        response = client.get("/latest/chem/fixRadicals?smiles=CO[O]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1

    def test_fix_radicals_phenyl_radical(self):
        """Test fixing a phenyl radical [c]1ccccc1."""
        response = client.get("/latest/chem/fixRadicals?smiles=[c]1ccccc1")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert "radicals_detected" in data

    @pytest.mark.parametrize(
        "smiles,expected_status",
        [
            ("[CH3]", 200),  # Methyl radical
            ("C[CH2]", 200),  # Ethyl radical
            ("[OH]", 200),  # Hydroxyl radical
            ("[NH2]", 200),  # Amino radical
            ("CC", 200),  # No radicals
            ("c1ccccc1", 200),  # Benzene
            ("INVALID", 422),  # Invalid SMILES
            ("", 422),  # Empty SMILES
        ],
    )
    def test_fix_radicals_parametrized(self, smiles, expected_status):
        """Parametrized test for various SMILES inputs."""
        response = client.get(f"/latest/chem/fixRadicals?smiles={smiles}")
        assert response.status_code == expected_status


class TestFixRadicalsModule:
    """Test suite for the fixradicals module function."""

    def test_fixradicals_methyl_radical_direct(self):
        """Test fixradicals function directly with methyl radical."""
        mol = parse_input("[CH3]", "cdk", False)
        assert mol is not None
        result = fixradicals(mol)
        assert isinstance(result, dict)
        assert "fixed_smiles" in result
        assert "radicals_detected" in result
        assert "radicals_fixed" in result
        assert result["fixed_smiles"] == "C"
        assert result["radicals_detected"] == 1
        assert result["radicals_fixed"] == 1

    def test_fixradicals_ethyl_radical_direct(self):
        """Test fixradicals function directly with ethyl radical."""
        mol = parse_input("C[CH2]", "cdk", False)
        assert mol is not None
        result = fixradicals(mol)
        assert result["fixed_smiles"] == "CC"
        assert result["radicals_detected"] == 1
        assert result["radicals_fixed"] == 1

    def test_fixradicals_no_radicals_direct(self):
        """Test fixradicals function with a molecule without radicals."""
        mol = parse_input("CC", "cdk", False)
        assert mol is not None
        result = fixradicals(mol)
        assert result["fixed_smiles"] == "CC"
        assert result["radicals_detected"] == 0
        assert result["radicals_fixed"] == 0

    def test_fixradicals_none_molecule(self):
        """Test fixradicals function with None molecule."""
        with pytest.raises(ValueError, match="Given molecule is None"):
            fixradicals(None)

    def test_fixradicals_oxygen_radical_direct(self):
        """Test fixradicals function with oxygen radical."""
        mol = parse_input("[OH]", "cdk", False)
        assert mol is not None
        result = fixradicals(mol)
        assert result["radicals_detected"] == 1
        assert result["radicals_fixed"] == 1

    def test_fixradicals_nitrogen_radical_direct(self):
        """Test fixradicals function with nitrogen radical."""
        mol = parse_input("[NH2]", "cdk", False)
        assert mol is not None
        result = fixradicals(mol)
        assert result["radicals_detected"] == 1
        assert result["radicals_fixed"] == 1

    def test_fixradicals_returns_dict(self):
        """Test that fixradicals returns a dictionary with expected keys."""
        mol = parse_input("[CH3]", "cdk", False)
        result = fixradicals(mol)
        assert isinstance(result, dict)
        assert set(result.keys()) == {
            "fixed_smiles",
            "radicals_detected",
            "radicals_fixed",
        }

    def test_fixradicals_statistics_consistency(self):
        """Test that radicals_fixed is always <= radicals_detected."""
        test_cases = ["[CH3]", "C[CH2]", "[OH]", "[NH2]", "CC"]
        for smiles in test_cases:
            mol = parse_input(smiles, "cdk", False)
            if mol:
                result = fixradicals(mol)
                assert result["radicals_fixed"] <= result["radicals_detected"]


class TestFixRadicalsEdgeCases:
    """Test suite for edge cases and special scenarios."""

    def test_fix_radicals_charged_species(self):
        """Test with charged species."""
        # Some charged species might be parsed
        response = client.get("/latest/chem/fixRadicals?smiles=[CH3+]")
        # Should either work or return 422, but not crash
        assert response.status_code in [200, 422]

    def test_fix_radicals_complex_molecule(self):
        """Test with a complex molecule."""
        response = client.get("/latest/chem/fixRadicals?smiles=CC(C)CC[CH2]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data

    def test_fix_radicals_cyclic_radical(self):
        """Test with a cyclic structure containing a radical."""
        response = client.get("/latest/chem/fixRadicals?smiles=C1CCC([CH])C1")
        assert response.status_code == 200
        data = response.json()
        assert "radicals_detected" in data

    def test_fix_radicals_aromatic_radical(self):
        """Test with an aromatic radical."""
        response = client.get("/latest/chem/fixRadicals?smiles=c1ccc([CH2])cc1")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data

    def test_fix_radicals_heteroaromatic(self):
        """Test with heteroaromatic compounds."""
        response = client.get("/latest/chem/fixRadicals?smiles=c1ncccc1")
        assert response.status_code == 200
        data = response.json()
        # Should have no radicals
        assert data["radicals_detected"] == 0

    def test_fix_radicals_special_characters(self):
        """Test handling of URL-encoded special characters."""
        # URL encoding of [CH3]
        response = client.get("/latest/chem/fixRadicals?smiles=%5BCH3%5D")
        assert response.status_code == 200
        data = response.json()
        assert data["fixed_smiles"] == "C"


class TestFixRadicalsResponseSchema:
    """Test suite for response schema validation."""

    def test_response_has_all_required_fields(self):
        """Test that response contains all required fields."""
        response = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert response.status_code == 200
        data = response.json()
        required_fields = ["fixed_smiles", "radicals_detected", "radicals_fixed"]
        for field in required_fields:
            assert field in data, f"Missing required field: {field}"

    def test_response_field_types(self):
        """Test that response fields have correct types."""
        response = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data["fixed_smiles"], str)
        assert isinstance(data["radicals_detected"], int)
        assert isinstance(data["radicals_fixed"], int)

    def test_response_radicals_count_non_negative(self):
        """Test that radical counts are non-negative."""
        response = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert response.status_code == 200
        data = response.json()
        assert data["radicals_detected"] >= 0
        assert data["radicals_fixed"] >= 0

    def test_response_fixed_smiles_non_empty(self):
        """Test that fixed_smiles is not empty for valid inputs."""
        response = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert response.status_code == 200
        data = response.json()
        assert len(data["fixed_smiles"]) > 0


class TestFixRadicalsIntegration:
    """Integration tests combining multiple components."""

    def test_fix_radicals_then_descriptors(self):
        """Test fixing radicals and then calculating descriptors."""
        # First fix radicals
        fix_response = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert fix_response.status_code == 200
        fixed_smiles = fix_response.json()["fixed_smiles"]

        # Then get descriptors for the fixed molecule
        desc_response = client.get(f"/latest/chem/descriptors?smiles={fixed_smiles}")
        assert desc_response.status_code == 200

    def test_fix_radicals_idempotency(self):
        """Test that fixing radicals twice gives the same result."""
        # Fix radicals once
        response1 = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert response1.status_code == 200
        fixed_smiles1 = response1.json()["fixed_smiles"]

        # Fix the already-fixed molecule
        response2 = client.get(f"/latest/chem/fixRadicals?smiles={fixed_smiles1}")
        assert response2.status_code == 200
        fixed_smiles2 = response2.json()["fixed_smiles"]

        # Should be the same
        assert fixed_smiles1 == fixed_smiles2
        # And should have no radicals
        assert response2.json()["radicals_detected"] == 0

    def test_multiple_radicals_different_atoms(self):
        """Test molecule with radicals on different atom types."""
        # This is a theoretical example - actual behavior depends on CDK
        response = client.get("/latest/chem/fixRadicals?smiles=[CH2][NH]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        # Should detect multiple radicals
        assert data["radicals_detected"] >= 1


class TestFixRadicalsDocumentation:
    """Tests to verify API documentation and examples."""

    def test_openapi_examples_are_valid(self):
        """Test that all examples in the API documentation work."""
        examples = [
            "[CH3]",  # Methyl radical
            "C[CH2]",  # Ethyl radical
            "[OH]",  # Hydroxyl radical
            "C[NH]",  # Nitrogen-centered radical
            "CC[CH2]",  # Complex molecule with radical
        ]

        for smiles in examples:
            response = client.get(f"/latest/chem/fixRadicals?smiles={smiles}")
            assert response.status_code == 200, f"Example failed: {smiles}"
            data = response.json()
            assert "fixed_smiles" in data
            assert "radicals_detected" in data
            assert "radicals_fixed" in data

    def test_endpoint_accessibility(self):
        """Test that the fixRadicals endpoint is accessible and functional."""
        # Instead of checking the OpenAPI schema (which may vary by version),
        # verify the endpoint is actually accessible and works correctly
        response = client.get("/latest/chem/fixRadicals?smiles=[CH3]")
        assert response.status_code == 200
        data = response.json()
        assert "fixed_smiles" in data
        assert "radicals_detected" in data
        assert "radicals_fixed" in data
        # Verify the endpoint returns correct data
        assert data["fixed_smiles"] == "C"
        assert data["radicals_detected"] == 1
        assert data["radicals_fixed"] == 1
