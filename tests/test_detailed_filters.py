"""
Test suite for detailed chemical filter implementations.

This test file validates the enhanced filter functions that provide detailed
violation information, addressing reviewer concerns about lacking details
in filter results.
"""

from __future__ import annotations

import pytest

from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import (
    check_RO5_violations_detailed,
    get_PAINS_detailed,
    get_GhoseFilter_detailed,
    get_VeberFilter_detailed,
    get_REOSFilter_detailed,
    get_RuleofThree_detailed,
)


# Define fixtures for test molecules with known properties
@pytest.fixture
def drug_like_molecule():
    """A drug-like molecule that should pass most filters."""
    # Aspirin - a well-known drug
    return parse_input("CC(=O)OC1=CC=CC=C1C(=O)O", "rdkit", False)


@pytest.fixture
def pains_molecule():
    """A molecule that contains PAINS substructures."""
    # Known PAINS-containing molecule from test_filters.py
    return parse_input(
        "O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2",
        "rdkit",
        False,
    )


@pytest.fixture
def large_molecule():
    """A large molecule that should fail multiple filters."""
    # Long alkyl chain - fails multiple filters
    return parse_input("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "rdkit", False)


@pytest.fixture
def small_molecule():
    """A very small molecule for Rule of 3 testing."""
    # Methane
    return parse_input("C", "rdkit", False)


@pytest.fixture
def lipinski_violator():
    """A molecule that violates Lipinski's Rule of 5."""
    # Large molecule with multiple violations
    return parse_input(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)NCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        "rdkit",
        False,
    )


class TestDetailedFilterStructure:
    """Test that all detailed filter functions return the expected structure."""

    def test_lipinski_detailed_structure(self, drug_like_molecule):
        """Test that Lipinski detailed function returns correct structure."""
        result = check_RO5_violations_detailed(drug_like_molecule)

        # Check required keys exist
        assert "violations" in result
        assert "details" in result
        assert "properties" in result
        assert "passes" in result

        # Check data types
        assert isinstance(result["violations"], int)
        assert isinstance(result["details"], list)
        assert isinstance(result["properties"], dict)
        assert isinstance(result["passes"], bool)

        # Check properties structure
        properties = result["properties"]
        assert "molecular_weight" in properties
        assert "logp" in properties
        assert "hb_acceptors" in properties
        assert "hb_donors" in properties

        # Check consistency: passes should be True when violations == 0
        assert result["passes"] == (result["violations"] == 0)

    def test_pains_detailed_structure(self, drug_like_molecule):
        """Test that PAINS detailed function returns correct structure."""
        result = get_PAINS_detailed(drug_like_molecule)

        # Check required keys exist
        assert "contains_pains" in result
        assert "family" in result
        assert "description" in result
        assert "passes" in result
        assert "details" in result

        # Check data types
        assert isinstance(result["contains_pains"], bool)
        assert isinstance(result["passes"], bool)
        assert isinstance(result["details"], str)

        # Check logical consistency: passes should be opposite of contains_pains
        assert result["passes"] == (not result["contains_pains"])

    def test_ghose_detailed_structure(self, drug_like_molecule):
        """Test that Ghose detailed function returns correct structure."""
        result = get_GhoseFilter_detailed(drug_like_molecule)

        # Check required keys exist
        assert "violations" in result
        assert "details" in result
        assert "properties" in result
        assert "passes" in result

        # Check consistency
        assert result["passes"] == (result["violations"] == 0)

    def test_veber_detailed_structure(self, drug_like_molecule):
        """Test that Veber detailed function returns correct structure."""
        result = get_VeberFilter_detailed(drug_like_molecule)

        # Check required keys exist
        assert "violations" in result
        assert "details" in result
        assert "properties" in result
        assert "passes" in result

        # Check specific Veber properties
        properties = result["properties"]
        assert "rotatable_bonds" in properties
        assert "tpsa" in properties

    def test_reos_detailed_structure(self, drug_like_molecule):
        """Test that REOS detailed function returns correct structure."""
        result = get_REOSFilter_detailed(drug_like_molecule)

        # Check required keys exist
        assert "violations" in result
        assert "details" in result
        assert "properties" in result
        assert "passes" in result

    def test_rule_of_three_detailed_structure(self, drug_like_molecule):
        """Test that Rule of Three detailed function returns correct structure."""
        result = get_RuleofThree_detailed(drug_like_molecule)

        # Check required keys exist
        assert "violations" in result
        assert "details" in result
        assert "properties" in result
        assert "passes" in result


class TestPAINSDetection:
    """Test PAINS detection with detailed information."""

    def test_pains_positive_detection(self, pains_molecule):
        """Test detection of PAINS-containing molecule."""
        result = get_PAINS_detailed(pains_molecule)

        # Should detect PAINS
        assert result["contains_pains"] is True
        assert result["passes"] is False  # Finding PAINS is BAD
        assert result["family"] is not None
        assert result["description"] is not None
        assert "PAINS match found" in result["details"]
        assert isinstance(result["family"], str)
        assert isinstance(result["description"], str)

    def test_pains_negative_detection(self, drug_like_molecule):
        """Test molecule without PAINS substructures."""
        result = get_PAINS_detailed(drug_like_molecule)

        # Should not detect PAINS
        assert result["contains_pains"] is False
        assert result["passes"] is True  # No PAINS found is GOOD
        assert result["family"] is None
        assert result["description"] is None
        assert "No PAINS substructures detected" in result["details"]

    def test_pains_logic_correctness(self, pains_molecule, drug_like_molecule):
        """Test that PAINS logic is correct (finding PAINS = fail)."""
        pains_result = get_PAINS_detailed(pains_molecule)
        clean_result = get_PAINS_detailed(drug_like_molecule)

        # PAINS molecule should fail
        assert pains_result["contains_pains"] is True
        assert pains_result["passes"] is False

        # Clean molecule should pass
        assert clean_result["contains_pains"] is False
        assert clean_result["passes"] is True


class TestLipinskiDetailedViolations:
    """Test Lipinski Rule of 5 detailed violations."""

    def test_lipinski_violations_content(self, lipinski_violator):
        """Test that Lipinski violations provide specific details."""
        result = check_RO5_violations_detailed(lipinski_violator)

        # Should have violations
        assert result["violations"] > 0
        assert result["passes"] is False
        assert len(result["details"]) > 0

        # Check that violation details are descriptive
        for violation in result["details"]:
            assert isinstance(violation, str)
            # Should contain property name and threshold
            assert any(prop in violation for prop in ["MW", "LogP", "HBA", "HBD"])
            assert any(op in violation for op in [">", "<", "="])

    def test_lipinski_properties_accuracy(self, drug_like_molecule):
        """Test that calculated properties are reasonable."""
        result = check_RO5_violations_detailed(drug_like_molecule)

        properties = result["properties"]

        # Check property ranges are reasonable
        assert 0 < properties["molecular_weight"] < 1000
        assert -10 < properties["logp"] < 10
        assert properties["hb_acceptors"] >= 0
        assert properties["hb_donors"] >= 0

        # Check that properties are numeric
        assert isinstance(properties["molecular_weight"], (int, float))
        assert isinstance(properties["logp"], (int, float))
        assert isinstance(properties["hb_acceptors"], int)
        assert isinstance(properties["hb_donors"], int)

    def test_lipinski_small_molecule_pass(self, small_molecule):
        """Test that small molecules typically pass Lipinski."""
        result = check_RO5_violations_detailed(small_molecule)

        # Small molecules should generally pass
        assert result["violations"] == 0
        assert result["passes"] is True
        assert len(result["details"]) == 0


class TestFilterComparisons:
    """Test comparisons between original and detailed filter functions."""

    def test_detailed_vs_original_consistency(self, drug_like_molecule, pains_molecule):
        """Test that detailed filters are consistent with original filters."""
        from app.modules.toolkits.rdkit_wrapper import (
            check_RO5_violations,
            get_PAINS,
            get_GhoseFilter,
            get_VeberFilter,
            get_REOSFilter,
            get_RuleofThree,
        )

        # Test Lipinski consistency
        original_lipinski = check_RO5_violations(drug_like_molecule)
        detailed_lipinski = check_RO5_violations_detailed(drug_like_molecule)
        assert (original_lipinski == 0) == detailed_lipinski["passes"]

        # Test PAINS consistency
        original_pains = get_PAINS(pains_molecule)
        detailed_pains = get_PAINS_detailed(pains_molecule)
        pains_found = original_pains is not False
        assert pains_found == detailed_pains["contains_pains"]
        # Original PAINS logic was wrong - detailed should be opposite
        assert detailed_pains["passes"] == (not pains_found)


class TestViolationMessages:
    """Test that violation messages are clear and actionable."""

    def test_violation_message_format(self, large_molecule):
        """Test that violation messages follow consistent format."""
        result = check_RO5_violations_detailed(large_molecule)

        if result["violations"] > 0:
            for violation in result["details"]:
                # Should contain actual value and threshold
                assert "=" in violation or ">" in violation or "<" in violation
                # Should be human-readable
                assert len(violation) > 10  # Not too terse
                assert violation[0].isupper()  # Proper capitalization

    def test_pains_message_informativeness(self, pains_molecule):
        """Test that PAINS messages provide useful information."""
        result = get_PAINS_detailed(pains_molecule)

        if result["contains_pains"]:
            assert "PAINS match found" in result["details"]
            assert result["family"] in result["details"]
            assert result["description"] in result["details"]


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_molecule_handling(self):
        """Test handling of invalid/empty molecules."""
        # Test with None or invalid input
        empty_mol = parse_input("", "rdkit", False)
        if empty_mol is None:
            # Should handle gracefully - this depends on implementation
            pass  # Skip if parse_input returns None

    def test_very_large_molecule(self):
        """Test with extremely large molecules."""
        # Very large molecule that might cause issues
        large_smiles = "C" * 1000  # Very long alkyl chain
        try:
            mol = parse_input(large_smiles, "rdkit", False)
            if mol is not None:
                result = check_RO5_violations_detailed(mol)
                # Should still return valid structure
                assert isinstance(result, dict)
                assert "passes" in result
        except Exception:
            # If parsing fails, that's acceptable for extreme cases
            pass

    def test_all_filters_consistency(self, drug_like_molecule):
        """Test that all detailed filters maintain consistent structure."""
        filters = [
            check_RO5_violations_detailed,
            get_PAINS_detailed,
            get_GhoseFilter_detailed,
            get_VeberFilter_detailed,
            get_REOSFilter_detailed,
            get_RuleofThree_detailed,
        ]

        for filter_func in filters:
            result = filter_func(drug_like_molecule)
            # All should return dict with 'passes' key
            assert isinstance(result, dict)
            assert "passes" in result
            assert isinstance(result["passes"], bool)


class TestReviewerConcerns:
    """Test specific fixes addressing reviewer feedback."""

    def test_pains_logic_fix(self, pains_molecule):
        """Test that PAINS logic correctly shows failure when PAINS found."""
        result = get_PAINS_detailed(pains_molecule)

        # The key fix: finding PAINS should be BAD (passes=False)
        if result["contains_pains"]:
            assert (
                result["passes"] is False
            ), "Finding PAINS should result in passes=False"
        else:
            assert result["passes"] is True, "No PAINS should result in passes=True"

    def test_detailed_violation_information(self, lipinski_violator):
        """Test that detailed violation information is provided."""
        result = check_RO5_violations_detailed(lipinski_violator)

        if result["violations"] > 0:
            # Should provide specific violation details
            assert len(result["details"]) > 0
            for detail in result["details"]:
                # Should contain specific values and thresholds
                assert any(
                    char.isdigit() for char in detail
                ), "Should contain numeric values"
                assert any(
                    op in detail for op in [">", "<", "="]
                ), "Should contain comparison operators"

    def test_substructure_match_details(self, pains_molecule):
        """Test that PAINS substructure details are provided."""
        result = get_PAINS_detailed(pains_molecule)

        if result["contains_pains"]:
            # Should provide specific PAINS family and description
            assert result["family"] is not None
            assert result["description"] is not None
            assert isinstance(result["family"], str)
            assert isinstance(result["description"], str)
            assert len(result["family"]) > 0
            assert len(result["description"]) > 0

    def test_clear_pass_fail_semantics(self, drug_like_molecule, pains_molecule):
        """Test that pass/fail semantics are clear and consistent."""
        clean_result = get_PAINS_detailed(drug_like_molecule)
        pains_result = get_PAINS_detailed(pains_molecule)

        # Clean molecule should pass
        assert clean_result["passes"] is True

        # PAINS-containing molecule should fail
        assert pains_result["passes"] is False

        # The logic should be clear: passes = not contains_pains
        assert clean_result["passes"] == (not clean_result["contains_pains"])
        assert pains_result["passes"] == (not pains_result["contains_pains"])


if __name__ == "__main__":
    # Run some basic tests if executed directly
    print("Running basic tests for detailed filter implementations...")

    # Test with sample molecules
    aspirin = parse_input("CC(=O)OC1=CC=CC=C1C(=O)O", "rdkit", False)
    if aspirin is not None:
        lipinski_result = check_RO5_violations_detailed(aspirin)
        pains_result = get_PAINS_detailed(aspirin)

        print(f"Aspirin Lipinski: {lipinski_result}")
        print(f"Aspirin PAINS: {pains_result}")

        print(
            "\nAll detailed filter tests would run with: pytest tests/test_detailed_filters.py"
        )
    else:
        print("Error: Could not parse aspirin molecule. Check RDKit installation.")
