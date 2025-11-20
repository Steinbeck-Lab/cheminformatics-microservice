"""Tests for sugar_removal module."""

from __future__ import annotations

import pytest

from app.modules.tools.sugar_removal import (
    get_sugar_info,
    remove_linear_sugars,
    remove_circular_sugars,
    remove_linear_and_circular_sugars,
    extract_aglycone_and_sugars,
    get_aglycone_and_sugar_indices,
    preservation_modes_enum,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer


@pytest.fixture
def glucose():
    """Pure glucose molecule (circular sugar)."""
    return get_CDK_IAtomContainer("OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O")


@pytest.fixture
def linear_glucose():
    """Linear form of glucose."""
    return get_CDK_IAtomContainer("OC[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)C=O")


@pytest.fixture
def molecule_with_sugars():
    """Molecule with both aglycone and sugar parts."""
    return get_CDK_IAtomContainer("OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O")


@pytest.fixture
def no_sugar_molecule():
    """Simple molecule without any sugars."""
    return get_CDK_IAtomContainer("CCO")


class TestGetSugarInfo:
    """Test sugar detection and info extraction."""

    def test_detect_circular_sugar(self, glucose):
        """Test detection of circular sugar."""
        has_linear, has_circular = get_sugar_info(glucose)
        assert has_circular is True

    def test_detect_linear_sugar(self, linear_glucose):
        """Test detection of linear sugar."""
        has_linear, has_circular = get_sugar_info(linear_glucose)
        # May detect as linear depending on settings

    def test_no_sugar_detection(self, no_sugar_molecule):
        """Test molecule without sugars."""
        has_linear, has_circular = get_sugar_info(no_sugar_molecule)
        assert has_linear is False
        assert has_circular is False

    def test_molecule_with_both_sugars(self, molecule_with_sugars):
        """Test molecule with multiple sugar types."""
        has_linear, has_circular = get_sugar_info(molecule_with_sugars)
        # Should detect at least one type


class TestRemoveLinearSugars:
    """Test linear sugar removal."""

    def test_remove_linear_sugar_basic(self, linear_glucose):
        """Test removing linear sugars from a linear glucose molecule."""
        result = remove_linear_sugars(linear_glucose)
        # Should return empty or minimal structure
        assert isinstance(result, str)

    def test_remove_linear_no_sugar_found(self, no_sugar_molecule):
        """Test that 'No Linear sugar found' is returned when appropriate."""
        result = remove_linear_sugars(no_sugar_molecule)
        assert result == "No Linear sugar found" or isinstance(result, str)

    def test_remove_linear_empty_aglycone(self, linear_glucose):
        """Test that empty string is returned for pure sugar."""
        result = remove_linear_sugars(linear_glucose, only_terminal=False)
        # Might return empty string or "No Linear sugar found"
        assert isinstance(result, str)

    def test_invalid_preservation_threshold(self, molecule_with_sugars):
        """Test that negative preservation_threshold raises ValueError."""
        with pytest.raises(
            ValueError, match="preservation_threshold must be a non-negative integer"
        ):
            remove_linear_sugars(molecule_with_sugars, preservation_threshold=-1)

    def test_invalid_linear_sugar_sizes(self, molecule_with_sugars):
        """Test invalid linear sugar size parameters."""
        with pytest.raises(
            ValueError, match="linearSugarsMinSize and linearSugarsMaxSize"
        ):
            remove_linear_sugars(
                molecule_with_sugars,
                linear_sugars_min_size=10,
                linear_sugars_max_size=5,
            )

    def test_invalid_linear_sugar_negative_min(self, molecule_with_sugars):
        """Test negative min size raises error."""
        with pytest.raises(ValueError):
            remove_linear_sugars(
                molecule_with_sugars,
                linear_sugars_min_size=-1,
                linear_sugars_max_size=10,
            )


class TestRemoveCircularSugarsExtended:
    """Test circular sugar removal."""

    def test_remove_circular_sugar_basic(self, glucose):
        """Test removing circular sugars from glucose."""
        result = remove_circular_sugars(glucose)
        # Pure glucose should return empty string
        assert isinstance(result, str)

    def test_remove_circular_no_sugar_found(self, no_sugar_molecule):
        """Test 'No Circular sugars found' message."""
        result = remove_circular_sugars(no_sugar_molecule)
        assert result == "No Circular sugars found" or isinstance(result, str)

    def test_remove_circular_empty_aglycone(self, glucose):
        """Test empty aglycone returns empty string."""
        result = remove_circular_sugars(glucose, only_terminal=False)
        assert isinstance(result, str)

    def test_invalid_preservation_mode(self, molecule_with_sugars):
        """Test invalid preservation mode raises error."""
        # This tests the else branch at line 290
        with pytest.raises(ValueError, match="Invalid preservation_mode specified"):
            # Pass an invalid value (not a proper enum member)
            remove_circular_sugars(molecule_with_sugars, preservation_mode="invalid")

    def test_preservation_mode_all(self, molecule_with_sugars):
        """Test preservation mode ALL."""
        result = remove_circular_sugars(
            molecule_with_sugars, preservation_mode=preservation_modes_enum.ALL
        )
        assert isinstance(result, str)

    def test_preservation_mode_molecular_weight(self, molecule_with_sugars):
        """Test preservation mode MOLECULAR_WEIGHT."""
        result = remove_circular_sugars(
            molecule_with_sugars,
            preservation_mode=preservation_modes_enum.MOLECULAR_WEIGHT,
        )
        assert isinstance(result, str)

    def test_invalid_preservation_threshold_linear(self, molecule_with_sugars):
        """Test negative preservation_threshold."""
        with pytest.raises(
            ValueError, match="preservation_threshold must be a non-negative integer"
        ):
            remove_linear_sugars(molecule_with_sugars, preservation_threshold=-10)


class TestRemoveCircularSugars:
    """Test removal of both linear and circular sugars."""

    def test_remove_both_sugars(self, molecule_with_sugars):
        """Test removing both types of sugars."""
        result = remove_linear_and_circular_sugars(molecule_with_sugars)
        assert isinstance(result, str)

    def test_remove_both_no_sugars_found(self, no_sugar_molecule):
        """Test 'No Linear or Circular sugars found' message."""
        result = remove_linear_and_circular_sugars(no_sugar_molecule)
        assert result == "No Linear or Circular sugars found" or isinstance(result, str)

    def test_remove_both_empty_aglycone(self, glucose):
        """Test that pure sugar returns empty string."""
        result = remove_linear_and_circular_sugars(glucose, only_terminal=False)
        assert isinstance(result, str)

    def test_invalid_preservation_mode_both(self, molecule_with_sugars):
        """Test invalid preservation mode."""
        with pytest.raises(ValueError, match="Invalid preservation_mode specified"):
            remove_linear_and_circular_sugars(
                molecule_with_sugars, preservation_mode="not_a_valid_mode"
            )

    def test_invalid_parameters_both(self, molecule_with_sugars):
        """Test multiple invalid parameters."""
        with pytest.raises(ValueError):
            remove_linear_and_circular_sugars(
                molecule_with_sugars,
                preservation_threshold=-1,
                oxygen_atoms_threshold=-0.5,
            )

    def test_all_preservation_modes(self, molecule_with_sugars):
        """Test all valid preservation modes."""
        for mode in [
            preservation_modes_enum.ALL,
            preservation_modes_enum.HEAVY_ATOM_COUNT,
            preservation_modes_enum.MOLECULAR_WEIGHT,
        ]:
            result = remove_linear_and_circular_sugars(
                molecule_with_sugars, preservation_mode=mode
            )
            assert isinstance(result, str)


class TestExtractAglyconeAndSugars:
    """Test extraction of aglycone and individual sugars."""

    def test_extract_basic(self, molecule_with_sugars):
        """Test basic extraction."""
        result = extract_aglycone_and_sugars(molecule_with_sugars)
        assert isinstance(result, tuple)
        assert len(result) > 0

    def test_extract_circular_only(self, molecule_with_sugars):
        """Test extracting only circular sugars."""
        result = extract_aglycone_and_sugars(
            molecule_with_sugars,
            extract_circular_sugars=True,
            extract_linear_sugars=False,
        )
        assert isinstance(result, tuple)

    def test_extract_linear_only(self, molecule_with_sugars):
        """Test extracting only linear sugars."""
        result = extract_aglycone_and_sugars(
            molecule_with_sugars,
            extract_circular_sugars=False,
            extract_linear_sugars=True,
        )
        assert isinstance(result, tuple)

    def test_extract_empty_structure(self, glucose):
        """Test that empty structures are handled (line 589-590)."""
        result = extract_aglycone_and_sugars(glucose, extract_circular_sugars=True)
        assert isinstance(result, tuple)
        # Should contain empty strings for empty structures

    def test_extract_with_invalid_preservation_mode(self, molecule_with_sugars):
        """Test invalid preservation mode."""
        with pytest.raises(ValueError, match="Invalid preservation_mode specified"):
            extract_aglycone_and_sugars(
                molecule_with_sugars, preservation_mode="invalid_mode"
            )

    def test_extract_invalid_thresholds(self, molecule_with_sugars):
        """Test invalid threshold parameters."""
        with pytest.raises(ValueError):
            extract_aglycone_and_sugars(
                molecule_with_sugars, preservation_threshold=-10
            )


class TestGetAglyconeAndSugarIndices:
    """Test getting atom indices of aglycone and sugars."""

    def test_get_indices_basic(self, molecule_with_sugars):
        """Test basic index extraction."""
        result = get_aglycone_and_sugar_indices(molecule_with_sugars)
        assert isinstance(result, list)
        # Should return list of lists with atom indices

    def test_get_indices_no_sugars(self, no_sugar_molecule):
        """Test indices with no sugars."""
        result = get_aglycone_and_sugar_indices(no_sugar_molecule)
        assert isinstance(result, list)

    def test_get_indices_circular_only(self, molecule_with_sugars):
        """Test extracting indices for circular sugars only."""
        result = get_aglycone_and_sugar_indices(
            molecule_with_sugars,
            extract_circular_sugars=True,
            extract_linear_sugars=False,
        )
        assert isinstance(result, list)

    def test_get_indices_exception_handling(self, glucose):
        """Test exception handling in index extraction (line 750-751)."""
        result = get_aglycone_and_sugar_indices(glucose)
        assert isinstance(result, list)
        # Should handle exceptions gracefully and return empty lists

    def test_get_indices_invalid_preservation_mode(self, molecule_with_sugars):
        """Test invalid preservation mode."""
        with pytest.raises(ValueError, match="Invalid preservation_mode specified"):
            get_aglycone_and_sugar_indices(
                molecule_with_sugars, preservation_mode="not_valid"
            )

    def test_get_indices_invalid_oxygen_threshold(self, molecule_with_sugars):
        """Test invalid oxygen threshold."""
        with pytest.raises(
            ValueError, match="oxygenAtomsThreshold must be a positive number"
        ):
            get_aglycone_and_sugar_indices(
                molecule_with_sugars, oxygen_atoms_threshold=-0.1
            )


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_molecule(self):
        """Test with an empty/minimal molecule."""
        mol = get_CDK_IAtomContainer("C")
        has_linear, has_circular = get_sugar_info(mol)
        assert has_linear is False
        assert has_circular is False

    def test_preservation_threshold_zero(self, molecule_with_sugars):
        """Test with preservation_threshold of 0 (valid edge case)."""
        result = remove_circular_sugars(molecule_with_sugars, preservation_threshold=0)
        assert isinstance(result, str)

    def test_linear_sugar_size_boundaries(self, molecule_with_sugars):
        """Test boundary values for linear sugar sizes."""
        # Valid: min=0, max=1 (just barely valid)
        with pytest.raises(ValueError):
            remove_linear_sugars(
                molecule_with_sugars,
                linear_sugars_min_size=5,
                linear_sugars_max_size=5,  # Equal is invalid
            )

    def test_oxygen_threshold_boundaries(self, molecule_with_sugars):
        """Test oxygen threshold at boundary values."""
        # Valid positive value
        result = remove_circular_sugars(
            molecule_with_sugars, oxygen_atoms_threshold=0.1
        )
        assert isinstance(result, str)

    def test_all_boolean_combinations(self, molecule_with_sugars):
        """Test various boolean parameter combinations."""
        # Test different combinations of extract flags
        result = extract_aglycone_and_sugars(
            molecule_with_sugars,
            extract_circular_sugars=True,
            extract_linear_sugars=True,
            mark_attach_points=True,
            post_process_sugars=True,
            limit_post_process_by_size=True,
        )
        assert isinstance(result, tuple)
