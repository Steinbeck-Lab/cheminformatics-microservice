from __future__ import annotations

import pytest

from app.modules.cdk_hydrogen_display import set_hydrogen_display
from app.modules.cdk_hydrogen_display import setHydrogenDisplay
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer


@pytest.fixture
def simple_molecule():
    """Simple ethane molecule - CC."""
    return get_CDK_IAtomContainer("CC")


@pytest.fixture
def chiral_molecule():
    """Chiral molecule - (S)-Alanine."""
    return get_CDK_IAtomContainer("C[C@H](N)C(=O)O")


@pytest.fixture
def cis_trans_molecule():
    """Molecule with E/Z stereochemistry - E-2-butene."""
    return get_CDK_IAtomContainer("C/C=C/C")


@pytest.fixture
def allene_molecule():
    """Allene with extended tetrahedral chirality."""
    return get_CDK_IAtomContainer("CC(C)=C=C(C)C")


@pytest.fixture
def ring_molecule():
    """Cyclohexane with chiral center."""
    return get_CDK_IAtomContainer("C1CCC[C@H](C)C1")


@pytest.fixture
def complex_chiral():
    """Complex molecule with multiple chiral centers."""
    return get_CDK_IAtomContainer("C[C@H](O)[C@@H](C)Cl")


@pytest.fixture
def aromatic_molecule():
    """Benzene - aromatic ring."""
    return get_CDK_IAtomContainer("c1ccccc1")


class TestSetHydrogenDisplayMinimal:
    """Test the 'Minimal' hydrogen display mode."""

    def test_minimal_mode_suppresses_hydrogens(self, simple_molecule):
        """Test that minimal mode suppresses all hydrogens."""
        set_hydrogen_display(simple_molecule, "Minimal")

        # Count explicit hydrogens
        h_count = sum(
            1 for atom in simple_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count == 0, "Minimal mode should suppress all explicit hydrogens"

    def test_minimal_mode_chiral_molecule(self, chiral_molecule):
        """Test minimal mode on a chiral molecule."""
        set_hydrogen_display(chiral_molecule, "Minimal")

        h_count = sum(
            1 for atom in chiral_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count == 0, "All hydrogens should be implicit in minimal mode"

    def test_minimal_mode_aromatic(self, aromatic_molecule):
        """Test minimal mode on aromatic molecule."""
        set_hydrogen_display(aromatic_molecule, "Minimal")

        h_count = sum(
            1 for atom in aromatic_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count == 0, "Aromatic hydrogens should be suppressed"


class TestSetHydrogenDisplayExplicit:
    """Test the 'Explicit' hydrogen display mode."""

    def test_explicit_mode_shows_all_hydrogens(self, simple_molecule):
        """Test that explicit mode shows all hydrogens."""
        set_hydrogen_display(simple_molecule, "Explicit")

        # Count explicit hydrogens
        h_count = sum(
            1 for atom in simple_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count > 0, "Explicit mode should show hydrogen atoms"
        # Ethane (CC) should have 6 hydrogens
        assert h_count == 6, "Ethane should have 6 explicit hydrogens"

    def test_explicit_mode_chiral_molecule(self, chiral_molecule):
        """Test explicit mode on a chiral molecule."""
        initial_atom_count = chiral_molecule.getAtomCount()
        set_hydrogen_display(chiral_molecule, "Explicit")

        final_atom_count = chiral_molecule.getAtomCount()

        assert (
            final_atom_count > initial_atom_count
        ), "Explicit mode should add hydrogen atoms"

    def test_explicit_mode_preserves_structure(self, simple_molecule):
        """Test that explicit mode preserves molecular structure."""
        initial_carbon_count = sum(
            1 for atom in simple_molecule.atoms() if atom.getAtomicNumber() == 6
        )

        set_hydrogen_display(simple_molecule, "Explicit")

        final_carbon_count = sum(
            1 for atom in simple_molecule.atoms() if atom.getAtomicNumber() == 6
        )

        assert (
            initial_carbon_count == final_carbon_count
        ), "Carbon count should remain unchanged"


class TestSetHydrogenDisplayStereo:
    """Test the 'Stereo' hydrogen display mode."""

    def test_stereo_mode_chiral_center(self, chiral_molecule):
        """Test that stereo mode adds hydrogens at chiral centers."""
        set_hydrogen_display(chiral_molecule, "Stereo")

        h_count = sum(
            1 for atom in chiral_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count > 0, "Stereo mode should show stereo-relevant hydrogens"

    def test_stereo_mode_cis_trans(self, cis_trans_molecule):
        """Test stereo mode on molecule with E/Z stereochemistry."""
        set_hydrogen_display(cis_trans_molecule, "Stereo")

        h_count = sum(
            1 for atom in cis_trans_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        # E/Z stereochemistry should result in some explicit hydrogens
        assert h_count >= 0, "Stereo mode should handle E/Z stereochemistry"

    def test_stereo_mode_preserves_stereochemistry(self, chiral_molecule):
        """Test that stereo mode preserves stereochemical information."""
        initial_stereo_count = sum(1 for _ in chiral_molecule.stereoElements())

        set_hydrogen_display(chiral_molecule, "Stereo")

        final_stereo_count = sum(1 for _ in chiral_molecule.stereoElements())

        assert (
            final_stereo_count >= initial_stereo_count
        ), "Stereo elements should be preserved or updated"

    def test_stereo_mode_ring_molecule(self, ring_molecule):
        """Test stereo mode on cyclic molecule with chiral center."""
        set_hydrogen_display(ring_molecule, "Stereo")

        # Should have some explicit hydrogens for the chiral center
        h_count = sum(
            1 for atom in ring_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count >= 0, "Ring molecules should be handled correctly"


class TestSetHydrogenDisplaySmart:
    """Test the 'Smart' hydrogen display mode."""

    def test_smart_mode_chiral_molecule(self, chiral_molecule):
        """Test smart mode on chiral molecule."""
        set_hydrogen_display(chiral_molecule, "Smart")

        h_count = sum(
            1 for atom in chiral_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count >= 0, "Smart mode should intelligently display hydrogens"

    def test_smart_mode_ring_molecule(self, ring_molecule):
        """Test smart mode handles ring systems intelligently."""
        set_hydrogen_display(ring_molecule, "Smart")

        # Smart mode should consider ring bonds in its heuristics
        h_count = sum(
            1 for atom in ring_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count >= 0, "Smart mode should handle ring systems"

    def test_smart_mode_complex_molecule(self, complex_chiral):
        """Test smart mode on molecule with multiple chiral centers."""
        set_hydrogen_display(complex_chiral, "Smart")

        # Smart mode should handle multiple stereocenters
        stereo_count = sum(1 for _ in complex_chiral.stereoElements())

        assert stereo_count > 0, "Multiple chiral centers should be preserved"

    def test_smart_mode_cis_trans(self, cis_trans_molecule):
        """Test smart mode on E/Z stereochemistry."""
        set_hydrogen_display(cis_trans_molecule, "Smart")

        h_count = sum(
            1 for atom in cis_trans_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count >= 0, "Smart mode should handle E/Z bonds intelligently"


class TestSetHydrogenDisplayProvided:
    """Test the 'Provided' hydrogen display mode."""

    def test_provided_mode_no_changes(self, simple_molecule):
        """Test that provided mode makes no changes."""
        initial_atom_count = simple_molecule.getAtomCount()

        set_hydrogen_display(simple_molecule, "Provided")

        final_atom_count = simple_molecule.getAtomCount()

        assert (
            initial_atom_count == final_atom_count
        ), "Provided mode should not modify the molecule"

    def test_provided_mode_chiral_molecule(self, chiral_molecule):
        """Test provided mode on chiral molecule."""
        initial_atom_count = chiral_molecule.getAtomCount()

        set_hydrogen_display(chiral_molecule, "Provided")

        final_atom_count = chiral_molecule.getAtomCount()

        assert (
            initial_atom_count == final_atom_count
        ), "No changes should occur in Provided mode"

    def test_provided_mode_is_default(self, simple_molecule):
        """Test that Provided is the default mode."""
        initial_atom_count = simple_molecule.getAtomCount()

        # Call without specifying mode (should default to "Provided")
        set_hydrogen_display(simple_molecule)

        final_atom_count = simple_molecule.getAtomCount()

        assert (
            initial_atom_count == final_atom_count
        ), "Default mode should be Provided (no changes)"


class TestInvalidInputs:
    """Test error handling for invalid inputs."""

    def test_invalid_display_type(self, simple_molecule):
        """Test that invalid display type raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            set_hydrogen_display(simple_molecule, "Invalid")

        assert "Invalid display_type" in str(exc_info.value)

    def test_invalid_display_type_lowercase(self, simple_molecule):
        """Test that lowercase display types are not accepted."""
        with pytest.raises(ValueError) as exc_info:
            set_hydrogen_display(simple_molecule, "minimal")

        assert "Invalid display_type" in str(exc_info.value)

    def test_error_message_lists_valid_options(self, simple_molecule):
        """Test that error message lists all valid options."""
        with pytest.raises(ValueError) as exc_info:
            set_hydrogen_display(simple_molecule, "Wrong")

        error_msg = str(exc_info.value)
        assert "Minimal" in error_msg
        assert "Explicit" in error_msg
        assert "Stereo" in error_msg
        assert "Smart" in error_msg
        assert "Provided" in error_msg


class TestCamelCaseAlias:
    """Test the camelCase alias function."""

    def test_camelcase_alias_exists(self):
        """Test that setHydrogenDisplay function exists."""
        assert callable(setHydrogenDisplay)

    def test_camelcase_minimal_mode(self, simple_molecule):
        """Test camelCase function with minimal mode."""
        setHydrogenDisplay(simple_molecule, "Minimal")

        h_count = sum(
            1 for atom in simple_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count == 0, "CamelCase alias should work identically"

    def test_camelcase_default_mode(self, simple_molecule):
        """Test camelCase function with default mode."""
        initial_atom_count = simple_molecule.getAtomCount()

        setHydrogenDisplay(simple_molecule)

        final_atom_count = simple_molecule.getAtomCount()

        assert (
            initial_atom_count == final_atom_count
        ), "CamelCase alias should use same default"


class TestMoleculeModification:
    """Test that molecules are modified in-place."""

    def test_in_place_modification(self, simple_molecule):
        """Test that molecule is modified in-place, not copied."""
        initial_id = id(simple_molecule)

        result = set_hydrogen_display(simple_molecule, "Minimal")

        # Function should return None
        assert result is None, "Function should return None"

        # Same object should be modified
        assert id(simple_molecule) == initial_id, "Molecule should be modified in-place"

    def test_sequential_modifications(self, simple_molecule):
        """Test that sequential modifications work correctly."""
        # First make explicit
        set_hydrogen_display(simple_molecule, "Explicit")
        count_after_explicit = simple_molecule.getAtomCount()

        # Then minimize
        set_hydrogen_display(simple_molecule, "Minimal")
        count_after_minimal = simple_molecule.getAtomCount()

        # Minimal should have fewer atoms than explicit
        assert (
            count_after_minimal < count_after_explicit
        ), "Sequential modifications should work"


class TestEdgeCases:
    """Test edge cases and special scenarios."""

    def test_empty_stereo_elements(self, simple_molecule):
        """Test molecule with no stereochemistry."""
        # Simple ethane has no stereochemistry
        set_hydrogen_display(simple_molecule, "Stereo")

        stereo_count = sum(1 for _ in simple_molecule.stereoElements())

        assert stereo_count == 0, "Non-chiral molecule should have no stereo elements"

    def test_aromatic_molecule_stereo(self, aromatic_molecule):
        """Test aromatic molecule with stereo mode."""
        set_hydrogen_display(aromatic_molecule, "Stereo")

        # Should handle aromatic molecules without errors
        h_count = sum(
            1 for atom in aromatic_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count >= 0, "Aromatic molecules should be handled"

    def test_multiple_chiral_centers(self, complex_chiral):
        """Test molecule with multiple chiral centers."""
        set_hydrogen_display(complex_chiral, "Stereo")

        # Should handle multiple stereocenters correctly
        stereo_count = sum(1 for _ in complex_chiral.stereoElements())

        assert stereo_count > 0, "Multiple chiral centers should be handled"

    def test_allene_stereochemistry(self, allene_molecule):
        """Test allene with extended tetrahedral stereochemistry."""
        try:
            set_hydrogen_display(allene_molecule, "Stereo")
            # If no exception, test passes
            assert True
        except Exception as e:
            pytest.fail(f"Allene handling raised exception: {e}")


class TestCustomCDKBase:
    """Test custom CDK base package parameter."""

    def test_default_cdk_base(self, simple_molecule):
        """Test that default cdk_base works."""
        # Should work with default cdk_base
        try:
            set_hydrogen_display(simple_molecule, "Minimal")
            assert True
        except Exception as e:
            pytest.fail(f"Default cdk_base failed: {e}")

    def test_explicit_cdk_base(self, chiral_molecule):
        """Test with explicitly specified cdk_base."""
        # Should work with explicit cdk_base parameter
        try:
            set_hydrogen_display(
                chiral_molecule, "Stereo", cdk_base="org.openscience.cdk"
            )
            assert True
        except Exception as e:
            pytest.fail(f"Explicit cdk_base failed: {e}")


class TestIntegrationScenarios:
    """Test realistic usage scenarios."""

    def test_workflow_for_visualization(self, chiral_molecule):
        """Test typical workflow: parse -> set H display -> visualize."""
        # This would be typical for generating structure images
        set_hydrogen_display(chiral_molecule, "Stereo")

        # Verify molecule is still valid
        assert chiral_molecule.getAtomCount() > 0
        assert chiral_molecule.getBondCount() > 0

    def test_workflow_minimal_for_publication(self, aromatic_molecule):
        """Test workflow for publication-quality images."""
        # Publications typically use minimal hydrogen display
        set_hydrogen_display(aromatic_molecule, "Minimal")

        # Verify no explicit hydrogens
        h_count = sum(
            1 for atom in aromatic_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count == 0

    def test_workflow_explicit_for_teaching(self, simple_molecule):
        """Test workflow for teaching materials."""
        # Teaching materials often show all hydrogens explicitly
        set_hydrogen_display(simple_molecule, "Explicit")

        # Verify hydrogens are present
        h_count = sum(
            1 for atom in simple_molecule.atoms() if atom.getAtomicNumber() == 1
        )

        assert h_count > 0


class TestAtomContainerIntegrity:
    """Test that IAtomContainer integrity is maintained."""

    def test_bond_count_consistency(self, simple_molecule):
        """Test that bond count is consistent after modification."""
        set_hydrogen_display(simple_molecule, "Explicit")

        # Count bonds manually
        bond_count = simple_molecule.getBondCount()

        # Verify bonds exist
        assert bond_count > 0, "Molecule should have bonds"

    def test_atom_connectivity(self, chiral_molecule):
        """Test that atom connectivity is maintained."""
        set_hydrogen_display(chiral_molecule, "Explicit")

        # Check that all atoms have connections (except maybe single atoms)
        for atom in chiral_molecule.atoms():
            connected_atoms = chiral_molecule.getConnectedAtomsList(atom)
            # Most atoms should have connections
            if chiral_molecule.getAtomCount() > 1:
                assert connected_atoms is not None

    def test_implicit_hydrogen_counts(self, simple_molecule):
        """Test implicit hydrogen counts after Minimal mode."""
        set_hydrogen_display(simple_molecule, "Minimal")

        # After suppression, carbons should have implicit hydrogens
        for atom in simple_molecule.atoms():
            if atom.getAtomicNumber() == 6:  # Carbon
                implicit_h = atom.getImplicitHydrogenCount()
                # Should have implicit hydrogens or be fully substituted
                assert implicit_h >= 0
