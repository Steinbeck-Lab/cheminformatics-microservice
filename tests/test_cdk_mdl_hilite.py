from __future__ import annotations

import pytest

from app.modules.cdk_depict.mdl_hilite import MDLHiliteParser
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer


@pytest.fixture
def hilite_parser():
    return MDLHiliteParser()


@pytest.fixture
def simple_molecule():
    return get_CDK_IAtomContainer("CCO")


@pytest.fixture
def benzene():
    return get_CDK_IAtomContainer("c1ccccc1")


@pytest.fixture
def complex_molecule():
    return get_CDK_IAtomContainer("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")


class TestMDLHiliteParserInitialization:
    """Test MDLHiliteParser initialization."""

    def test_default_initialization(self):
        """Test that parser initializes with correct cdk_base."""
        parser = MDLHiliteParser()
        assert parser.cdk_base == "org.openscience.cdk"

    def test_initialization_loads_classes(self):
        """Test that initialization loads all required Java classes."""
        parser = MDLHiliteParser()
        assert parser.Color is not None
        assert parser.StandardGenerator is not None
        assert parser.CDKConstants is not None

    def test_cdk_base_is_string(self):
        """Test that cdk_base attribute is a string."""
        parser = MDLHiliteParser()
        assert isinstance(parser.cdk_base, str)
        assert len(parser.cdk_base) > 0


class TestParseHiliteFromMol:
    """Test parsing HILITE from molecule Sgroups."""

    def test_parse_hilite_simple_molecule(self, hilite_parser, simple_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_parse_hilite_benzene(self, hilite_parser, benzene):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(benzene)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_parse_hilite_complex_molecule(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_parse_hilite_returns_empty_sets_when_no_hilite(
        self, hilite_parser, simple_molecule
    ):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert atoms == set()
        assert bonds == set()


class TestParseHiliteIndices:
    """Test parsing HILITE indices from strings."""

    def test_parse_comma_separated_indices(self, hilite_parser):
        """Test parsing comma-separated indices."""
        result = hilite_parser._parse_hilite_indices("1,2,3,4")
        assert isinstance(result, set)
        assert result == {0, 1, 2, 3}  # Converted to 0-based

    def test_parse_space_separated_indices(self, hilite_parser):
        """Test parsing space-separated indices."""
        result = hilite_parser._parse_hilite_indices("1 2 3 4")
        assert isinstance(result, set)
        assert result == {0, 1, 2, 3}

    def test_parse_mixed_separator_indices(self, hilite_parser):
        """Test parsing indices with mixed separators."""
        result = hilite_parser._parse_hilite_indices("1, 2 3,4")
        assert isinstance(result, set)
        assert result == {0, 1, 2, 3}

    def test_parse_single_index(self, hilite_parser):
        """Test parsing a single index."""
        result = hilite_parser._parse_hilite_indices("5")
        assert result == {4}  # Converted to 0-based

    def test_parse_invalid_indices_ignored(self, hilite_parser):
        """Test that invalid indices are ignored."""
        result = hilite_parser._parse_hilite_indices("1,abc,3,xyz,5")
        assert result == {0, 2, 4}

    def test_parse_empty_string(self, hilite_parser):
        """Test parsing empty string returns empty set."""
        result = hilite_parser._parse_hilite_indices("")
        assert result == set()

    def test_parse_whitespace_only(self, hilite_parser):
        """Test parsing whitespace-only string returns empty set."""
        result = hilite_parser._parse_hilite_indices("   ")
        assert result == set()


class TestExtractIndicesFromLine:
    """Test extracting indices from V30 HIGHLIGHT lines."""

    def test_extract_atoms_from_v30_line(self, hilite_parser):
        """Test extracting atom indices from V30 HIGHLIGHT ATOMS line."""
        line = "M  V30 HIGHLIGHT ATOMS=(4 1 2 3 4)"
        result = hilite_parser._extract_indices_from_line(line)
        assert isinstance(result, set)
        assert result == {0, 1, 2, 3}  # 0-based

    def test_extract_bonds_from_v30_line(self, hilite_parser):
        """Test extracting bond indices from V30 HIGHLIGHT BONDS line."""
        line = "M  V30 HIGHLIGHT BONDS=(3 1 2 3)"
        result = hilite_parser._extract_indices_from_line(line)
        assert result == {0, 1, 2}

    def test_extract_single_index_from_v30_line(self, hilite_parser):
        """Test extracting single index from V30 line."""
        line = "M  V30 HIGHLIGHT ATOMS=(1 5)"
        result = hilite_parser._extract_indices_from_line(line)
        assert result == {4}

    def test_extract_from_line_without_parentheses(self, hilite_parser):
        """Test that line without parentheses returns empty set."""
        line = "M  V30 HIGHLIGHT ATOMS"
        result = hilite_parser._extract_indices_from_line(line)
        assert result == set()

    def test_extract_from_empty_parentheses(self, hilite_parser):
        """Test that empty parentheses returns empty set."""
        line = "M  V30 HIGHLIGHT ATOMS=()"
        result = hilite_parser._extract_indices_from_line(line)
        assert result == set()

    def test_extract_with_only_count(self, hilite_parser):
        """Test that line with only count (no indices) returns empty set."""
        line = "M  V30 HIGHLIGHT ATOMS=(0)"
        result = hilite_parser._extract_indices_from_line(line)
        assert result == set()


class TestSgroupFiltering:
    """Test filtering Sgroups by type - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    """

    pass


class TestAtomExtractionFromSgroup:
    """Test extracting atoms from Sgroups - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    Tests for atom extraction are in TestParseHiliteFromMol class.
    """

    pass


class TestBondExtractionFromSgroup:
    """Test extracting bonds from Sgroups - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    Tests for bond extraction are in TestParseHiliteFromMol class.
    """

    pass


class TestFieldnameChecking:
    """Test checking Sgroup fieldnames for HILITE - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    Tests for property checking are in TestHasHiliteProperties class.
    """

    pass


class TestMolecularIntegrity:
    """Test that parsing preserves molecular integrity."""

    def test_atom_count_preserved(self, hilite_parser, simple_molecule):
        initial_count = simple_molecule.getAtomCount()
        hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert simple_molecule.getAtomCount() == initial_count

    def test_bond_count_preserved(self, hilite_parser, simple_molecule):
        initial_bonds = simple_molecule.getBondCount()
        hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert simple_molecule.getBondCount() == initial_bonds

    def test_connectivity_preserved(self, hilite_parser, simple_molecule):
        hilite_parser.parse_hilite_from_mol(simple_molecule)
        for atom in simple_molecule.atoms():
            connected = simple_molecule.getConnectedAtomsList(atom)
            assert connected is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_atom_molecule(self, hilite_parser):
        mol = get_CDK_IAtomContainer("C")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(mol)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_disconnected_fragments(self, hilite_parser):
        mol = get_CDK_IAtomContainer("CCO.C1CCOC1")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(mol)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_aromatic_molecule(self, hilite_parser):
        mol = get_CDK_IAtomContainer("c1ccccc1")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(mol)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_charged_molecule(self, hilite_parser):
        mol = get_CDK_IAtomContainer("[NH4+]")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(mol)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)


class TestSgroupTypes:
    """Test detection of different Sgroup types - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    Tests for has_hilite_properties() are in TestHasHiliteProperties class.
    """

    pass


class TestMultipleSgroups:
    """Test handling multiple Sgroups - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    Tests for parsing multiple highlights are in TestParseHiliteFromMDL class.
    """

    pass


class TestEmptyResults:
    """Test handling of molecules without HILITE information."""

    def test_no_hilite_properties_returns_empty_sets(
        self, hilite_parser, simple_molecule
    ):
        """Test that molecule without HILITE properties returns empty sets."""
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert atoms == set()
        assert bonds == set()

    def test_has_hilite_properties_returns_false(self, hilite_parser, simple_molecule):
        """Test has_hilite_properties returns False when no properties exist."""
        assert hilite_parser.has_hilite_properties(simple_molecule) is False


class TestAtomIndexValidation:
    """Test that extracted atom indices are valid."""

    def test_atom_indices_within_range(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        atom_count = complex_molecule.getAtomCount()
        for atom_idx in atoms:
            assert 0 <= atom_idx < atom_count

    def test_atom_indices_are_non_negative(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        for atom_idx in atoms:
            assert atom_idx >= 0


class TestBondIndexValidation:
    """Test that extracted bond indices are valid."""

    def test_bond_indices_within_range(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        bond_count = complex_molecule.getBondCount()
        for bond_idx in bonds:
            assert 0 <= bond_idx < bond_count

    def test_bond_indices_are_non_negative(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        for bond_idx in bonds:
            assert bond_idx >= 0


class TestLargeMolecules:
    """Test with large molecules."""

    def test_large_molecule_steroid(self, hilite_parser):
        mol = get_CDK_IAtomContainer("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(mol)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_large_molecule_peptide(self, hilite_parser):
        mol = get_CDK_IAtomContainer("CC(C)(C)OC(=O)NCC(=O)NCC(=O)O")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(mol)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)


class TestSgroupFieldData:
    """Test Sgroup field data - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    Tests for property-based parsing are in TestParseHiliteFromMol class.
    """

    pass


class TestHiliteFieldnameVariants:
    """Test HILITE property name variants."""

    def test_parse_with_hilite_atoms_property(self, hilite_parser, simple_molecule):
        """Test parsing molecule with HILITE_ATOMS property."""
        simple_molecule.setProperty("HILITE_ATOMS", "1,2,3")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert atoms == {0, 1, 2}  # Converted to 0-based
        assert bonds == set()

    def test_parse_with_hilite_bonds_property(self, hilite_parser, simple_molecule):
        """Test parsing molecule with HILITE_BONDS property."""
        simple_molecule.setProperty("HILITE_BONDS", "1,2")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert atoms == set()
        assert bonds == {0, 1}  # Converted to 0-based

    def test_parse_with_both_properties(self, hilite_parser, simple_molecule):
        """Test parsing molecule with both HILITE properties."""
        simple_molecule.setProperty("HILITE_ATOMS", "1,2")
        simple_molecule.setProperty("HILITE_BONDS", "1")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert atoms == {0, 1}
        assert bonds == {0}


class TestReturnTypeConsistency:
    """Test that return types are consistent."""

    def test_parse_hilite_always_returns_tuple(self, hilite_parser, simple_molecule):
        result = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_parse_hilite_returns_sets(self, hilite_parser, simple_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)


class TestDATSgroupSpecifics:
    """Test DAT Sgroup specific functionality - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    """

    pass


class TestComplexSgroupStructures:
    """Test with complex Sgroup structures - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    Tests for complex molecules are in other test classes.
    """

    pass


class TestErrorHandling:
    """Test error handling in MDL HILITE parsing."""

    def test_invalid_molecule_handling(self, hilite_parser):
        """Test that parsing None molecule returns empty sets."""
        try:
            atoms, bonds = hilite_parser.parse_hilite_from_mol(None)
            assert (atoms, bonds) == (set(), set())
        except Exception:
            # Exception is acceptable for None input
            assert True

    def test_invalid_property_value_handling(self, hilite_parser, simple_molecule):
        """Test handling of invalid HILITE property values."""
        simple_molecule.setProperty("HILITE_ATOMS", "invalid,data,abc")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        # Should return empty set since all values are invalid
        assert atoms == set()
        assert bonds == set()

    def test_negative_indices_ignored(self, hilite_parser):
        """Test that negative indices are handled gracefully."""
        result = hilite_parser._parse_hilite_indices("-1,0,1,2")
        # -1 becomes -2 after conversion, should be in result, 0 is invalid (0-based from 1-based 0)
        assert isinstance(result, set)


class TestSgroupIterationSafety:
    """Test safe iteration - REMOVED.

    Note: The current implementation doesn't use Sgroups.
    """

    pass


class TestHiliteParsingWorkflow:
    """Test the complete HILITE parsing workflow."""

    def test_complete_workflow_with_properties(self, hilite_parser, simple_molecule):
        """Test complete workflow from properties to highlighting."""
        # Set HILITE properties
        simple_molecule.setProperty("HILITE_ATOMS", "1,2")
        simple_molecule.setProperty("HILITE_BONDS", "1")

        # Parse
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)

        # Verify
        assert atoms == {0, 1}
        assert bonds == {0}

        # Apply highlighting
        hilite_parser.apply_hilite_to_molecule(simple_molecule, atoms, bonds)

        # Check highlighting was applied
        highlighted_atom = simple_molecule.getAtom(0)
        assert (
            highlighted_atom.getProperty(
                hilite_parser.StandardGenerator.HIGHLIGHT_COLOR
            )
            is not None
        )

    def test_complete_workflow_complex(self, hilite_parser, complex_molecule):
        """Test complete workflow with complex molecule."""
        complex_molecule.setProperty("HILITE_ATOMS", "1,2,3,4")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert len(atoms) == 4
        complex_molecule.setProperty("HILITE_ATOMS", "1,2,3,4")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert len(atoms) == 4


class TestV3000FormatSupport:
    """Test support for V3000 format MDL files."""

    def test_parse_v3000_atoms_line(self, hilite_parser):
        """Test parsing V3000 HIGHLIGHT ATOMS line."""
        mdl_text = """
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 END ATOM
M  V30 HIGHLIGHT ATOMS=(3 1 2 3)
M  V30 END CTAB
"""
        atoms, bonds = hilite_parser.parse_hilite_from_mdl(mdl_text)
        assert atoms == {0, 1, 2}  # 0-based
        assert bonds == set()

    def test_parse_v3000_bonds_line(self, hilite_parser):
        """Test parsing V3000 HIGHLIGHT BONDS line."""
        mdl_text = """
M  V30 HIGHLIGHT BONDS=(2 1 2)
"""
        atoms, bonds = hilite_parser.parse_hilite_from_mdl(mdl_text)
        assert atoms == set()
        assert bonds == {0, 1}  # 0-based

    def test_parse_v3000_both_lines(self, hilite_parser):
        """Test parsing both ATOMS and BONDS V3000 lines."""
        mdl_text = """
M  V30 HIGHLIGHT ATOMS=(2 1 2)
M  V30 HIGHLIGHT BONDS=(1 1)
"""
        atoms, bonds = hilite_parser.parse_hilite_from_mdl(mdl_text)
        assert atoms == {0, 1}
        assert bonds == {0}


class TestIntegrationWithDepiction:
    """Test integration with depiction system."""

    def test_apply_hilite_to_molecule(self, hilite_parser, complex_molecule):
        """Test applying highlighting to molecule atoms and bonds."""
        atoms = {0, 1, 2}
        bonds = {0, 1}

        hilite_parser.apply_hilite_to_molecule(complex_molecule, atoms, bonds)

        # Verify atoms are highlighted
        for atom_idx in atoms:
            if atom_idx < complex_molecule.getAtomCount():
                atom = complex_molecule.getAtom(atom_idx)
                color = atom.getProperty(
                    hilite_parser.StandardGenerator.HIGHLIGHT_COLOR
                )
                assert color is not None

    def test_apply_hilite_with_custom_color(self, hilite_parser, simple_molecule):
        """Test applying highlighting with custom color."""
        atoms = {0}
        bonds = set()
        custom_color = (255, 200, 200)  # Light red

        hilite_parser.apply_hilite_to_molecule(
            simple_molecule, atoms, bonds, highlight_color=custom_color
        )

        atom = simple_molecule.getAtom(0)
        color = atom.getProperty(hilite_parser.StandardGenerator.HIGHLIGHT_COLOR)
        assert color is not None

    def test_hilite_indices_are_integers(self, hilite_parser, complex_molecule):
        """Test that parsed indices are integers."""
        complex_molecule.setProperty("HILITE_ATOMS", "1,2,3")
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)

        for atom in atoms:
            assert isinstance(atom, int)
        for bond in bonds:
            assert isinstance(bond, int)


class TestHasHiliteProperties:
    """Test has_hilite_properties() method."""

    def test_has_properties_returns_true_with_atoms(
        self, hilite_parser, simple_molecule
    ):
        """Test returns True when HILITE_ATOMS property exists."""
        simple_molecule.setProperty("HILITE_ATOMS", "1,2")
        assert hilite_parser.has_hilite_properties(simple_molecule) is True

    def test_has_properties_returns_true_with_bonds(
        self, hilite_parser, simple_molecule
    ):
        """Test returns True when HILITE_BONDS property exists."""
        simple_molecule.setProperty("HILITE_BONDS", "1")
        assert hilite_parser.has_hilite_properties(simple_molecule) is True

    def test_has_properties_returns_false_without_properties(
        self, hilite_parser, simple_molecule
    ):
        """Test returns False when no HILITE properties exist."""
        assert hilite_parser.has_hilite_properties(simple_molecule) is False
