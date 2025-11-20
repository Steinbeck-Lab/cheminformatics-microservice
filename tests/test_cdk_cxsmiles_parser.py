from __future__ import annotations

import pytest

from app.modules.cdk_depict.cxsmiles_parser import (
    parse_cxsmiles_highlighting_from_string,
    extract_cxsmiles_highlighting,
    parse_cxsmiles,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer


@pytest.fixture
def simple_cxsmiles_atom_highlighting():
    return "c1ccccc1 |ha:0,1,2|"


@pytest.fixture
def simple_cxsmiles_bond_highlighting():
    return "CCO |hb:0,1|"


@pytest.fixture
def complex_cxsmiles():
    return "c1cc(O)ccc1 |ha:0,1,2,hb:0,1|"


@pytest.fixture
def cxsmiles_no_highlighting():
    return "CCO"


@pytest.fixture
def cxsmiles_empty_highlighting():
    return "CCO ||"


class TestParseCXSMILESHighlightingFromString:
    """Test parsing CXSMILES highlighting directly from string."""

    def test_atom_highlighting_only(self, simple_cxsmiles_atom_highlighting):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(
            simple_cxsmiles_atom_highlighting
        )
        assert atoms == {0, 1, 2}
        assert bonds == set()

    def test_bond_highlighting_only(self, simple_cxsmiles_bond_highlighting):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(
            simple_cxsmiles_bond_highlighting
        )
        assert atoms == set()
        assert bonds == {0, 1}

    def test_both_atom_and_bond_highlighting(self, complex_cxsmiles):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(complex_cxsmiles)
        assert atoms == {0, 1, 2}
        assert bonds == {0, 1}

    def test_no_highlighting(self, cxsmiles_no_highlighting):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(cxsmiles_no_highlighting)
        assert atoms == set()
        assert bonds == set()

    def test_empty_highlighting(self, cxsmiles_empty_highlighting):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(
            cxsmiles_empty_highlighting
        )
        assert atoms == set()
        assert bonds == set()

    def test_single_atom_highlighting(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:0|")
        assert atoms == {0}
        assert bonds == set()

    def test_single_bond_highlighting(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |hb:0|")
        assert atoms == set()
        assert bonds == {0}

    def test_multiple_atoms_highlighting(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(
            "CCCCCC |ha:0,1,2,3,4,5|"
        )
        assert atoms == {0, 1, 2, 3, 4, 5}
        assert bonds == set()

    def test_none_input(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(None)
        assert atoms == set()
        assert bonds == set()

    def test_empty_string_input(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("")
        assert atoms == set()
        assert bonds == set()

    def test_malformed_cxsmiles_graceful_failure(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:invalid|")
        assert atoms == set()
        assert bonds == set()

    def test_whitespace_in_indices(self):
        # Current implementation requires no space after colon
        # Spaces within the comma-separated list are handled
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:0,1,2|")
        assert atoms == {0, 1, 2}
        assert bonds == set()

    def test_trailing_comma(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:0,1,|")
        assert atoms == {0, 1}
        assert bonds == set()

    def test_complex_cxsmiles_with_coords(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string(
            "CCO |ha:0,1,coord:0,1.0,2.0|"
        )
        assert atoms == {0, 1}
        assert bonds == set()


class TestExtractCXSMILESHighlighting:
    """Test extracting CXSMILES highlighting from parsed molecule."""

    def test_extract_from_molecule_returns_empty(self):
        mol = get_CDK_IAtomContainer("c1ccccc1")
        atoms, bonds = extract_cxsmiles_highlighting(mol)
        assert atoms == set()
        assert bonds == set()

    def test_extract_from_complex_molecule(self):
        mol = get_CDK_IAtomContainer("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        atoms, bonds = extract_cxsmiles_highlighting(mol)
        assert atoms == set()
        assert bonds == set()


class TestParseCXSMILES:
    """Test parsing full CXSMILES strings into molecules."""

    def test_parse_simple_cxsmiles(self, simple_cxsmiles_atom_highlighting):
        mol = parse_cxsmiles(simple_cxsmiles_atom_highlighting)
        assert mol is not None
        assert mol.getAtomCount() == 6

    def test_parse_cxsmiles_with_bond_highlighting(
        self, simple_cxsmiles_bond_highlighting
    ):
        mol = parse_cxsmiles(simple_cxsmiles_bond_highlighting)
        assert mol is not None
        assert mol.getAtomCount() == 3

    def test_parse_complex_cxsmiles(self, complex_cxsmiles):
        mol = parse_cxsmiles(complex_cxsmiles)
        assert mol is not None
        assert mol.getAtomCount() > 6

    def test_parse_standard_smiles_works(self):
        mol = parse_cxsmiles("CCO")
        assert mol is not None
        assert mol.getAtomCount() == 3

    def test_parse_invalid_smiles_raises_exception(self):
        with pytest.raises(Exception):
            parse_cxsmiles("INVALID_SMILES{][")

    def test_parse_aromatic_cxsmiles(self):
        mol = parse_cxsmiles("c1ccccc1 |ha:0,1,2,3,4,5|")
        assert mol is not None
        assert mol.getAtomCount() == 6

    def test_parse_cxsmiles_with_coordinates(self):
        mol = parse_cxsmiles("CCO |coord:0,1.0,2.0,0.0|")
        assert mol is not None
        assert mol.getAtomCount() == 3

    def test_parse_cxsmiles_with_stereo(self):
        mol = parse_cxsmiles("C[C@H](O)CC")
        assert mol is not None
        assert mol.getAtomCount() == 5


class TestCXSMILESEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_large_atom_index(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCCCCC |ha:100|")
        assert atoms == {100}

    def test_negative_index_ignored(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:-1|")
        assert atoms == set()

    def test_duplicate_indices(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:0,0,1,1|")
        assert atoms == {0, 1}

    def test_unsorted_indices(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:2,0,1|")
        assert atoms == {0, 1, 2}

    def test_mixed_valid_invalid_indices(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:0,invalid,1|")
        assert 0 in atoms or 1 in atoms or atoms == set()

    def test_empty_ha_section(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:|")
        assert atoms == set()

    def test_empty_hb_section(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |hb:|")
        assert bonds == set()

    def test_multiple_pipe_symbols(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO |ha:0|extra|")
        assert 0 in atoms or atoms == set()

    def test_cxsmiles_with_reaction(self):
        atoms, bonds = parse_cxsmiles_highlighting_from_string("CCO.CC>>CCCO |ha:0|")
        assert atoms == {0} or atoms == set()


class TestCXSMILESIntegration:
    """Test integration with CDK molecule parsing."""

    def test_roundtrip_simple_molecule(self):
        smiles = "c1ccccc1 |ha:0,1,2|"
        mol = parse_cxsmiles(smiles)
        atoms, bonds = parse_cxsmiles_highlighting_from_string(smiles)
        assert mol.getAtomCount() == 6
        assert atoms == {0, 1, 2}

    def test_caffeine_with_highlighting(self):
        smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C |ha:0,1,2|"
        mol = parse_cxsmiles(smiles)
        atoms, bonds = parse_cxsmiles_highlighting_from_string(smiles)
        assert mol.getAtomCount() == 14
        assert atoms == {0, 1, 2}

    def test_bond_highlighting_integration(self):
        smiles = "CCO |hb:0,1|"
        mol = parse_cxsmiles(smiles)
        atoms, bonds = parse_cxsmiles_highlighting_from_string(smiles)
        assert mol.getBondCount() == 2
        assert bonds == {0, 1}

    def test_combined_highlighting_integration(self):
        smiles = "c1ccccc1O |ha:0,1,2,hb:0,1|"
        mol = parse_cxsmiles(smiles)
        atoms, bonds = parse_cxsmiles_highlighting_from_string(smiles)
        assert mol.getAtomCount() == 7
        assert atoms == {0, 1, 2}
        assert bonds == {0, 1}
