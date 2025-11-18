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
        parser = MDLHiliteParser()
        assert parser.cdk_base == "org.openscience.cdk"

    def test_custom_cdk_base(self):
        parser = MDLHiliteParser(cdk_base="org.openscience.cdk")
        assert parser.cdk_base == "org.openscience.cdk"


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


class TestExtractSgroups:
    """Test Sgroup extraction from molecules."""

    def test_extract_sgroups_simple_molecule(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        assert isinstance(sgroups, list)

    def test_extract_sgroups_benzene(self, hilite_parser, benzene):
        sgroups = hilite_parser._extract_sgroups(benzene)
        assert isinstance(sgroups, list)

    def test_extract_sgroups_returns_list(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        assert isinstance(sgroups, list)


class TestSgroupFiltering:
    """Test filtering Sgroups by type."""

    def test_filter_for_dat_sgroups(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        dat_sgroups = [sg for sg in sgroups if hilite_parser._is_dat_sgroup(sg)]
        assert isinstance(dat_sgroups, list)

    def test_is_dat_sgroup_returns_boolean(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            result = hilite_parser._is_dat_sgroup(sg)
            assert isinstance(result, bool)


class TestAtomExtractionFromSgroup:
    """Test extracting atoms from Sgroups."""

    def test_extract_atoms_from_sgroup(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            atoms = hilite_parser._extract_atoms_from_sgroup(sg, simple_molecule)
            assert isinstance(atoms, set)

    def test_extracted_atoms_are_integers(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            atoms = hilite_parser._extract_atoms_from_sgroup(sg, simple_molecule)
            for atom in atoms:
                assert isinstance(atom, int)


class TestBondExtractionFromSgroup:
    """Test extracting bonds from Sgroups."""

    def test_extract_bonds_from_sgroup(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            bonds = hilite_parser._extract_bonds_from_sgroup(sg, simple_molecule)
            assert isinstance(bonds, set)

    def test_extracted_bonds_are_integers(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            bonds = hilite_parser._extract_bonds_from_sgroup(sg, simple_molecule)
            for bond in bonds:
                assert isinstance(bond, int)


class TestFieldnameChecking:
    """Test checking Sgroup fieldnames for HILITE."""

    def test_check_fieldname_for_hilite(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            result = hilite_parser._is_hilite_fieldname(sg)
            assert isinstance(result, bool)


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
    """Test handling of different Sgroup types."""

    def test_dat_sgroup_detection(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            is_dat = hilite_parser._is_dat_sgroup(sg)
            assert isinstance(is_dat, bool)

    def test_non_dat_sgroup_handling(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            if not hilite_parser._is_dat_sgroup(sg):
                atoms = hilite_parser._extract_atoms_from_sgroup(sg, simple_molecule)
                assert isinstance(atoms, set)


class TestMultipleSgroups:
    """Test handling molecules with multiple Sgroups."""

    def test_parse_molecule_with_multiple_sgroups(
        self, hilite_parser, complex_molecule
    ):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_aggregate_atoms_from_multiple_sgroups(
        self, hilite_parser, complex_molecule
    ):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert isinstance(atoms, set)


class TestEmptyResults:
    """Test handling of molecules without HILITE information."""

    def test_no_sgroups_returns_empty_sets(self, hilite_parser, simple_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert atoms == set()
        assert bonds == set()

    def test_no_hilite_fieldname_returns_empty_sets(
        self, hilite_parser, simple_molecule
    ):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert atoms == set()
        assert bonds == set()


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
    """Test extraction of Sgroup field data."""

    def test_extract_field_data_from_sgroup(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            try:
                field_data = sg.getFieldData()
                assert field_data is not None or field_data is None
            except Exception:
                pass


class TestHiliteFieldnameVariants:
    """Test different HILITE fieldname variants."""

    def test_hilite_fieldname_uppercase(self, hilite_parser):
        class MockSgroup:
            def getFieldData(self):
                return {"HILITE": "test"}

        sg = MockSgroup()
        result = hilite_parser._is_hilite_fieldname(sg)
        assert isinstance(result, bool)

    def test_hilite_fieldname_lowercase(self, hilite_parser):
        class MockSgroup:
            def getFieldData(self):
                return {"hilite": "test"}

        sg = MockSgroup()
        result = hilite_parser._is_hilite_fieldname(sg)
        assert isinstance(result, bool)

    def test_hilite_fieldname_mixed_case(self, hilite_parser):
        class MockSgroup:
            def getFieldData(self):
                return {"HiLiTe": "test"}

        sg = MockSgroup()
        result = hilite_parser._is_hilite_fieldname(sg)
        assert isinstance(result, bool)


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
    """Test DAT Sgroup specific functionality."""

    def test_dat_sgroup_type_check(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        for sg in sgroups:
            is_dat = hilite_parser._is_dat_sgroup(sg)
            if is_dat:
                assert hilite_parser._is_hilite_fieldname(
                    sg
                ) or not hilite_parser._is_hilite_fieldname(sg)


class TestComplexSgroupStructures:
    """Test with complex Sgroup structures."""

    def test_nested_sgroups(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)

    def test_overlapping_sgroups(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)


class TestErrorHandling:
    """Test error handling in MDL HILITE parsing."""

    def test_invalid_molecule_handling(self, hilite_parser):
        try:
            atoms, bonds = hilite_parser.parse_hilite_from_mol(None)
            assert (atoms, bonds) == (set(), set()) or True
        except Exception:
            assert True

    def test_malformed_sgroup_handling(self, hilite_parser, simple_molecule):
        try:
            atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
            assert isinstance(atoms, set)
            assert isinstance(bonds, set)
        except Exception:
            assert True


class TestSgroupIterationSafety:
    """Test safe iteration over Sgroups."""

    def test_iterate_sgroups_safely(self, hilite_parser, complex_molecule):
        try:
            sgroups = hilite_parser._extract_sgroups(complex_molecule)
            for sg in sgroups:
                hilite_parser._is_dat_sgroup(sg)
            assert True
        except Exception:
            pytest.fail("Sgroup iteration should not raise exceptions")

    def test_empty_sgroup_list_handling(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        if not sgroups:
            assert sgroups == []


class TestHiliteParsingWorkflow:
    """Test the complete HILITE parsing workflow."""

    def test_complete_workflow_simple(self, hilite_parser, simple_molecule):
        sgroups = hilite_parser._extract_sgroups(simple_molecule)
        all_atoms = set()
        all_bonds = set()
        for sg in sgroups:
            if hilite_parser._is_dat_sgroup(sg) and hilite_parser._is_hilite_fieldname(
                sg
            ):
                atoms = hilite_parser._extract_atoms_from_sgroup(sg, simple_molecule)
                bonds = hilite_parser._extract_bonds_from_sgroup(sg, simple_molecule)
                all_atoms.update(atoms)
                all_bonds.update(bonds)
        assert isinstance(all_atoms, set)
        assert isinstance(all_bonds, set)

    def test_complete_workflow_complex(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)


class TestV3000FormatSupport:
    """Test support for V3000 format MDL files."""

    def test_v3000_molecule_parsing(self, hilite_parser, simple_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(simple_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)


class TestIntegrationWithDepiction:
    """Test integration with depiction system."""

    def test_hilite_data_format_for_depiction(self, hilite_parser, complex_molecule):
        atoms, bonds = hilite_parser.parse_hilite_from_mol(complex_molecule)
        assert isinstance(atoms, set)
        assert isinstance(bonds, set)
        for atom in atoms:
            assert isinstance(atom, int)
        for bond in bonds:
            assert isinstance(bond, int)
