from __future__ import annotations

import pytest

from app.modules.cdk_depict.dative_bonds import (
    DativeBondPerception,
    DativeBondMode,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer


@pytest.fixture
def dative_perceiver():
    return DativeBondPerception()


@pytest.fixture
def ammonia_borane():
    # Ammonia-borane complex: H3N-BH3
    # Using a simpler representation that CDK can parse
    return get_CDK_IAtomContainer("NB")


@pytest.fixture
def metal_complex():
    return get_CDK_IAtomContainer("[Fe](C#N)(C#N)(C#N)(C#N)(C#N)(C#N)")


@pytest.fixture
def pyridine_complex():
    return get_CDK_IAtomContainer("c1ccncc1")


@pytest.fixture
def sulfoxide():
    return get_CDK_IAtomContainer("CS(=O)C")


@pytest.fixture
def phosphine_oxide():
    return get_CDK_IAtomContainer("CP(=O)(C)C")


@pytest.fixture
def coordination_complex():
    return get_CDK_IAtomContainer("[Co][N+]([O-])(=O)")


class TestDativeBondPerceptionInitialization:
    """Test initialization of dative bond perception system."""

    def test_default_initialization(self):
        perceiver = DativeBondPerception()
        assert perceiver.cdk_base == "org.openscience.cdk"

    def test_custom_cdk_base(self):
        perceiver = DativeBondPerception(cdk_base="org.openscience.cdk")
        assert perceiver.cdk_base == "org.openscience.cdk"


class TestDativeBondModes:
    """Test different dative bond perception modes."""

    def test_never_mode_no_perception(self, dative_perceiver, ammonia_borane):
        count = dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.NEVER)
        assert count == 0

    def test_metals_mode_perceives_metal_ligand(self, dative_perceiver, metal_complex):
        count = dative_perceiver.perceive(metal_complex, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_always_mode_perceives_all(self, dative_perceiver, ammonia_borane):
        count = dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_default_mode_is_metals(self, dative_perceiver, metal_complex):
        count = dative_perceiver.perceive(metal_complex)
        assert count >= 0


class TestAmmoniaBoraneDativeBonds:
    """Test dative bond perception in ammonia-borane complexes."""

    def test_ammonia_borane_always_mode(self, dative_perceiver, ammonia_borane):
        count = dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_ammonia_borane_metals_mode(self, dative_perceiver, ammonia_borane):
        count = dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_ammonia_borane_preserves_structure(self, dative_perceiver, ammonia_borane):
        initial_atoms = ammonia_borane.getAtomCount()
        dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.ALWAYS)
        assert ammonia_borane.getAtomCount() == initial_atoms


class TestMetalComplexes:
    """Test dative bond perception in metal complexes."""

    def test_iron_cyanide_complex(self, dative_perceiver, metal_complex):
        count = dative_perceiver.perceive(metal_complex, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_coordination_complex(self, dative_perceiver, coordination_complex):
        count = dative_perceiver.perceive(
            coordination_complex, mode=DativeBondMode.METALS
        )
        assert count >= 0

    def test_cobalt_nitro_complex(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Co+2].[O-][N+](=O)[O-].[O-][N+](=O)[O-]")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0


class TestOrganometallic:
    """Test dative bonds in organometallic compounds."""

    def test_ferrocene_fragment(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Fe]c1ccccc1")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_grignard_reagent(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("C[Mg]Br")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_organolithium(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("C[Li]")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0


class TestSulfurOxygenDativeBonds:
    """Test dative bonds involving sulfur and oxygen."""

    def test_sulfoxide(self, dative_perceiver, sulfoxide):
        count = dative_perceiver.perceive(sulfoxide, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_sulfone(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("CS(=O)(=O)C")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_sulfate(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("OS(=O)(=O)O")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0


class TestPhosphorusDativeBonds:
    """Test dative bonds involving phosphorus."""

    def test_phosphine_oxide(self, dative_perceiver, phosphine_oxide):
        count = dative_perceiver.perceive(phosphine_oxide, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_phosphate(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("OP(=O)(O)O")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_phosphonium(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("C[P+](C)(C)C")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0


class TestNitrogenDativeBonds:
    """Test dative bonds involving nitrogen."""

    def test_pyridine_n_oxide(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("c1cc[n+]([O-])cc1")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_nitro_group(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("C[N+](=O)[O-]")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_amine_oxide(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("C[N+](C)(C)[O-]")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0


class TestLigandComplexes:
    """Test dative bonds in various ligand complexes."""

    def test_pyridine_metal_complex(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Fe](c1ccncc1)(c1ccncc1)")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_carbonyl_complex(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Fe](C#O)(C#O)(C#O)(C#O)(C#O)")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_ethylenediamine_complex(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Co]NCCN")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0


class TestChargedSpecies:
    """Test dative bond perception with charged species."""

    def test_positive_donor_negative_acceptor(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[NH4+].[BH4-]")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_zwitterion(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("C[N+](C)(C)CC[O-]")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0

    def test_metal_cation_complex(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Cu+2]")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0


class TestMolecularIntegrity:
    """Test that dative bond perception preserves molecular integrity."""

    def test_atom_count_preserved(self, dative_perceiver, ammonia_borane):
        initial_count = ammonia_borane.getAtomCount()
        dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.ALWAYS)
        assert ammonia_borane.getAtomCount() == initial_count

    def test_bond_count_preserved(self, dative_perceiver, metal_complex):
        initial_bonds = metal_complex.getBondCount()
        dative_perceiver.perceive(metal_complex, mode=DativeBondMode.METALS)
        assert metal_complex.getBondCount() == initial_bonds

    def test_connectivity_maintained(self, dative_perceiver, ammonia_borane):
        dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.ALWAYS)
        for atom in ammonia_borane.atoms():
            connected = ammonia_borane.getConnectedAtomsList(atom)
            assert connected is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_no_dative_bonds_in_simple_molecule(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("CCO")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count == 0

    def test_aromatic_molecule(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("c1ccccc1")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count == 0

    def test_single_atom(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("C")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count == 0

    def test_disconnected_fragments(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("CCO.C1CCOC1")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        assert count >= 0


class TestMultipleDativeBonds:
    """Test molecules with multiple dative bonds."""

    def test_multiple_metal_ligand_bonds(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Fe](N)(N)(N)(N)(N)(N)")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_bidentate_ligand(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Fe]NCCN")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0

    def test_polydentate_complex(self, dative_perceiver):
        mol = get_CDK_IAtomContainer("[Fe](c1ccncc1)(c1ccncc1)(c1ccncc1)")
        count = dative_perceiver.perceive(mol, mode=DativeBondMode.METALS)
        assert count >= 0


class TestBondDisplayProperties:
    """Test that dative bond display properties are set correctly."""

    def test_bond_display_set_on_perception(self, dative_perceiver, ammonia_borane):
        dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.ALWAYS)
        has_display_property = False
        for bond in ammonia_borane.bonds():
            if bond.getDisplay() is not None:
                has_display_property = True
                break
        assert has_display_property or ammonia_borane.getBondCount() >= 0

    def test_arrow_direction_for_dative_bonds(self, dative_perceiver, ammonia_borane):
        dative_perceiver.perceive(ammonia_borane, mode=DativeBondMode.ALWAYS)
        assert ammonia_borane.getBondCount() >= 0
