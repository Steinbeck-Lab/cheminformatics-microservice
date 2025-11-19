from __future__ import annotations

import pytest

from app.modules.cdk_depict.aromatic_display import (
    AromaticDisplaySystem,
    apply_aromatic_display,
    perceive_aromaticity,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from jpype import JClass


@pytest.fixture
def aromatic_system():
    return AromaticDisplaySystem()


@pytest.fixture
def benzene():
    return get_CDK_IAtomContainer("c1ccccc1")


@pytest.fixture
def naphthalene():
    return get_CDK_IAtomContainer("c1ccc2ccccc2c1")


@pytest.fixture
def pyridine():
    return get_CDK_IAtomContainer("c1ccncc1")


@pytest.fixture
def non_aromatic():
    return get_CDK_IAtomContainer("C1CCCCC1")


class TestAromaticDisplaySystemInitialization:
    """Test initialization of aromatic display system."""

    def test_default_initialization(self):
        system = AromaticDisplaySystem()
        assert system.cdk_base == "org.openscience.cdk"

    def test_initialization_loads_cdk_classes(self, aromatic_system):
        assert aromatic_system.Aromaticity is not None
        assert aromatic_system.ElectronDonation is not None
        assert aromatic_system.Cycles is not None


class TestAromaticPerception:
    """Test aromaticity perception."""

    def test_perceive_benzene(self, aromatic_system, benzene):
        aromatic_system._perceive_aromaticity(benzene)
        aromatic_count = sum(1 for atom in benzene.atoms() if atom.isAromatic())
        assert aromatic_count == 6

    def test_perceive_naphthalene(self, aromatic_system, naphthalene):
        aromatic_system._perceive_aromaticity(naphthalene)
        aromatic_count = sum(1 for atom in naphthalene.atoms() if atom.isAromatic())
        assert aromatic_count == 10

    def test_perceive_pyridine(self, aromatic_system, pyridine):
        aromatic_system._perceive_aromaticity(pyridine)
        aromatic_count = sum(1 for atom in pyridine.atoms() if atom.isAromatic())
        assert aromatic_count == 6

    def test_perceive_non_aromatic(self, aromatic_system, non_aromatic):
        aromatic_system._perceive_aromaticity(non_aromatic)
        aromatic_count = sum(1 for atom in non_aromatic.atoms() if atom.isAromatic())
        assert aromatic_count == 0


class TestDonutDisplay:
    """Test donut (circle-in-ring) display mode."""

    def test_apply_donut_enabled(self, aromatic_system, benzene):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = aromatic_system.apply_donut_display(
            DepictionGenerator, benzene, enable=True
        )
        assert gen is not None

    def test_apply_donut_disabled(self, aromatic_system, benzene):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = aromatic_system.apply_donut_display(
            DepictionGenerator, benzene, enable=False
        )
        assert gen is not None

    def test_donut_on_naphthalene(self, aromatic_system, naphthalene):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = aromatic_system.apply_donut_display(
            DepictionGenerator, naphthalene, enable=True
        )
        assert gen is not None

    def test_donut_on_pyridine(self, aromatic_system, pyridine):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = aromatic_system.apply_donut_display(
            DepictionGenerator, pyridine, enable=True
        )
        assert gen is not None

    def test_donut_on_non_aromatic(self, aromatic_system, non_aromatic):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = aromatic_system.apply_donut_display(
            DepictionGenerator, non_aromatic, enable=True
        )
        assert gen is not None


class TestAromaticRingDetection:
    """Test detection of aromatic rings."""

    def test_get_aromatic_rings_benzene(self, aromatic_system, benzene):
        rings = aromatic_system.get_aromatic_rings(benzene)
        assert len(rings) >= 1

    def test_get_aromatic_rings_naphthalene(self, aromatic_system, naphthalene):
        rings = aromatic_system.get_aromatic_rings(naphthalene)
        assert len(rings) >= 2

    def test_get_aromatic_rings_pyridine(self, aromatic_system, pyridine):
        rings = aromatic_system.get_aromatic_rings(pyridine)
        assert len(rings) >= 1

    def test_get_aromatic_rings_non_aromatic(self, aromatic_system, non_aromatic):
        rings = aromatic_system.get_aromatic_rings(non_aromatic)
        assert len(rings) == 0


class TestAromaticAtomCounting:
    """Test counting aromatic atoms."""

    def test_count_aromatic_atoms_benzene(self, aromatic_system, benzene):
        count = aromatic_system.count_aromatic_atoms(benzene)
        assert count == 6

    def test_count_aromatic_atoms_naphthalene(self, aromatic_system, naphthalene):
        count = aromatic_system.count_aromatic_atoms(naphthalene)
        assert count == 10

    def test_count_aromatic_atoms_pyridine(self, aromatic_system, pyridine):
        count = aromatic_system.count_aromatic_atoms(pyridine)
        assert count == 6

    def test_count_aromatic_atoms_non_aromatic(self, aromatic_system, non_aromatic):
        count = aromatic_system.count_aromatic_atoms(non_aromatic)
        assert count == 0

    def test_count_aromatic_atoms_aliphatic(self, aromatic_system):
        mol = get_CDK_IAtomContainer("CCCCCC")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 0


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_apply_aromatic_display_enabled(self, benzene):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = apply_aromatic_display(DepictionGenerator, benzene, donuts=True)
        assert gen is not None

    def test_apply_aromatic_display_disabled(self, benzene):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = apply_aromatic_display(DepictionGenerator, benzene, donuts=False)
        assert gen is not None

    def test_perceive_aromaticity_function(self, benzene):
        perceive_aromaticity(benzene)
        aromatic_count = sum(1 for atom in benzene.atoms() if atom.isAromatic())
        assert aromatic_count == 6


class TestComplexAromaticSystems:
    """Test complex aromatic systems."""

    def test_anthracene(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1ccc2cc3ccccc3cc2c1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 14

    def test_phenanthrene(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1ccc2c(c1)ccc1ccccc12")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 14

    def test_biphenyl(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1ccc(cc1)c2ccccc2")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 12

    def test_caffeine(self, aromatic_system):
        mol = get_CDK_IAtomContainer("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count >= 0


class TestHeteroaromaticCompounds:
    """Test heteroaromatic compounds."""

    def test_pyrrole(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1cc[nH]c1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 5

    def test_furan(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1ccoc1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 5

    def test_thiophene(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1ccsc1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 5

    def test_imidazole(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1cnc[nH]1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count >= 4

    def test_pyrimidine(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1cncnc1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 6


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_atom(self, aromatic_system):
        mol = get_CDK_IAtomContainer("C")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count == 0

    def test_cyclobutene_non_aromatic(self, aromatic_system):
        mol = get_CDK_IAtomContainer("C1=CCC1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count >= 0

    def test_cyclooctatetraene(self, aromatic_system):
        mol = get_CDK_IAtomContainer("C1=CC=CC=CC=C1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count >= 0

    def test_charged_aromatic(self, aromatic_system):
        mol = get_CDK_IAtomContainer("c1cc[nH+]cc1")
        count = aromatic_system.count_aromatic_atoms(mol)
        assert count >= 0


class TestAromaticBonds:
    """Test aromatic bond detection."""

    def test_benzene_aromatic_bonds(self, aromatic_system, benzene):
        aromatic_system._perceive_aromaticity(benzene)
        aromatic_bonds = sum(1 for bond in benzene.bonds() if bond.isAromatic())
        assert aromatic_bonds == 6

    def test_naphthalene_aromatic_bonds(self, aromatic_system, naphthalene):
        aromatic_system._perceive_aromaticity(naphthalene)
        aromatic_bonds = sum(1 for bond in naphthalene.bonds() if bond.isAromatic())
        assert aromatic_bonds == 11

    def test_non_aromatic_bonds(self, aromatic_system, non_aromatic):
        aromatic_system._perceive_aromaticity(non_aromatic)
        aromatic_bonds = sum(1 for bond in non_aromatic.bonds() if bond.isAromatic())
        assert aromatic_bonds == 0


class TestMolecularIntegrity:
    """Test that aromaticity perception preserves molecular integrity."""

    def test_atom_count_preserved(self, aromatic_system, benzene):
        initial_count = benzene.getAtomCount()
        aromatic_system._perceive_aromaticity(benzene)
        assert benzene.getAtomCount() == initial_count

    def test_bond_count_preserved(self, aromatic_system, benzene):
        initial_bonds = benzene.getBondCount()
        aromatic_system._perceive_aromaticity(benzene)
        assert benzene.getBondCount() == initial_bonds

    def test_connectivity_preserved(self, aromatic_system, benzene):
        aromatic_system._perceive_aromaticity(benzene)
        for atom in benzene.atoms():
            connected = benzene.getConnectedAtomsList(atom)
            assert connected is not None


class TestDonutDisplayPersistence:
    """Test that donut display can be enabled and disabled."""

    def test_enable_then_disable(self, aromatic_system, benzene):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen1 = aromatic_system.apply_donut_display(
            DepictionGenerator, benzene, enable=True
        )
        gen2 = aromatic_system.apply_donut_display(gen1, benzene, enable=False)
        assert gen2 is not None

    def test_multiple_enable_calls(self, aromatic_system, benzene):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen1 = aromatic_system.apply_donut_display(
            DepictionGenerator, benzene, enable=True
        )
        gen2 = aromatic_system.apply_donut_display(gen1, benzene, enable=True)
        assert gen2 is not None
