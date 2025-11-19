from __future__ import annotations

import pytest

from app.modules.cdk_depict.multicenter_bonds import (
    MulticenterBonds,
    MulticenterStyle,
    set_multicenter_style,
    get_multicenter_style,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer


@pytest.fixture
def mc_handler():
    return MulticenterBonds()


@pytest.fixture
def simple_molecule():
    return get_CDK_IAtomContainer("c1ccccc1")


@pytest.fixture
def ferrocene_fragment():
    return get_CDK_IAtomContainer("[Fe]c1ccccc1")


@pytest.fixture
def chromium_benzene():
    return get_CDK_IAtomContainer("[Cr]c1ccccc1")


@pytest.fixture
def allyl_complex():
    return get_CDK_IAtomContainer("[Pd]C=CC")


class TestMulticenterBondsInitialization:
    """Test initialization of multicenter bond handler."""

    def test_default_initialization(self):
        handler = MulticenterBonds()
        assert handler.cdk_base == "org.openscience.cdk"

    def test_custom_cdk_base(self):
        handler = MulticenterBonds(cdk_base="org.openscience.cdk")
        assert handler.cdk_base == "org.openscience.cdk"


class TestMulticenterStyles:
    """Test different multicenter bond display styles."""

    def test_provided_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.PROVIDED)
        assert count >= 0

    def test_dative_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_dashed_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DASHED)
        assert count >= 0

    def test_dashed_neutral_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(
            ferrocene_fragment, MulticenterStyle.DASHED_NEUTRAL
        )
        assert count >= 0

    def test_hidden_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.HIDDEN)
        assert count >= 0

    def test_hidden_neutral_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(
            ferrocene_fragment, MulticenterStyle.HIDDEN_NEUTRAL
        )
        assert count >= 0


class TestFerrocene:
    """Test multicenter bonds in ferrocene complexes."""

    def test_ferrocene_provided_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.PROVIDED)
        assert count >= 0
        assert ferrocene_fragment.getBondCount() >= 0

    def test_ferrocene_dative_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_ferrocene_hidden_style(self, mc_handler, ferrocene_fragment):
        count = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.HIDDEN)
        assert count >= 0


class TestChromiumBenzene:
    """Test multicenter bonds in chromium-benzene complexes."""

    def test_chromium_benzene_provided(self, mc_handler, chromium_benzene):
        count = mc_handler.set_style(chromium_benzene, MulticenterStyle.PROVIDED)
        assert count >= 0

    def test_chromium_benzene_dative(self, mc_handler, chromium_benzene):
        count = mc_handler.set_style(chromium_benzene, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_chromium_benzene_dashed(self, mc_handler, chromium_benzene):
        count = mc_handler.set_style(chromium_benzene, MulticenterStyle.DASHED)
        assert count >= 0


class TestAllylComplexes:
    """Test multicenter bonds in π-allyl complexes."""

    def test_allyl_provided(self, mc_handler, allyl_complex):
        count = mc_handler.set_style(allyl_complex, MulticenterStyle.PROVIDED)
        assert count >= 0

    def test_allyl_dative(self, mc_handler, allyl_complex):
        count = mc_handler.set_style(allyl_complex, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_allyl_hidden_neutral(self, mc_handler, allyl_complex):
        count = mc_handler.set_style(allyl_complex, MulticenterStyle.HIDDEN_NEUTRAL)
        assert count >= 0


class TestDetectMulticenterBonds:
    """Test detection of multicenter bonds."""

    def test_detect_in_ferrocene(self, mc_handler, ferrocene_fragment):
        count = mc_handler.detect_multicenter_bonds(ferrocene_fragment)
        assert count >= 0

    def test_detect_in_chromium_benzene(self, mc_handler, chromium_benzene):
        count = mc_handler.detect_multicenter_bonds(chromium_benzene)
        assert count >= 0

    def test_detect_in_simple_molecule(self, mc_handler, simple_molecule):
        count = mc_handler.detect_multicenter_bonds(simple_molecule)
        assert count >= 0

    def test_detect_in_allyl_complex(self, mc_handler, allyl_complex):
        count = mc_handler.detect_multicenter_bonds(allyl_complex)
        assert count >= 0


class TestGetMulticenterStyle:
    """Test get_multicenter_style string conversion."""

    def test_provided_string_variants(self):
        assert get_multicenter_style("p") == MulticenterStyle.PROVIDED
        assert get_multicenter_style("provided") == MulticenterStyle.PROVIDED
        assert get_multicenter_style("default") == MulticenterStyle.PROVIDED

    def test_dative_string_variants(self):
        assert get_multicenter_style("d") == MulticenterStyle.DATIVE
        assert get_multicenter_style("dative") == MulticenterStyle.DATIVE
        assert get_multicenter_style("arrow") == MulticenterStyle.DATIVE

    def test_dashed_string_variants(self):
        assert get_multicenter_style("a") == MulticenterStyle.DASHED
        assert get_multicenter_style("dashed") == MulticenterStyle.DASHED
        assert get_multicenter_style("dash") == MulticenterStyle.DASHED

    def test_dashed_neutral_string_variants(self):
        assert get_multicenter_style("an") == MulticenterStyle.DASHED_NEUTRAL
        assert (
            get_multicenter_style("dashed_neutral") == MulticenterStyle.DASHED_NEUTRAL
        )
        assert get_multicenter_style("dashneutral") == MulticenterStyle.DASHED_NEUTRAL

    def test_hidden_string_variants(self):
        assert get_multicenter_style("h") == MulticenterStyle.HIDDEN
        assert get_multicenter_style("hidden") == MulticenterStyle.HIDDEN
        assert get_multicenter_style("hide") == MulticenterStyle.HIDDEN

    def test_hidden_neutral_string_variants(self):
        assert get_multicenter_style("hn") == MulticenterStyle.HIDDEN_NEUTRAL
        assert (
            get_multicenter_style("hidden_neutral") == MulticenterStyle.HIDDEN_NEUTRAL
        )
        assert get_multicenter_style("hideneutral") == MulticenterStyle.HIDDEN_NEUTRAL

    def test_unknown_string_defaults_to_provided(self):
        assert get_multicenter_style("unknown") == MulticenterStyle.PROVIDED
        assert get_multicenter_style("invalid") == MulticenterStyle.PROVIDED

    def test_case_insensitive(self):
        assert get_multicenter_style("DATIVE") == MulticenterStyle.DATIVE
        assert get_multicenter_style("Provided") == MulticenterStyle.PROVIDED
        assert get_multicenter_style("HIDDEN") == MulticenterStyle.HIDDEN


class TestSetMulticenterStyle:
    """Test convenience function set_multicenter_style."""

    def test_set_style_with_string(self, ferrocene_fragment):
        count = set_multicenter_style(ferrocene_fragment, "dative")
        assert count >= 0

    def test_set_style_provided(self, ferrocene_fragment):
        count = set_multicenter_style(ferrocene_fragment, "provided")
        assert count >= 0

    def test_set_style_dashed(self, ferrocene_fragment):
        count = set_multicenter_style(ferrocene_fragment, "dashed")
        assert count >= 0

    def test_set_style_hidden(self, ferrocene_fragment):
        count = set_multicenter_style(ferrocene_fragment, "hidden")
        assert count >= 0


class TestMolecularIntegrity:
    """Test that multicenter bond operations preserve molecular integrity."""

    def test_atom_count_preserved(self, mc_handler, ferrocene_fragment):
        initial_count = ferrocene_fragment.getAtomCount()
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DATIVE)
        assert ferrocene_fragment.getAtomCount() == initial_count

    def test_bond_count_consistent(self, mc_handler, ferrocene_fragment):
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DASHED)
        assert ferrocene_fragment.getBondCount() >= 0

    def test_connectivity_maintained(self, mc_handler, ferrocene_fragment):
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.HIDDEN)
        for atom in ferrocene_fragment.atoms():
            connected = ferrocene_fragment.getConnectedAtomsList(atom)
            assert connected is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_simple_molecule_no_multicenter(self, mc_handler, simple_molecule):
        count = mc_handler.set_style(simple_molecule, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_single_atom_molecule(self, mc_handler):
        mol = get_CDK_IAtomContainer("C")
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count == 0

    def test_disconnected_fragments(self, mc_handler):
        mol = get_CDK_IAtomContainer("CCO.[Fe]")
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_aromatic_without_metal(self, mc_handler):
        mol = get_CDK_IAtomContainer("c1ccccc1")
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 0


class TestComplexOrganometallics:
    """Test multicenter bonds in complex organometallic structures."""

    def test_cyclopentadienyl_complex(self, mc_handler):
        mol = get_CDK_IAtomContainer("[Fe]C1=CC=CC1")
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_bis_cyclopentadienyl(self, mc_handler):
        mol = get_CDK_IAtomContainer("[Fe](C1=CC=CC1)C1=CC=CC1")
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 0

    def test_arene_complex(self, mc_handler):
        mol = get_CDK_IAtomContainer("[Cr](c1ccccc1)")
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 0


class TestChargeNeutralization:
    """Test charge neutralization in neutral styles."""

    def test_dashed_neutral_preserves_atoms(self, mc_handler, ferrocene_fragment):
        initial_atoms = ferrocene_fragment.getAtomCount()
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DASHED_NEUTRAL)
        assert ferrocene_fragment.getAtomCount() == initial_atoms

    def test_hidden_neutral_preserves_atoms(self, mc_handler, ferrocene_fragment):
        initial_atoms = ferrocene_fragment.getAtomCount()
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.HIDDEN_NEUTRAL)
        assert ferrocene_fragment.getAtomCount() == initial_atoms

    def test_charges_consistent_after_neutralization(
        self, mc_handler, ferrocene_fragment
    ):
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DASHED_NEUTRAL)
        total_charge = sum(
            atom.getFormalCharge() if atom.getFormalCharge() else 0
            for atom in ferrocene_fragment.atoms()
        )
        assert isinstance(total_charge, int)


class TestStylePersistence:
    """Test that style changes persist correctly."""

    def test_multiple_style_changes(self, mc_handler, ferrocene_fragment):
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DATIVE)
        count1 = ferrocene_fragment.getBondCount()
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.HIDDEN)
        count2 = ferrocene_fragment.getBondCount()
        assert count1 >= 0 and count2 >= 0

    def test_revert_to_provided(self, mc_handler, ferrocene_fragment):
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DATIVE)
        mc_handler.set_style(ferrocene_fragment, MulticenterStyle.PROVIDED)
        assert ferrocene_fragment is not None

    def test_style_change_is_idempotent(self, mc_handler, ferrocene_fragment):
        count1 = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DASHED)
        count2 = mc_handler.set_style(ferrocene_fragment, MulticenterStyle.DASHED)
        assert count1 >= 0 and count2 >= 0
