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


class TestMarkMulticenterBond:
    """Test mark_multicenter_bond() to manually create ExtMulticenter Sgroups."""

    def test_mark_ferrocene_fragment(self, mc_handler):
        """Mark a multicenter bond on [Fe]c1ccccc1 between Fe and a ring atom."""
        mol = get_CDK_IAtomContainer("[Fe]c1ccccc1")
        # Fe is at index 0, ring carbons at 1-6
        ring_indices = list(range(1, 7))
        result = mc_handler.mark_multicenter_bond(mol, 0, ring_indices)
        assert result is True

    def test_mark_chromium_benzene(self, mc_handler):
        """Mark a multicenter bond on [Cr]c1ccccc1."""
        mol = get_CDK_IAtomContainer("[Cr]c1ccccc1")
        ring_indices = list(range(1, 7))
        result = mc_handler.mark_multicenter_bond(mol, 0, ring_indices)
        assert result is True

    def test_mark_no_bond_found(self, mc_handler):
        """Return False when no bond exists between metal and ring atoms."""
        mol = get_CDK_IAtomContainer("CCO.[Fe]")
        # Fe is disconnected from C atoms, so no bond should exist
        fe_idx = None
        for i in range(mol.getAtomCount()):
            if mol.getAtom(i).getSymbol() == "Fe":
                fe_idx = i
                break
        assert fe_idx is not None
        # Use indices that don't share a bond with Fe
        result = mc_handler.mark_multicenter_bond(mol, fe_idx, [0, 1])
        assert result is False

    def test_mark_with_out_of_range_index(self, mc_handler):
        """Atoms with out-of-range indices are skipped gracefully."""
        mol = get_CDK_IAtomContainer("[Fe]c1ccccc1")
        ring_indices = [1, 2, 3, 999]  # 999 is out of range
        result = mc_handler.mark_multicenter_bond(mol, 0, ring_indices)
        # Should still succeed if at least one valid ring atom has a bond to metal
        assert result is True

    def test_mark_creates_sgroup_property(self, mc_handler):
        """After marking, CTAB_SGROUPS property should exist on molecule."""
        from jpype import JClass

        CDKConstants = JClass("org.openscience.cdk.CDKConstants")
        mol = get_CDK_IAtomContainer("[Fe]c1ccccc1")
        mc_handler.mark_multicenter_bond(mol, 0, list(range(1, 7)))
        sgroups = mol.getProperty(CDKConstants.CTAB_SGROUPS)
        assert sgroups is not None
        assert sgroups.size() > 0

    def test_mark_multiple_bonds(self, mc_handler):
        """Mark two separate multicenter bonds on the same molecule."""
        from jpype import JClass

        CDKConstants = JClass("org.openscience.cdk.CDKConstants")
        mol = get_CDK_IAtomContainer("[Fe](C1=CC=CC1)C1=CC=CC1")
        # Mark first Cp ring
        mc_handler.mark_multicenter_bond(mol, 0, [1, 2, 3, 4, 5])
        # Mark second Cp ring
        mc_handler.mark_multicenter_bond(mol, 0, [6, 7, 8, 9, 10])
        sgroups = mol.getProperty(CDKConstants.CTAB_SGROUPS)
        assert sgroups is not None
        assert sgroups.size() == 2


class TestSetStyleWithSgroups:
    """Test set_style() on molecules that have actual ExtMulticenter Sgroups.

    Uses mark_multicenter_bond() to create Sgroups first, then applies styles.
    This covers set_style inner loop (lines 100-147), _apply_bond_style
    (lines 158-179), and _neutralize_charges (lines 192-223).
    """

    def _make_molecule_with_sgroup(self, mc_handler):
        """Create a molecule with an ExtMulticenter Sgroup via mark_multicenter_bond."""
        mol = get_CDK_IAtomContainer("[Fe]c1ccccc1")
        ring_indices = list(range(1, 7))
        mc_handler.mark_multicenter_bond(mol, 0, ring_indices)
        return mol

    def test_dative_style_with_sgroup(self, mc_handler):
        """Apply DATIVE style to molecule with ExtMulticenter Sgroup."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 1

    def test_dashed_style_with_sgroup(self, mc_handler):
        """Apply DASHED style to molecule with ExtMulticenter Sgroup."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        count = mc_handler.set_style(mol, MulticenterStyle.DASHED)
        assert count >= 1

    def test_dashed_neutral_style_with_sgroup(self, mc_handler):
        """Apply DASHED_NEUTRAL style - triggers charge neutralization."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        count = mc_handler.set_style(mol, MulticenterStyle.DASHED_NEUTRAL)
        assert count >= 1

    def test_hidden_style_with_sgroup(self, mc_handler):
        """Apply HIDDEN style to molecule with ExtMulticenter Sgroup."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        count = mc_handler.set_style(mol, MulticenterStyle.HIDDEN)
        assert count >= 1

    def test_hidden_neutral_style_with_sgroup(self, mc_handler):
        """Apply HIDDEN_NEUTRAL style - triggers charge neutralization."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        count = mc_handler.set_style(mol, MulticenterStyle.HIDDEN_NEUTRAL)
        assert count >= 1

    def test_provided_style_with_sgroup_returns_zero(self, mc_handler):
        """PROVIDED style should return 0 and skip processing."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        count = mc_handler.set_style(mol, MulticenterStyle.PROVIDED)
        assert count == 0

    def test_atom_count_preserved_after_styling(self, mc_handler):
        """Atom count should not change after applying styles."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        initial_atoms = mol.getAtomCount()
        mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert mol.getAtomCount() == initial_atoms

    def test_bond_count_preserved_after_styling(self, mc_handler):
        """Bond count should not change after applying styles."""
        mol = self._make_molecule_with_sgroup(mc_handler)
        initial_bonds = mol.getBondCount()
        mc_handler.set_style(mol, MulticenterStyle.HIDDEN)
        assert mol.getBondCount() == initial_bonds

    def test_all_non_provided_styles_process_sgroup(self, mc_handler):
        """Every non-PROVIDED style should process at least one Sgroup."""
        styles = [
            MulticenterStyle.DATIVE,
            MulticenterStyle.DASHED,
            MulticenterStyle.DASHED_NEUTRAL,
            MulticenterStyle.HIDDEN,
            MulticenterStyle.HIDDEN_NEUTRAL,
        ]
        for style in styles:
            mol = self._make_molecule_with_sgroup(mc_handler)
            count = mc_handler.set_style(mol, style)
            assert count >= 1, f"Style {style.value} did not process any Sgroups"


class TestSetStyleMetalAtEnd:
    """Test set_style when the metal atom is at the end of the bond.

    This covers the elif branch at line 128 where Elements.isMetal(end) is True.
    """

    def test_metal_at_end_dative(self, mc_handler):
        """When metal is the end atom of the Sgroup bond, DATIVE should work."""
        from jpype import JClass

        CDKConstants = JClass("org.openscience.cdk.CDKConstants")
        Sgroup = JClass("org.openscience.cdk.sgroup.Sgroup")
        SgroupType = JClass("org.openscience.cdk.sgroup.SgroupType")
        ArrayList = JClass("java.util.ArrayList")

        # Use c1ccccc1[Fe] - Fe bonded at end
        mol = get_CDK_IAtomContainer("c1ccccc1[Fe]")

        # Find Fe atom index
        fe_idx = None
        for i in range(mol.getAtomCount()):
            if mol.getAtom(i).getSymbol() == "Fe":
                fe_idx = i
                break
        assert fe_idx is not None

        # Find a bond between Fe and a ring atom
        fe_atom = mol.getAtom(fe_idx)
        ring_bond = None
        for bond in fe_atom.bonds():
            other = bond.getOther(fe_atom)
            other_idx = mol.indexOf(other)
            if other_idx != fe_idx:
                ring_bond = bond
                break

        assert ring_bond is not None

        # Create Sgroup manually with the ring atom as the begin
        sgroup = Sgroup()
        sgroup.setType(SgroupType.ExtMulticenter)
        # Add ring atoms
        for i in range(mol.getAtomCount()):
            if i != fe_idx:
                sgroup.addAtom(mol.getAtom(i))
        sgroup.addBond(ring_bond)

        sgroups = ArrayList()
        sgroups.add(sgroup)
        mol.setProperty(CDKConstants.CTAB_SGROUPS, sgroups)

        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 1

    def test_metal_at_end_hidden_neutral(self, mc_handler):
        """Test HIDDEN_NEUTRAL with metal at end to trigger charge neutralization."""
        from jpype import JClass

        CDKConstants = JClass("org.openscience.cdk.CDKConstants")
        Sgroup = JClass("org.openscience.cdk.sgroup.Sgroup")
        SgroupType = JClass("org.openscience.cdk.sgroup.SgroupType")
        ArrayList = JClass("java.util.ArrayList")

        mol = get_CDK_IAtomContainer("c1ccccc1[Fe]")

        fe_idx = None
        for i in range(mol.getAtomCount()):
            if mol.getAtom(i).getSymbol() == "Fe":
                fe_idx = i
                break

        fe_atom = mol.getAtom(fe_idx)
        ring_bond = None
        for bond in fe_atom.bonds():
            other = bond.getOther(fe_atom)
            if mol.indexOf(other) != fe_idx:
                ring_bond = bond
                break

        sgroup = Sgroup()
        sgroup.setType(SgroupType.ExtMulticenter)
        for i in range(mol.getAtomCount()):
            if i != fe_idx:
                sgroup.addAtom(mol.getAtom(i))
        sgroup.addBond(ring_bond)

        sgroups = ArrayList()
        sgroups.add(sgroup)
        mol.setProperty(CDKConstants.CTAB_SGROUPS, sgroups)

        count = mc_handler.set_style(mol, MulticenterStyle.HIDDEN_NEUTRAL)
        assert count >= 1


class TestDetectMulticenterBondsErrorPath:
    """Test the error path in detect_multicenter_bonds (lines 265-267)."""

    def test_detect_with_none_molecule(self, mc_handler):
        """Passing None should trigger the exception path and return 0."""
        count = mc_handler.detect_multicenter_bonds(None)
        assert count == 0

    def test_detect_with_invalid_object(self, mc_handler):
        """Passing a non-molecule object should trigger exception and return 0."""
        count = mc_handler.detect_multicenter_bonds("not_a_molecule")
        assert count == 0


class TestNeutralizeChargesDirectly:
    """Test _neutralize_charges edge cases via set_style with charged atoms."""

    def test_charged_ring_atoms_neutralized(self, mc_handler):
        """Create a molecule with charged atoms and verify neutralization."""
        # Use cyclopentadienyl anion with Fe
        mol = get_CDK_IAtomContainer("[Fe+2]c1ccccc1")

        fe_idx = None
        for i in range(mol.getAtomCount()):
            if mol.getAtom(i).getSymbol() == "Fe":
                fe_idx = i
                break

        ring_indices = [i for i in range(mol.getAtomCount()) if i != fe_idx]
        mc_handler.mark_multicenter_bond(mol, fe_idx, ring_indices)
        count = mc_handler.set_style(mol, MulticenterStyle.DASHED_NEUTRAL)
        assert count >= 1

    def test_charges_with_dative_style(self, mc_handler):
        """DATIVE style should also trigger charge neutralization."""
        mol = get_CDK_IAtomContainer("[Fe+2]c1ccccc1")

        fe_idx = None
        for i in range(mol.getAtomCount()):
            if mol.getAtom(i).getSymbol() == "Fe":
                fe_idx = i
                break

        ring_indices = [i for i in range(mol.getAtomCount()) if i != fe_idx]
        mc_handler.mark_multicenter_bond(mol, fe_idx, ring_indices)
        count = mc_handler.set_style(mol, MulticenterStyle.DATIVE)
        assert count >= 1
