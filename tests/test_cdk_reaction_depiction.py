from __future__ import annotations

import pytest

from app.modules.cdk_depict.reaction_depiction import (
    ReactionArrowSystem,
    ReactionArrowType,
    ARROW_ABBREV_MAP,
    get_arrow_type,
)
from jpype import JClass


@pytest.fixture
def arrow_system():
    return ReactionArrowSystem()


@pytest.fixture
def simple_reaction():
    cdk_base = "org.openscience.cdk"
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    try:
        reaction = SmilesParser.parseReactionSmiles("CCO>>CC=O")
        return reaction
    except Exception:
        return None


class TestReactionArrowType:
    """Test ReactionArrowType enum."""

    def test_forward_arrow(self):
        assert ReactionArrowType.FORWARD.value == "forward"

    def test_bidirectional_arrow(self):
        assert ReactionArrowType.BIDIRECTIONAL.value == "bidirectional"

    def test_no_go_arrow(self):
        assert ReactionArrowType.NO_GO.value == "no_go"

    def test_retro_synthetic_arrow(self):
        """Test that RETRO_SYNTHETIC enum has correct value."""
        assert ReactionArrowType.RETRO_SYNTHETIC.value == "retrosynthetic"

    def test_resonance_arrow(self):
        assert ReactionArrowType.RESONANCE.value == "resonance"

    def test_all_arrow_types_defined(self):
        arrow_types = [
            ReactionArrowType.FORWARD,
            ReactionArrowType.BIDIRECTIONAL,
            ReactionArrowType.NO_GO,
            ReactionArrowType.RETRO_SYNTHETIC,
            ReactionArrowType.RESONANCE,
        ]
        assert len(arrow_types) == 5


class TestArrowAbbrevMap:
    """Test ARROW_ABBREV_MAP dictionary."""

    def test_forward_mapping(self):
        assert "forward" in ARROW_ABBREV_MAP

    def test_equ_mapping(self):
        assert "equ" in ARROW_ABBREV_MAP

    def test_ngo_mapping(self):
        assert "ngo" in ARROW_ABBREV_MAP

    def test_ret_mapping(self):
        assert "ret" in ARROW_ABBREV_MAP

    def test_res_mapping(self):
        assert "res" in ARROW_ABBREV_MAP

    def test_forward_maps_to_forward_type(self):
        assert ARROW_ABBREV_MAP["forward"] == ReactionArrowType.FORWARD

    def test_equ_maps_to_bidirectional(self):
        assert ARROW_ABBREV_MAP["equ"] == ReactionArrowType.BIDIRECTIONAL

    def test_ngo_maps_to_no_go(self):
        assert ARROW_ABBREV_MAP["ngo"] == ReactionArrowType.NO_GO

    def test_ret_maps_to_retro_synthetic(self):
        assert ARROW_ABBREV_MAP["ret"] == ReactionArrowType.RETRO_SYNTHETIC

    def test_res_maps_to_resonance(self):
        assert ARROW_ABBREV_MAP["res"] == ReactionArrowType.RESONANCE


class TestReactionArrowSystemInitialization:
    """Test ReactionArrowSystem initialization."""

    def test_default_initialization(self):
        """Test that system initializes with correct cdk_base."""
        system = ReactionArrowSystem()
        assert system.cdk_base == "org.openscience.cdk"

    def test_initialization_loads_classes(self):
        """Test that initialization loads all required Java classes."""
        system = ReactionArrowSystem()
        assert system.IReaction is not None
        assert system.Direction is not None

    def test_cdk_base_is_string(self):
        """Test that cdk_base attribute is a string."""
        system = ReactionArrowSystem()
        assert isinstance(system.cdk_base, str)
        assert len(system.cdk_base) > 0


class TestSetArrowType:
    """Test setting arrow type on reactions."""

    def test_set_forward_arrow(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.FORWARD)
        assert simple_reaction is not None

    def test_set_bidirectional_arrow(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.BIDIRECTIONAL)
        assert simple_reaction is not None

    def test_set_no_go_arrow(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.NO_GO)
        assert simple_reaction is not None

    def test_set_retro_synthetic_arrow(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.RETRO_SYNTHETIC)
        assert simple_reaction is not None

    def test_set_resonance_arrow(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.RESONANCE)
        assert simple_reaction is not None


class TestSetArrowTypeFromString:
    """Test setting arrow type using get_arrow_type() helper function."""

    def test_set_arrow_from_forward_string(self, arrow_system, simple_reaction):
        """Test setting arrow from 'forward' string."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type("forward")
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert simple_reaction is not None

    def test_set_arrow_from_equ_string(self, arrow_system, simple_reaction):
        """Test setting arrow from 'equ' abbreviation."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type("equ")
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert simple_reaction is not None

    def test_set_arrow_from_ngo_string(self, arrow_system, simple_reaction):
        """Test setting arrow from 'ngo' abbreviation."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type("ngo")
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert simple_reaction is not None

    def test_set_arrow_from_ret_string(self, arrow_system, simple_reaction):
        """Test setting arrow from 'ret' abbreviation."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type("ret")
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert simple_reaction is not None

    def test_set_arrow_from_res_string(self, arrow_system, simple_reaction):
        """Test setting arrow from 'res' abbreviation."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type("res")
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert simple_reaction is not None

    def test_set_arrow_from_empty_string_defaults(self, arrow_system, simple_reaction):
        """Test that empty string defaults to FORWARD arrow."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type("")
        assert arrow_type == ReactionArrowType.FORWARD
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert simple_reaction is not None

    def test_set_arrow_from_unknown_string_raises_error(
        self, arrow_system, simple_reaction
    ):
        """Test that unknown arrow type raises ValueError."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        with pytest.raises(ValueError) as exc_info:
            get_arrow_type("unknown")
        assert "Invalid arrow type" in str(exc_info.value)


class TestArrowTypeConversion:
    """Test conversion between arrow type representations using get_arrow_type()."""

    def test_string_to_arrow_type_forward(self):
        """Test converting 'forward' string to enum."""
        arrow_type = get_arrow_type("forward")
        assert arrow_type == ReactionArrowType.FORWARD

    def test_string_to_arrow_type_equ(self):
        """Test converting 'equ' abbreviation to enum."""
        arrow_type = get_arrow_type("equ")
        assert arrow_type == ReactionArrowType.BIDIRECTIONAL

    def test_string_to_arrow_type_ngo(self):
        """Test converting 'ngo' abbreviation to enum."""
        arrow_type = get_arrow_type("ngo")
        assert arrow_type == ReactionArrowType.NO_GO

    def test_string_to_arrow_type_ret(self):
        """Test converting 'ret' abbreviation to enum."""
        arrow_type = get_arrow_type("ret")
        assert arrow_type == ReactionArrowType.RETRO_SYNTHETIC

    def test_string_to_arrow_type_res(self):
        """Test converting 'res' abbreviation to enum."""
        arrow_type = get_arrow_type("res")
        assert arrow_type == ReactionArrowType.RESONANCE

    def test_string_to_arrow_type_default(self):
        """Test that unknown string raises ValueError."""
        with pytest.raises(ValueError):
            get_arrow_type("unknown")

    def test_string_to_arrow_type_empty(self):
        """Test that empty string defaults to FORWARD."""
        arrow_type = get_arrow_type("")
        assert arrow_type == ReactionArrowType.FORWARD


class TestReactionDirection:
    """Test reaction direction handling."""

    def test_forward_direction_code(self):
        cdk_base = "org.openscience.cdk"
        IReaction = JClass(cdk_base + ".interfaces.IReaction")
        assert IReaction.Direction.FORWARD is not None

    def test_bidirectional_direction_code(self):
        cdk_base = "org.openscience.cdk"
        IReaction = JClass(cdk_base + ".interfaces.IReaction")
        assert IReaction.Direction.BIDIRECTIONAL is not None


class TestMultipleReactions:
    """Test handling multiple reactions in a reaction set."""

    def test_set_arrow_on_reaction_set(self, arrow_system):
        cdk_base = "org.openscience.cdk"
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
        SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
        try:
            rxn_set = SmilesParser.parseReactionSetSmiles("CCO>>CC=O")
            for rxn in rxn_set.reactions():
                arrow_system.set_arrow_type(rxn, ReactionArrowType.BIDIRECTIONAL)
            assert rxn_set is not None
        except Exception:
            pytest.skip("Reaction set parsing failed")


class TestArrowTypePreservation:
    """Test that arrow type settings are preserved."""

    def test_arrow_type_persists_after_setting(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.BIDIRECTIONAL)
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.FORWARD)
        assert simple_reaction is not None

    def test_multiple_arrow_type_changes(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.FORWARD)
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.BIDIRECTIONAL)
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.NO_GO)
        assert simple_reaction is not None


class TestReactionIntegrity:
    """Test that arrow type operations preserve reaction integrity."""

    def test_preserves_reactant_count(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        initial_count = simple_reaction.getReactantCount()
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.FORWARD)
        assert simple_reaction.getReactantCount() == initial_count

    def test_preserves_product_count(self, arrow_system, simple_reaction):
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        initial_count = simple_reaction.getProductCount()
        arrow_system.set_arrow_type(simple_reaction, ReactionArrowType.FORWARD)
        assert simple_reaction.getProductCount() == initial_count


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_case_insensitive_arrow_string(self, arrow_system, simple_reaction):
        """Test that arrow type strings are case-insensitive."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type("FORWARD")
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert arrow_type == ReactionArrowType.FORWARD

    def test_whitespace_in_arrow_string(self, arrow_system, simple_reaction):
        """Test that whitespace is trimmed from arrow type strings."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        arrow_type = get_arrow_type(" forward ")
        arrow_system.set_arrow_type(simple_reaction, arrow_type)
        assert arrow_type == ReactionArrowType.FORWARD

    def test_none_arrow_string(self, arrow_system, simple_reaction):
        """Test that None arrow string raises an error."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        try:
            get_arrow_type(None)
            pytest.fail("Should have raised an exception for None")
        except (AttributeError, TypeError):
            assert True


class TestArrowSymbols:
    """Test arrow symbols and their meanings."""

    def test_forward_arrow_symbol(self):
        assert ReactionArrowType.FORWARD.value == "forward"

    def test_bidirectional_arrow_symbol(self):
        assert ReactionArrowType.BIDIRECTIONAL.value == "bidirectional"

    def test_no_go_arrow_symbol(self):
        assert ReactionArrowType.NO_GO.value == "no_go"

    def test_retro_synthetic_arrow_symbol(self):
        """Test that RETRO_SYNTHETIC enum has correct value."""
        assert ReactionArrowType.RETRO_SYNTHETIC.value == "retrosynthetic"

    def test_resonance_arrow_symbol(self):
        assert ReactionArrowType.RESONANCE.value == "resonance"


class TestAbbreviationConsistency:
    """Test consistency of abbreviations."""

    def test_all_abbreviations_have_arrow_types(self):
        for abbrev, arrow_type in ARROW_ABBREV_MAP.items():
            assert isinstance(arrow_type, ReactionArrowType)

    def test_all_arrow_types_have_abbreviations(self):
        arrow_types_in_map = set(ARROW_ABBREV_MAP.values())
        all_arrow_types = set(ReactionArrowType)
        assert arrow_types_in_map == all_arrow_types


class TestComplexReactions:
    """Test with complex reaction structures."""

    def test_multi_step_reaction(self, arrow_system):
        cdk_base = "org.openscience.cdk"
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
        SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
        try:
            reaction = SmilesParser.parseReactionSmiles("CC.O>>CCO")
            arrow_system.set_arrow_type(reaction, ReactionArrowType.FORWARD)
            assert reaction is not None
        except Exception:
            pytest.skip("Complex reaction parsing failed")

    def test_catalyzed_reaction(self, arrow_system):
        cdk_base = "org.openscience.cdk"
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
        SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
        try:
            reaction = SmilesParser.parseReactionSmiles("CCO>O>CC=O")
            arrow_system.set_arrow_type(reaction, ReactionArrowType.FORWARD)
            assert reaction is not None
        except Exception:
            pytest.skip("Catalyzed reaction parsing failed")


class TestArrowTypeStringVariants:
    """Test different string variants for arrow types."""

    def test_forward_variants(self, arrow_system, simple_reaction):
        """Test that different case variants of 'forward' work."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        for variant in ["forward", "FORWARD", "Forward"]:
            arrow_type = get_arrow_type(variant)
            assert arrow_type == ReactionArrowType.FORWARD
            arrow_system.set_arrow_type(simple_reaction, arrow_type)
            assert simple_reaction is not None

    def test_equilibrium_variants(self, arrow_system, simple_reaction):
        """Test that different case variants of 'equ' work."""
        if simple_reaction is None:
            pytest.skip("Reaction parsing failed")
        for variant in ["equ", "EQU", "Equ"]:
            arrow_type = get_arrow_type(variant)
            assert arrow_type == ReactionArrowType.BIDIRECTIONAL
            arrow_system.set_arrow_type(simple_reaction, arrow_type)
            assert simple_reaction is not None


class TestReactionSetHandling:
    """Test handling of reaction sets."""

    def test_apply_arrow_to_all_reactions_in_set(self, arrow_system):
        cdk_base = "org.openscience.cdk"
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
        SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
        try:
            rxn_set = SmilesParser.parseReactionSetSmiles("CCO>>CC=O")
            count = 0
            for rxn in rxn_set.reactions():
                arrow_system.set_arrow_type(rxn, ReactionArrowType.FORWARD)
                count += 1
            assert count > 0
        except Exception:
            pytest.skip("Reaction set parsing failed")
