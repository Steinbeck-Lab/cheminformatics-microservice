from __future__ import annotations

import pytest

from app.modules.cdk_depict.annotations import (
    AnnotationSystem,
    AnnotationMode,
    MAPPING_COLORS,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from jpype import JClass


@pytest.fixture
def annotation_system():
    return AnnotationSystem()


@pytest.fixture
def simple_molecule():
    return get_CDK_IAtomContainer("CCO")


@pytest.fixture
def chiral_molecule():
    return get_CDK_IAtomContainer("C[C@H](O)CC")


@pytest.fixture
def complex_molecule():
    return get_CDK_IAtomContainer("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")


class TestAnnotationSystemInitialization:
    """Test annotation system initialization."""

    def test_default_initialization(self):
        system = AnnotationSystem()
        assert system.cdk_base == "org.openscience.cdk"

    def test_annotation_system_has_required_attributes(self):
        system = AnnotationSystem()
        assert hasattr(system, "StandardGenerator")
        assert hasattr(system, "Color")
        assert hasattr(system, "CDKConstants")


class TestAnnotationModes:
    """Test different annotation modes."""

    def test_none_mode_no_annotations(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.NONE
        )
        assert gen is not None

    def test_number_mode(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.NUMBER
        )
        assert gen is not None

    def test_bondnumber_mode(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.BONDNUMBER
        )
        assert gen is not None

    def test_mapidx_mode(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.MAPIDX
        )
        assert gen is not None

    def test_atomvalue_mode(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.ATOMVALUE
        )
        assert gen is not None

    def test_colmap_mode(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.COLMAP
        )
        assert gen is not None

    def test_cip_mode(self, annotation_system, chiral_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, chiral_molecule, AnnotationMode.CIP
        )
        assert gen is not None


class TestAtomNumbering:
    """Test atom numbering annotations."""

    def test_apply_atom_numbers(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.NUMBER
        )
        assert gen is not None

    def test_atom_numbers_on_complex_molecule(
        self, annotation_system, complex_molecule
    ):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, complex_molecule, AnnotationMode.NUMBER
        )
        assert gen is not None

    def test_atom_numbers_with_aromatic(self, annotation_system):
        mol = get_CDK_IAtomContainer("c1ccccc1")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.NUMBER
        )
        assert gen is not None


class TestBondNumbering:
    """Test bond numbering annotations."""

    def test_apply_bond_numbers(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.BONDNUMBER
        )
        assert gen is not None

    def test_bond_numbers_on_complex_molecule(
        self, annotation_system, complex_molecule
    ):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, complex_molecule, AnnotationMode.BONDNUMBER
        )
        assert gen is not None


class TestMappingAnnotations:
    """Test atom mapping annotations."""

    def test_apply_map_indices(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.MAPIDX
        )
        assert gen is not None

    def test_colored_mapping(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.COLMAP
        )
        assert gen is not None

    def test_mapping_with_reaction_molecule(self, annotation_system):
        mol = get_CDK_IAtomContainer("C([H:1])([H:2])([H:3])[H:4]")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.MAPIDX
        )
        assert gen is not None


class TestAtomValueAnnotations:
    """Test atom value annotations."""

    def test_apply_atom_values(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.ATOMVALUE
        )
        assert gen is not None

    def test_atom_values_on_complex_molecule(self, annotation_system, complex_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, complex_molecule, AnnotationMode.ATOMVALUE
        )
        assert gen is not None


class TestCIPAnnotations:
    """Test CIP stereochemistry annotations."""

    def test_apply_cip_on_chiral_molecule(self, annotation_system, chiral_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, chiral_molecule, AnnotationMode.CIP
        )
        assert gen is not None

    def test_cip_on_simple_molecule(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.CIP
        )
        assert gen is not None

    def test_cip_on_ez_stereochemistry(self, annotation_system):
        mol = get_CDK_IAtomContainer("C/C=C/C")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.CIP
        )
        assert gen is not None

    def test_cip_on_multiple_chiral_centers(self, annotation_system):
        mol = get_CDK_IAtomContainer("C[C@H](O)[C@@H](C)Cl")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.CIP
        )
        assert gen is not None


class TestMappingColors:
    """Test mapping color constants."""

    def test_mapping_colors_defined(self):
        assert isinstance(MAPPING_COLORS, list)
        assert len(MAPPING_COLORS) > 0

    def test_mapping_colors_contain_indices(self):
        # MAPPING_COLORS is a list where index represents the mapping number
        for idx in range(1, 10):  # Start from 1 since 0 is None
            assert idx < len(MAPPING_COLORS)
            assert MAPPING_COLORS[idx] is not None

    def test_mapping_colors_are_tuples(self):
        # Each color (except index 0) should be an RGB tuple
        for idx in range(1, len(MAPPING_COLORS)):
            color = MAPPING_COLORS[idx]
            assert isinstance(color, tuple)
            assert len(color) == 3  # RGB
            # Verify values are in valid range
            for component in color:
                assert 0 <= component <= 255


class TestApplyAnnotationsIntegration:
    """Test apply_annotations method integration."""

    def test_apply_number_annotation(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.NUMBER
        )
        assert gen is not None

    def test_apply_bond_annotation(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.BONDNUMBER
        )
        assert gen is not None

    def test_apply_cip_annotation(self, annotation_system, chiral_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, chiral_molecule, AnnotationMode.CIP
        )
        assert gen is not None

    def test_apply_mapidx_annotation(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.MAPIDX
        )
        assert gen is not None

    def test_apply_colmap_annotation(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.COLMAP
        )
        assert gen is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_atom_molecule(self, annotation_system):
        mol = get_CDK_IAtomContainer("C")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.NUMBER
        )
        assert gen is not None

    def test_disconnected_fragments(self, annotation_system):
        mol = get_CDK_IAtomContainer("CCO.C1CCOC1")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.NUMBER
        )
        assert gen is not None

    def test_aromatic_molecule(self, annotation_system):
        mol = get_CDK_IAtomContainer("c1ccccc1")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.NUMBER
        )
        assert gen is not None

    def test_charged_molecule(self, annotation_system):
        mol = get_CDK_IAtomContainer("[NH4+]")
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = annotation_system.apply_annotations(
            DepictionGenerator, mol, AnnotationMode.NUMBER
        )
        assert gen is not None


class TestAnnotationPersistence:
    """Test that annotations can be applied and changed."""

    def test_multiple_annotation_modes(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen1 = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.NUMBER
        )
        gen2 = annotation_system.apply_annotations(
            gen1, simple_molecule, AnnotationMode.BONDNUMBER
        )
        assert gen2 is not None

    def test_annotation_change(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen1 = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.NUMBER
        )
        gen2 = annotation_system.apply_annotations(
            gen1, simple_molecule, AnnotationMode.CIP
        )
        assert gen2 is not None

    def test_remove_annotations(self, annotation_system, simple_molecule):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen1 = annotation_system.apply_annotations(
            DepictionGenerator, simple_molecule, AnnotationMode.NUMBER
        )
        gen2 = annotation_system.apply_annotations(
            gen1, simple_molecule, AnnotationMode.NONE
        )
        assert gen2 is not None


class TestGetAnnotationMode:
    """Test the get_annotation_mode helper function."""

    def test_valid_mode_number(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("number")
        assert mode == AnnotationMode.NUMBER

    def test_valid_mode_bondnumber(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("bondnumber")
        assert mode == AnnotationMode.BONDNUMBER

    def test_valid_mode_mapidx(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("mapidx")
        assert mode == AnnotationMode.MAPIDX

    def test_valid_mode_atomvalue(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("atomvalue")
        assert mode == AnnotationMode.ATOMVALUE

    def test_valid_mode_colmap(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("colmap")
        assert mode == AnnotationMode.COLMAP

    def test_valid_mode_cip(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("cip")
        assert mode == AnnotationMode.CIP

    def test_valid_mode_none(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("none")
        assert mode == AnnotationMode.NONE

    def test_valid_mode_with_whitespace(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("  number  ")
        assert mode == AnnotationMode.NUMBER

    def test_valid_mode_case_insensitive(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        mode = get_annotation_mode("NUMBER")
        assert mode == AnnotationMode.NUMBER

    def test_invalid_mode_raises_valueerror(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        with pytest.raises(ValueError) as exc_info:
            get_annotation_mode("invalid_mode")
        assert "Invalid annotation mode" in str(exc_info.value)
        assert "valid_mode" not in str(exc_info.value).lower() or "Valid modes" in str(
            exc_info.value
        )

    def test_invalid_mode_lists_valid_modes(self):
        from app.modules.cdk_depict.annotations import get_annotation_mode

        with pytest.raises(ValueError) as exc_info:
            get_annotation_mode("xyz")
        error_msg = str(exc_info.value)
        assert "none" in error_msg
        assert "number" in error_msg
        assert "cip" in error_msg


class TestApplyAnnotationsConvenience:
    """Test the module-level apply_annotations convenience function."""

    def test_apply_annotations_number_mode(self, simple_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, simple_molecule, mode="number")
        assert result is not None

    def test_apply_annotations_none_mode(self, simple_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, simple_molecule, mode="none")
        assert result is not None

    def test_apply_annotations_cip_mode(self, chiral_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, chiral_molecule, mode="cip")
        assert result is not None

    def test_apply_annotations_invalid_mode_returns_generator(self, simple_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, simple_molecule, mode="invalid_mode")
        # Should return unmodified generator on error
        assert result is not None

    def test_apply_annotations_mapidx_mode(self, simple_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, simple_molecule, mode="mapidx")
        assert result is not None

    def test_apply_annotations_bondnumber_mode(self, simple_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, simple_molecule, mode="bondnumber")
        assert result is not None

    def test_apply_annotations_atomvalue_mode(self, simple_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, simple_molecule, mode="atomvalue")
        assert result is not None

    def test_apply_annotations_colmap_mode(self, simple_molecule):
        from app.modules.cdk_depict.annotations import (
            apply_annotations as apply_annotations_func,
        )

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = apply_annotations_func(gen, simple_molecule, mode="colmap")
        assert result is not None


class TestAtomValueAnnotationsWithValues:
    """Test atom value annotations when atoms have actual values set."""

    def test_atom_values_with_comment_property(self, annotation_system):
        """Test annotating atoms that have COMMENT property set."""
        mol = get_CDK_IAtomContainer("CCO")
        CDKConstants = JClass("org.openscience.cdk.CDKConstants")
        # Set COMMENT property on atoms
        for i, atom in enumerate(mol.atoms()):
            atom.setProperty(CDKConstants.COMMENT, f"val{i}")

        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(gen, mol, AnnotationMode.ATOMVALUE)
        assert result is not None

    def test_atom_values_without_comment_property(self, annotation_system):
        """Test annotating atoms that have no COMMENT property."""
        mol = get_CDK_IAtomContainer("CCO")
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(gen, mol, AnnotationMode.ATOMVALUE)
        assert result is not None


class TestColorMappingWithMappedAtoms:
    """Test color mapping when atoms actually have mapping numbers."""

    def test_color_mapping_with_mapped_atoms(self, annotation_system):
        """Test colmap mode with atoms that have atom-atom mapping."""
        mol = get_CDK_IAtomContainer("[CH3:1][OH:2]")
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(gen, mol, AnnotationMode.COLMAP)
        assert result is not None

    def test_color_mapping_with_various_map_indices(self, annotation_system):
        """Test colmap mode with various mapping indices within MAPPING_COLORS range."""
        mol = get_CDK_IAtomContainer("[C:1]([H:2])([H:3])[O:4][H:5]")
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(gen, mol, AnnotationMode.COLMAP)
        assert result is not None


class TestAnnotationExceptionPaths:
    """Test error handling paths in annotation system."""

    def test_apply_annotations_with_error_returns_generator(self, annotation_system):
        """Test that apply_annotations returns generator even on internal error."""
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        # Pass None molecule - should trigger exception and return generator
        result = annotation_system.apply_annotations(
            gen, None, AnnotationMode.ATOMVALUE
        )
        assert result is not None

    def test_annotate_atom_values_exception_path(self, annotation_system):
        """Test _annotate_atom_values exception handling."""
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(
            gen, None, AnnotationMode.ATOMVALUE
        )
        assert result is not None

    def test_annotate_color_mapping_exception_path(self, annotation_system):
        """Test _annotate_color_mapping exception handling."""
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(gen, None, AnnotationMode.COLMAP)
        assert result is not None

    def test_annotate_cip_exception_path(self, annotation_system):
        """Test _annotate_cip exception handling."""
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(gen, None, AnnotationMode.CIP)
        assert result is not None

    def test_is_reaction_with_non_reaction_object(self, annotation_system):
        """Test _is_reaction with a non-reaction object."""
        mol = get_CDK_IAtomContainer("CCO")
        assert annotation_system._is_reaction(mol) is False

    def test_is_reaction_with_none(self, annotation_system):
        """Test _is_reaction with None."""
        assert annotation_system._is_reaction(None) is False


class TestReactionAnnotations:
    """Test annotation modes on reaction objects."""

    def _create_reaction(self):
        """Helper to create a simple CDK reaction."""
        Reaction = JClass("org.openscience.cdk.Reaction")
        reaction = Reaction()
        reactant = get_CDK_IAtomContainer("CCO")
        product = get_CDK_IAtomContainer("CC=O")
        reaction.addReactant(reactant)
        reaction.addProduct(product)
        return reaction

    def test_atomvalue_on_reaction(self, annotation_system):
        """Test atom value annotation on reaction objects."""
        reaction = self._create_reaction()
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(
            gen, reaction, AnnotationMode.ATOMVALUE
        )
        assert result is not None

    def test_colmap_on_reaction(self, annotation_system):
        """Test color mapping annotation on reaction objects."""
        reaction = self._create_reaction()
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(
            gen, reaction, AnnotationMode.COLMAP, is_reaction=True
        )
        assert result is not None

    def test_cip_on_reaction(self, annotation_system):
        """Test CIP annotation on reaction objects."""
        reaction = self._create_reaction()
        gen = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        result = annotation_system.apply_annotations(
            gen, reaction, AnnotationMode.CIP, is_reaction=True
        )
        assert result is not None

    def test_is_reaction_with_actual_reaction(self, annotation_system):
        """Test _is_reaction correctly identifies IReaction objects."""
        reaction = self._create_reaction()
        assert annotation_system._is_reaction(reaction) is True
