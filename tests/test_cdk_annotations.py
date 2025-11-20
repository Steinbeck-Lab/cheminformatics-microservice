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
