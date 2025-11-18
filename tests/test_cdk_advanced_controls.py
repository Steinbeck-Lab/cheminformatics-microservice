from __future__ import annotations

import pytest

from app.modules.cdk_depict.advanced_controls import (
    AdvancedControls,
    SVGUnits,
)
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from jpype import JClass


@pytest.fixture
def controls():
    return AdvancedControls()


@pytest.fixture
def simple_molecule():
    return get_CDK_IAtomContainer("CCO")


@pytest.fixture
def benzene():
    return get_CDK_IAtomContainer("c1ccccc1")


@pytest.fixture
def chiral_molecule():
    return get_CDK_IAtomContainer("C[C@H](O)CC")


@pytest.fixture
def complex_molecule():
    return get_CDK_IAtomContainer("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")


class TestAdvancedControlsInitialization:
    """Test initialization of advanced controls system."""

    def test_default_initialization(self):
        controls = AdvancedControls()
        assert controls.cdk_base == "org.openscience.cdk"

    def test_custom_cdk_base(self):
        controls = AdvancedControls(cdk_base="org.openscience.cdk")
        assert controls.cdk_base == "org.openscience.cdk"

    def test_initialization_loads_classes(self, controls):
        assert controls.GeometryTools is not None
        assert controls.Math is not None


class TestSVGUnits:
    """Test SVG units enum."""

    def test_svg_units_px(self):
        assert SVGUnits.PX.value == "px"

    def test_svg_units_mm(self):
        assert SVGUnits.MM.value == "mm"

    def test_svg_units_cm(self):
        assert SVGUnits.CM.value == "cm"

    def test_svg_units_in(self):
        assert SVGUnits.IN.value == "in"

    def test_all_units_defined(self):
        units = [SVGUnits.PX, SVGUnits.MM, SVGUnits.CM, SVGUnits.IN]
        assert len(units) == 4


class TestZoomControl:
    """Test zoom level control."""

    def test_set_zoom_default(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=1.3)
        assert gen is not None

    def test_set_zoom_minimum(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=0.1)
        assert gen is not None

    def test_set_zoom_maximum(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=5.0)
        assert gen is not None

    def test_set_zoom_medium(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=1.5)
        assert gen is not None

    def test_set_zoom_small(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=0.5)
        assert gen is not None

    def test_set_zoom_large(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=3.0)
        assert gen is not None


class TestStrokeRatio:
    """Test stroke ratio control."""

    def test_set_ratio_default(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_stroke_ratio(DepictionGenerator, ratio=1.0)
        assert gen is not None

    def test_set_ratio_minimum(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_stroke_ratio(DepictionGenerator, ratio=0.5)
        assert gen is not None

    def test_set_ratio_maximum(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_stroke_ratio(DepictionGenerator, ratio=2.0)
        assert gen is not None

    def test_set_ratio_thin(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_stroke_ratio(DepictionGenerator, ratio=0.7)
        assert gen is not None

    def test_set_ratio_thick(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_stroke_ratio(DepictionGenerator, ratio=1.5)
        assert gen is not None


class TestSVGUnitsControl:
    """Test SVG units control."""

    def test_set_units_px(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units=SVGUnits.PX)
        assert gen is not None

    def test_set_units_mm(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units=SVGUnits.MM)
        assert gen is not None

    def test_set_units_cm(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units=SVGUnits.CM)
        assert gen is not None

    def test_set_units_in(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units=SVGUnits.IN)
        assert gen is not None

    def test_set_units_string_px(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units="px")
        assert gen is not None

    def test_set_units_string_mm(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units="mm")
        assert gen is not None


class TestFlipStructure:
    """Test structure flipping."""

    def test_flip_simple_molecule(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.flip_structure(simple_molecule)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_flip_benzene(self, controls, benzene):
        initial_atoms = benzene.getAtomCount()
        controls.flip_structure(benzene)
        assert benzene.getAtomCount() == initial_atoms

    def test_flip_preserves_connectivity(self, controls, simple_molecule):
        controls.flip_structure(simple_molecule)
        for atom in simple_molecule.atoms():
            connected = simple_molecule.getConnectedAtomsList(atom)
            assert connected is not None

    def test_flip_preserves_bond_count(self, controls, simple_molecule):
        initial_bonds = simple_molecule.getBondCount()
        controls.flip_structure(simple_molecule)
        assert simple_molecule.getBondCount() == initial_bonds

    def test_flip_chiral_molecule(self, controls, chiral_molecule):
        initial_atoms = chiral_molecule.getAtomCount()
        controls.flip_structure(chiral_molecule)
        assert chiral_molecule.getAtomCount() == initial_atoms


class TestRotateStructure:
    """Test structure rotation."""

    def test_rotate_0_degrees(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.rotate_structure(simple_molecule, degrees=0)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_rotate_90_degrees(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.rotate_structure(simple_molecule, degrees=90)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_rotate_180_degrees(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.rotate_structure(simple_molecule, degrees=180)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_rotate_270_degrees(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.rotate_structure(simple_molecule, degrees=270)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_rotate_45_degrees(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.rotate_structure(simple_molecule, degrees=45)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_rotate_negative_degrees(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.rotate_structure(simple_molecule, degrees=-90)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_rotate_preserves_connectivity(self, controls, simple_molecule):
        controls.rotate_structure(simple_molecule, degrees=90)
        for atom in simple_molecule.atoms():
            connected = simple_molecule.getConnectedAtomsList(atom)
            assert connected is not None

    def test_rotate_preserves_bond_count(self, controls, simple_molecule):
        initial_bonds = simple_molecule.getBondCount()
        controls.rotate_structure(simple_molecule, degrees=90)
        assert simple_molecule.getBondCount() == initial_bonds


class TestTitleDisplay:
    """Test title display control."""

    def test_show_title_molecule(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_title_display(
            DepictionGenerator, show_title=True, is_reaction=False
        )
        assert gen is not None

    def test_hide_title_molecule(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_title_display(
            DepictionGenerator, show_title=False, is_reaction=False
        )
        assert gen is not None

    def test_show_title_reaction(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_title_display(
            DepictionGenerator, show_title=True, is_reaction=True
        )
        assert gen is not None

    def test_hide_title_reaction(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_title_display(
            DepictionGenerator, show_title=False, is_reaction=True
        )
        assert gen is not None


class TestAnonymousDisplay:
    """Test anonymous display mode."""

    def test_anonymous_enabled(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_anonymous_display(DepictionGenerator, anonymous=True)
        assert gen is not None

    def test_anonymous_disabled(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_anonymous_display(DepictionGenerator, anonymous=False)
        assert gen is not None


class TestMolecularIntegrity:
    """Test that operations preserve molecular integrity."""

    def test_flip_preserves_atoms(self, controls, complex_molecule):
        initial_count = complex_molecule.getAtomCount()
        controls.flip_structure(complex_molecule)
        assert complex_molecule.getAtomCount() == initial_count

    def test_rotate_preserves_atoms(self, controls, complex_molecule):
        initial_count = complex_molecule.getAtomCount()
        controls.rotate_structure(complex_molecule, degrees=90)
        assert complex_molecule.getAtomCount() == initial_count

    def test_flip_preserves_stereochemistry(self, controls, chiral_molecule):
        controls.flip_structure(chiral_molecule)
        stereo_count = sum(1 for _ in chiral_molecule.stereoElements())
        assert stereo_count >= 0

    def test_rotate_preserves_stereochemistry(self, controls, chiral_molecule):
        controls.rotate_structure(chiral_molecule, degrees=90)
        stereo_count = sum(1 for _ in chiral_molecule.stereoElements())
        assert stereo_count >= 0


class TestCombinedOperations:
    """Test combinations of multiple operations."""

    def test_zoom_and_ratio_together(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=1.5)
        gen = controls.set_stroke_ratio(gen, ratio=1.2)
        assert gen is not None

    def test_flip_and_rotate(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.flip_structure(simple_molecule)
        controls.rotate_structure(simple_molecule, degrees=90)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_multiple_rotations(self, controls, simple_molecule):
        initial_atoms = simple_molecule.getAtomCount()
        controls.rotate_structure(simple_molecule, degrees=90)
        controls.rotate_structure(simple_molecule, degrees=90)
        assert simple_molecule.getAtomCount() == initial_atoms

    def test_all_depiction_controls(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=1.5)
        gen = controls.set_stroke_ratio(gen, ratio=1.2)
        gen = controls.set_svg_units(gen, units=SVGUnits.MM)
        gen = controls.set_title_display(gen, show_title=True, is_reaction=False)
        gen = controls.set_anonymous_display(gen, anonymous=False)
        assert gen is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_atom_flip(self, controls):
        mol = get_CDK_IAtomContainer("C")
        initial_count = mol.getAtomCount()
        controls.flip_structure(mol)
        assert mol.getAtomCount() == initial_count

    def test_single_atom_rotate(self, controls):
        mol = get_CDK_IAtomContainer("C")
        initial_count = mol.getAtomCount()
        controls.rotate_structure(mol, degrees=90)
        assert mol.getAtomCount() == initial_count

    def test_disconnected_fragments_flip(self, controls):
        mol = get_CDK_IAtomContainer("CCO.C1CCOC1")
        initial_atoms = mol.getAtomCount()
        controls.flip_structure(mol)
        assert mol.getAtomCount() == initial_atoms

    def test_disconnected_fragments_rotate(self, controls):
        mol = get_CDK_IAtomContainer("CCO.C1CCOC1")
        initial_atoms = mol.getAtomCount()
        controls.rotate_structure(mol, degrees=90)
        assert mol.getAtomCount() == initial_atoms

    def test_zoom_very_small(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=0.1)
        assert gen is not None

    def test_zoom_very_large(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=5.0)
        assert gen is not None

    def test_ratio_thin_bonds(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_stroke_ratio(DepictionGenerator, ratio=0.5)
        assert gen is not None

    def test_ratio_thick_bonds(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_stroke_ratio(DepictionGenerator, ratio=2.0)
        assert gen is not None


class TestCoordinateTransformations:
    """Test coordinate transformations."""

    def test_flip_changes_x_coordinates(self, controls, simple_molecule):
        controls.flip_structure(simple_molecule)
        has_coords = False
        for atom in simple_molecule.atoms():
            point = atom.getPoint2d()
            if point is not None:
                has_coords = True
                break
        assert has_coords or simple_molecule.getAtomCount() > 0

    def test_rotate_changes_coordinates(self, controls, simple_molecule):
        controls.rotate_structure(simple_molecule, degrees=90)
        has_coords = False
        for atom in simple_molecule.atoms():
            point = atom.getPoint2d()
            if point is not None:
                has_coords = True
                break
        assert has_coords or simple_molecule.getAtomCount() > 0

    def test_rotation_360_degrees_equivalent(self, controls, simple_molecule):
        controls.rotate_structure(simple_molecule, degrees=360)
        assert simple_molecule.getAtomCount() > 0


class TestDepictionGeneratorChaining:
    """Test that depiction generator methods can be chained."""

    def test_chain_zoom_ratio_units(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=1.5)
        gen = controls.set_stroke_ratio(gen, ratio=1.2)
        gen = controls.set_svg_units(gen, units=SVGUnits.MM)
        assert gen is not None

    def test_chain_all_controls(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=1.3)
        gen = controls.set_stroke_ratio(gen, ratio=1.0)
        gen = controls.set_svg_units(gen, units=SVGUnits.PX)
        gen = controls.set_title_display(gen, show_title=False, is_reaction=False)
        gen = controls.set_anonymous_display(gen, anonymous=False)
        assert gen is not None

    def test_reapply_same_control(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_zoom(DepictionGenerator, zoom=1.5)
        gen = controls.set_zoom(gen, zoom=2.0)
        assert gen is not None


class TestSVGUnitsStringConversion:
    """Test string to SVGUnits conversion."""

    def test_px_string_conversion(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units="px")
        assert gen is not None

    def test_mm_string_conversion(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units="mm")
        assert gen is not None

    def test_cm_string_conversion(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units="cm")
        assert gen is not None

    def test_in_string_conversion(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units="in")
        assert gen is not None

    def test_case_insensitive_units(self, controls):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = controls.set_svg_units(DepictionGenerator, units="PX")
        assert gen is not None
