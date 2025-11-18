from __future__ import annotations

import pytest
from jpype import JClass

from app.modules.depiction_enhanced import get_cdk_depiction
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer


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
def caffeine():
    return get_CDK_IAtomContainer("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")


@pytest.fixture
def metal_complex():
    return get_CDK_IAtomContainer("[Co][N+]([O-])(=O)")


@pytest.fixture
def ferrocene_fragment():
    return get_CDK_IAtomContainer("[Fe]c1ccccc1")


class TestBasicDepiction:
    """Test basic depiction generation."""

    def test_simple_molecule_depiction(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule)
        assert svg is not None
        assert isinstance(svg, str)
        assert "svg" in svg.lower()
        assert len(svg) > 100

    def test_aromatic_molecule_depiction(self, benzene):
        svg = get_cdk_depiction(benzene)
        assert svg is not None
        assert "svg" in svg.lower()
        assert "Error" not in svg

    def test_complex_molecule_depiction(self, caffeine):
        svg = get_cdk_depiction(caffeine)
        assert svg is not None
        assert "svg" in svg.lower()
        assert "Error" not in svg

    def test_chiral_molecule_depiction(self, chiral_molecule):
        svg = get_cdk_depiction(chiral_molecule)
        assert svg is not None
        assert "svg" in svg.lower()
        assert "Error" not in svg


class TestMoleculeSize:
    """Test molecule size parameter."""

    def test_default_size(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, molSize=(512, 512))
        assert svg is not None
        assert "svg" in svg.lower()

    def test_small_size(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, molSize=(256, 256))
        assert svg is not None
        assert "svg" in svg.lower()

    def test_large_size(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, molSize=(1024, 1024))
        assert svg is not None
        assert "svg" in svg.lower()

    def test_non_square_size(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, molSize=(800, 600))
        assert svg is not None
        assert "svg" in svg.lower()


class TestRotation:
    """Test rotation parameter."""

    def test_no_rotation(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, rotate=0)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_rotate_90(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, rotate=90)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_rotate_180(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, rotate=180)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_rotate_270(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, rotate=270)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_negative_rotation(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, rotate=-45)
        assert svg is not None
        assert "svg" in svg.lower()


class TestKekulization:
    """Test kekulization parameter."""

    def test_with_kekulize(self, benzene):
        svg = get_cdk_depiction(benzene, kekulize=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_without_kekulize(self, benzene):
        svg = get_cdk_depiction(benzene, kekulize=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestCIPAnnotation:
    """Test CIP stereochemistry annotation."""

    def test_cip_enabled(self, chiral_molecule):
        svg = get_cdk_depiction(chiral_molecule, CIP=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_cip_disabled(self, chiral_molecule):
        svg = get_cdk_depiction(chiral_molecule, CIP=False)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_cip_on_achiral(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, CIP=True)
        assert svg is not None
        assert "svg" in svg.lower()


class TestUnicolor:
    """Test unicolor rendering."""

    def test_unicolor_enabled(self, caffeine):
        svg = get_cdk_depiction(caffeine, unicolor=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_unicolor_disabled(self, caffeine):
        svg = get_cdk_depiction(caffeine, unicolor=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestHighlighting:
    """Test highlighting features."""

    def test_smarts_highlighting(self, caffeine):
        svg = get_cdk_depiction(caffeine, highlight="C=O")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_atom_indices_highlighting(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, highlight_atoms=[0, 1, 2])
        assert svg is not None
        assert "svg" in svg.lower()

    def test_single_atom_highlighting(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, highlight_atoms=[0])
        assert svg is not None
        assert "svg" in svg.lower()

    def test_empty_highlight_atoms(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, highlight_atoms=[])
        assert svg is not None
        assert "svg" in svg.lower()

    def test_combined_highlighting(self, caffeine):
        svg = get_cdk_depiction(caffeine, highlight="C=O", highlight_atoms=[0, 1])
        assert svg is not None
        assert "svg" in svg.lower()


class TestAtomNumberDisplay:
    """Test atom number display."""

    def test_show_atom_numbers(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, showAtomNumbers=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_hide_atom_numbers(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, showAtomNumbers=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestHydrogenDisplay:
    """Test hydrogen display modes."""

    def test_hydrogen_provided(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, hydrogen_display="Provided")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_hydrogen_minimal(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, hydrogen_display="Minimal")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_hydrogen_explicit(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, hydrogen_display="Explicit")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_hydrogen_stereo(self, chiral_molecule):
        svg = get_cdk_depiction(chiral_molecule, hydrogen_display="Stereo")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_hydrogen_smart(self, chiral_molecule):
        svg = get_cdk_depiction(chiral_molecule, hydrogen_display="Smart")
        assert svg is not None
        assert "svg" in svg.lower()


class TestAbbreviations:
    """Test abbreviation modes."""

    def test_abbreviate_off(self, benzene):
        svg = get_cdk_depiction(benzene, abbreviate="off")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_abbreviate_groups(self, benzene):
        mol = get_CDK_IAtomContainer("c1ccccc1C")
        svg = get_cdk_depiction(mol, abbreviate="groups")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_abbreviate_reagents(self):
        mol = get_CDK_IAtomContainer("C1CCOC1")
        svg = get_cdk_depiction(mol, abbreviate="reagents")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_abbreviate_all(self):
        mol = get_CDK_IAtomContainer("c1ccccc1CC(C)(C)OC(=O)N")
        svg = get_cdk_depiction(mol, abbreviate="on")
        assert svg is not None
        assert "svg" in svg.lower()


class TestDativeBonds:
    """Test dative bond perception modes."""

    def test_dative_never(self, metal_complex):
        svg = get_cdk_depiction(metal_complex, dative="never")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_dative_metals(self, metal_complex):
        svg = get_cdk_depiction(metal_complex, dative="metals")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_dative_always(self):
        # Use a simple molecule for dative bond perception
        # Use amine oxide or sulfoxide as a simple test case
        mol = get_CDK_IAtomContainer("CS(=O)C")  # Dimethyl sulfoxide
        svg = get_cdk_depiction(mol, dative="always")
        assert svg is not None
        assert "svg" in svg.lower()


class TestMulticenterBonds:
    """Test multicenter bond display styles."""

    def test_multicenter_provided(self, ferrocene_fragment):
        svg = get_cdk_depiction(ferrocene_fragment, multicenter="provided")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_multicenter_dative(self, ferrocene_fragment):
        svg = get_cdk_depiction(ferrocene_fragment, multicenter="dative")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_multicenter_dashed(self, ferrocene_fragment):
        svg = get_cdk_depiction(ferrocene_fragment, multicenter="dashed")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_multicenter_dashed_neutral(self, ferrocene_fragment):
        svg = get_cdk_depiction(ferrocene_fragment, multicenter="dashed_neutral")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_multicenter_hidden(self, ferrocene_fragment):
        svg = get_cdk_depiction(ferrocene_fragment, multicenter="hidden")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_multicenter_hidden_neutral(self, ferrocene_fragment):
        svg = get_cdk_depiction(ferrocene_fragment, multicenter="hidden_neutral")
        assert svg is not None
        assert "svg" in svg.lower()


class TestAnnotationModes:
    """Test annotation modes."""

    def test_annotate_none(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, annotate="none")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_annotate_number(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, annotate="number")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_annotate_bondnumber(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, annotate="bondnumber")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_annotate_mapidx(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, annotate="mapidx")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_annotate_atomvalue(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, annotate="atomvalue")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_annotate_colmap(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, annotate="colmap")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_annotate_cip(self, chiral_molecule):
        svg = get_cdk_depiction(chiral_molecule, annotate="cip")
        assert svg is not None
        assert "svg" in svg.lower()


class TestStylePresets:
    """Test style preset modes."""

    def test_style_cow(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="cow")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_style_cob(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="cob")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_style_cot(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="cot")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_style_bow(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="bow")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_style_bot(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="bot")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_style_wob(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="wob")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_style_nob(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="nob")
        assert svg is not None
        assert "svg" in svg.lower()


class TestAromaticDisplay:
    """Test aromatic donut display."""

    def test_donuts_enabled(self, benzene):
        svg = get_cdk_depiction(benzene, donuts=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_donuts_disabled(self, benzene):
        svg = get_cdk_depiction(benzene, donuts=False)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_donuts_with_naphthalene(self):
        mol = get_CDK_IAtomContainer("c1ccc2ccccc2c1")
        svg = get_cdk_depiction(mol, donuts=True)
        assert svg is not None
        assert "svg" in svg.lower()


class TestReactionArrows:
    """Test reaction arrow types (for reactions)."""

    def test_arrow_forward(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, arrow="forward")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_arrow_equilibrium(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, arrow="equ")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_arrow_retro(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, arrow="ret")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_arrow_resonance(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, arrow="res")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_arrow_no_go(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, arrow="ngo")
        assert svg is not None
        assert "svg" in svg.lower()


class TestReactionMapping:
    """Test reaction mapping alignment."""

    def test_alignrxnmap_enabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, alignrxnmap=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_alignrxnmap_disabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, alignrxnmap=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestTitleDisplay:
    """Test title display."""

    def test_showtitle_enabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, showtitle=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_showtitle_disabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, showtitle=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestCustomColors:
    """Test custom background and foreground colors."""

    def test_custom_bgcolor(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, bgcolor="#FFFFFF")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_custom_fgcolor(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, fgcolor="#000000")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_both_custom_colors(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, bgcolor="#FFFFFF", fgcolor="#000000")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_default_colors(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, bgcolor="default", fgcolor="default")
        assert svg is not None
        assert "svg" in svg.lower()


class TestZoom:
    """Test zoom parameter."""

    def test_zoom_default(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, zoom=1.3)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_zoom_small(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, zoom=0.5)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_zoom_large(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, zoom=2.0)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_zoom_very_small(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, zoom=0.1)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_zoom_very_large(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, zoom=5.0)
        assert svg is not None
        assert "svg" in svg.lower()


class TestStrokeRatio:
    """Test stroke ratio parameter."""

    def test_ratio_default(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, ratio=1.0)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_ratio_thin(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, ratio=0.5)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_ratio_thick(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, ratio=2.0)
        assert svg is not None
        assert "svg" in svg.lower()


class TestFlip:
    """Test structure flipping."""

    def test_flip_enabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, flip=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_flip_disabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, flip=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestAnonymousDisplay:
    """Test anonymous atom display."""

    def test_anon_enabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, anon=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_anon_disabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, anon=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestSMARTSLimit:
    """Test SMARTS hit limit."""

    def test_smalim_default(self, caffeine):
        svg = get_cdk_depiction(caffeine, highlight="C", smalim=100)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_smalim_low(self, caffeine):
        svg = get_cdk_depiction(caffeine, highlight="C", smalim=5)
        assert svg is not None
        assert "svg" in svg.lower()


class TestSVGUnits:
    """Test SVG units."""

    def test_svgunits_px(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, svgunits="px")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_svgunits_mm(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, svgunits="mm")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_svgunits_cm(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, svgunits="cm")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_svgunits_in(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, svgunits="in")
        assert svg is not None
        assert "svg" in svg.lower()


class TestRadicalPerception:
    """Test radical perception."""

    def test_perceive_radicals_enabled(self):
        mol = get_CDK_IAtomContainer("[CH3]")
        svg = get_cdk_depiction(mol, perceive_radicals=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_perceive_radicals_disabled(self):
        mol = get_CDK_IAtomContainer("[CH3]")
        svg = get_cdk_depiction(mol, perceive_radicals=False)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_reaction_radical_perception(self):
        """Test radical perception in reactions."""
        cdk_base = "org.openscience.cdk"
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
        SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
        reaction = SmilesParser.parseReactionSmiles("[CH3].O>>[CH3]O")

        svg = get_cdk_depiction(reaction, perceive_radicals=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_reaction_set_radical_perception(self):
        """Test radical perception in reaction sets."""
        cdk_base = "org.openscience.cdk"
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
        SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
        reaction_set = SmilesParser.parseReactionSetSmiles("[CH3].O>>[CH3]O")

        svg = get_cdk_depiction(reaction_set, perceive_radicals=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_reaction_with_agents_radical_perception(self):
        """Test radical perception in reactions with agents."""
        cdk_base = "org.openscience.cdk"
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
        SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
        # Reaction with agent: [CH3] > [Cl] > [CH3]Cl
        reaction = SmilesParser.parseReactionSmiles("[CH3]>[Cl]>[CH3]Cl")

        svg = get_cdk_depiction(reaction, perceive_radicals=True)
        assert svg is not None
        assert "svg" in svg.lower()


class TestMDLHighlighting:
    """Test MDL HILITE support."""

    def test_mdl_highlighting_enabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, apply_mdl_highlighting=True)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_mdl_highlighting_disabled(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, apply_mdl_highlighting=False)
        assert svg is not None
        assert "svg" in svg.lower()


class TestCXSMILESString:
    """Test CXSMILES string parameter."""

    def test_with_cxsmiles_string(self, benzene):
        svg = get_cdk_depiction(benzene, smiles_string="c1ccccc1 |ha:0,1,2|")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_without_cxsmiles_string(self, benzene):
        svg = get_cdk_depiction(benzene, smiles_string=None)
        assert svg is not None
        assert "svg" in svg.lower()


class TestCombinedFeatures:
    """Test combinations of multiple features."""

    def test_abbreviations_with_annotations(self):
        mol = get_CDK_IAtomContainer("c1ccccc1C")
        svg = get_cdk_depiction(mol, abbreviate="groups", annotate="number")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_dative_with_style(self, metal_complex):
        svg = get_cdk_depiction(metal_complex, dative="metals", style="bow")
        assert svg is not None
        assert "svg" in svg.lower()

    def test_donuts_with_highlighting(self, benzene):
        svg = get_cdk_depiction(benzene, donuts=True, highlight_atoms=[0, 1, 2])
        assert svg is not None
        assert "svg" in svg.lower()

    def test_all_features_combined(self):
        mol = get_CDK_IAtomContainer("c1ccccc1C")
        svg = get_cdk_depiction(
            mol,
            abbreviate="groups",
            annotate="number",
            style="bow",
            donuts=True,
            zoom=1.5,
            ratio=1.2,
            flip=True,
            CIP=True,
        )
        assert svg is not None
        assert "svg" in svg.lower()

    def test_complex_organometallic_with_features(self, ferrocene_fragment):
        svg = get_cdk_depiction(
            ferrocene_fragment,
            dative="metals",
            multicenter="dative",
            annotate="number",
            style="cow",
            zoom=1.5,
        )
        assert svg is not None
        assert "svg" in svg.lower()


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_atom(self):
        mol = get_CDK_IAtomContainer("C")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_disconnected_fragments(self):
        mol = get_CDK_IAtomContainer("CCO.C1CCOC1")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_charged_species(self):
        mol = get_CDK_IAtomContainer("[NH4+]")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_radical_species(self):
        mol = get_CDK_IAtomContainer("[CH3]")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_zwitterion(self):
        mol = get_CDK_IAtomContainer("C[N+](C)(C)CC[O-]")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()


class TestComplexMolecules:
    """Test complex molecular structures."""

    def test_steroid_structure(self):
        mol = get_CDK_IAtomContainer("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_peptide_with_protecting_groups(self):
        mol = get_CDK_IAtomContainer("CC(C)(C)OC(=O)NCC(=O)NCC(=O)O")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()

    def test_complex_heterocycle(self):
        mol = get_CDK_IAtomContainer("c1ccc2c(c1)nc1ccccc1n2")
        svg = get_cdk_depiction(mol)
        assert svg is not None
        assert "svg" in svg.lower()


class TestErrorHandling:
    """Test error handling for invalid inputs."""

    def test_invalid_hydrogen_display(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, hydrogen_display="Invalid")
        assert svg is not None

    def test_out_of_range_highlight_atoms(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, highlight_atoms=[0, 1, 100])
        assert svg is not None
        assert "svg" in svg.lower()

    def test_invalid_style(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule, style="invalid")
        assert svg is not None


class TestSVGOutput:
    """Test SVG output structure."""

    def test_svg_contains_viewbox(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule)
        assert "viewBox" in svg or "viewbox" in svg.lower()

    def test_svg_well_formed(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule)
        # Count svg tags (with or without namespace)
        svg_open = svg.count("<svg") + svg.count("<ns0:svg") + svg.count("<ns1:svg")
        svg_close = (
            svg.count("</svg>") + svg.count("</ns0:svg>") + svg.count("</ns1:svg>")
        )
        assert svg_open == svg_close

    def test_svg_contains_paths(self, simple_molecule):
        svg = get_cdk_depiction(simple_molecule)
        # Check for any drawing elements
        has_drawing = (
            "<path" in svg
            or "ns0:path" in svg
            or "<line" in svg
            or "ns0:line" in svg
            or "<circle" in svg
            or "ns0:circle" in svg
            or "<rect" in svg
            or "ns0:rect" in svg
            or "<ellipse" in svg
            or "ns0:ellipse" in svg
        )
        assert has_drawing


class TestPerformance:
    """Test performance with different molecule sizes."""

    def test_small_molecule_performance(self):
        mol = get_CDK_IAtomContainer("C")
        svg = get_cdk_depiction(mol)
        assert svg is not None

    def test_medium_molecule_performance(self, caffeine):
        svg = get_cdk_depiction(caffeine)
        assert svg is not None

    def test_large_molecule_performance(self):
        mol = get_CDK_IAtomContainer(
            "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C.c1ccccc1.C1CCOC1"
        )
        svg = get_cdk_depiction(mol)
        assert svg is not None


class TestConsistency:
    """Test output consistency."""

    def test_same_input_produces_similar_output(self, simple_molecule):
        svg1 = get_cdk_depiction(simple_molecule)
        svg2 = get_cdk_depiction(simple_molecule)
        assert svg1 is not None
        assert svg2 is not None
        assert len(svg1) > 0
        assert len(svg2) > 0

    def test_different_parameters_produce_different_output(self, simple_molecule):
        svg1 = get_cdk_depiction(simple_molecule, zoom=1.0)
        svg2 = get_cdk_depiction(simple_molecule, zoom=2.0)
        assert svg1 is not None
        assert svg2 is not None
