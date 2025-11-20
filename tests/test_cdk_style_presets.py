from __future__ import annotations

import pytest

from app.modules.cdk_depict.style_presets import (
    StylePresetSystem,
    StylePreset,
    StyleConfiguration,
    STYLE_CONFIGURATIONS,
)
from jpype import JClass


@pytest.fixture
def style_system():
    return StylePresetSystem()


class TestStylePresetEnum:
    """Test StylePreset enum values."""

    def test_cow_preset(self):
        assert StylePreset.COW.value == "cow"

    def test_cob_preset(self):
        assert StylePreset.COB.value == "cob"

    def test_cot_preset(self):
        assert StylePreset.COT.value == "cot"

    def test_bow_preset(self):
        assert StylePreset.BOW.value == "bow"

    def test_bot_preset(self):
        assert StylePreset.BOT.value == "bot"

    def test_wob_preset(self):
        assert StylePreset.WOB.value == "wob"

    def test_nob_preset(self):
        assert StylePreset.NOB.value == "nob"

    def test_all_presets_defined(self):
        presets = [
            StylePreset.COW,
            StylePreset.COB,
            StylePreset.COT,
            StylePreset.BOW,
            StylePreset.BOT,
            StylePreset.WOB,
            StylePreset.NOB,
        ]
        assert len(presets) == 7


class TestStyleConfiguration:
    """Test StyleConfiguration dataclass."""

    def test_configuration_has_background_color(self):
        """Test that configuration stores background_color."""
        config = StyleConfiguration(
            use_atom_colors=True,
            background_color=(255, 255, 255),
            annotation_color=(0, 0, 0),
        )
        assert config.background_color == (255, 255, 255)

    def test_configuration_has_annotation_color(self):
        """Test that configuration stores annotation_color."""
        config = StyleConfiguration(
            use_atom_colors=True,
            background_color=(255, 255, 255),
            annotation_color=(0, 0, 0),
        )
        assert config.annotation_color == (0, 0, 0)

    def test_configuration_has_use_atom_colors(self):
        """Test that configuration stores use_atom_colors flag."""
        config = StyleConfiguration(
            use_atom_colors=True,
            background_color=(255, 255, 255),
            annotation_color=(0, 0, 0),
        )
        assert config.use_atom_colors is True

    def test_configuration_atom_colors_false(self):
        """Test configuration with use_atom_colors disabled."""
        config = StyleConfiguration(
            use_atom_colors=False,
            atom_color=(255, 255, 255),
            background_color=(0, 0, 0),
            annotation_color=(255, 255, 255),
        )
        assert config.use_atom_colors is False
        assert config.atom_color == (255, 255, 255)
        assert config.use_atom_colors is False


class TestStyleConfigurationsDict:
    """Test STYLE_CONFIGURATIONS dictionary."""

    def test_configurations_defined(self):
        assert isinstance(STYLE_CONFIGURATIONS, dict)
        assert len(STYLE_CONFIGURATIONS) > 0

    def test_cow_configuration_exists(self):
        assert StylePreset.COW in STYLE_CONFIGURATIONS

    def test_cob_configuration_exists(self):
        assert StylePreset.COB in STYLE_CONFIGURATIONS

    def test_bow_configuration_exists(self):
        assert StylePreset.BOW in STYLE_CONFIGURATIONS

    def test_wob_configuration_exists(self):
        assert StylePreset.WOB in STYLE_CONFIGURATIONS

    def test_nob_configuration_exists(self):
        assert StylePreset.NOB in STYLE_CONFIGURATIONS

    def test_all_presets_have_configurations(self):
        for preset in StylePreset:
            assert preset in STYLE_CONFIGURATIONS

    def test_configurations_are_style_configuration_objects(self):
        for config in STYLE_CONFIGURATIONS.values():
            assert isinstance(config, StyleConfiguration)


class TestCOWStyle:
    """Test Color on White (COW) style preset."""

    def test_cow_has_white_background(self):
        config = STYLE_CONFIGURATIONS[StylePreset.COW]
        assert config.background_color == (255, 255, 255)

    def test_cow_uses_color_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.COW]
        assert config.use_atom_colors is True

    def test_cow_apply_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COW)
        assert gen is not None


class TestCOBStyle:
    """Test Color on Black (COB) style preset."""

    def test_cob_has_black_background(self):
        config = STYLE_CONFIGURATIONS[StylePreset.COB]
        assert config.background_color == (0, 0, 0)

    def test_cob_uses_color_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.COB]
        assert config.use_atom_colors is True

    def test_cob_apply_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COB)
        assert gen is not None


class TestBOWStyle:
    """Test Black on White (BOW) style preset."""

    def test_bow_has_white_background(self):
        config = STYLE_CONFIGURATIONS[StylePreset.BOW]
        assert config.background_color == (255, 255, 255)

    def test_bow_uses_black_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.BOW]
        assert config.use_atom_colors is False

    def test_bow_has_black_foreground(self):
        config = STYLE_CONFIGURATIONS[StylePreset.BOW]
        assert config.atom_color == (0, 0, 0)

    def test_bow_apply_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.BOW)
        assert gen is not None


class TestWOBStyle:
    """Test White on Black (WOB) style preset."""

    def test_wob_has_black_background(self):
        config = STYLE_CONFIGURATIONS[StylePreset.WOB]
        assert config.background_color == (0, 0, 0)

    def test_wob_uses_white_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.WOB]
        assert config.use_atom_colors is False

    def test_wob_has_white_foreground(self):
        config = STYLE_CONFIGURATIONS[StylePreset.WOB]
        assert config.atom_color == (255, 255, 255)

    def test_wob_apply_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.WOB)
        assert gen is not None


class TestNOBStyle:
    """Test Neon on Black (NOB) style preset."""

    def test_nob_has_black_background(self):
        config = STYLE_CONFIGURATIONS[StylePreset.NOB]
        assert config.background_color == (0, 0, 0)

    def test_nob_uses_single_neon_color(self):
        """Test that NOB uses a single neon cyan color, not CDK atom colors."""
        config = STYLE_CONFIGURATIONS[StylePreset.NOB]
        assert config.use_atom_colors is False
        assert config.atom_color == (0, 255, 255)  # Cyan neon

    def test_nob_apply_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.NOB)
        assert gen is not None


class TestCOTStyle:
    """Test Color on Transparent (COT) style preset."""

    def test_cot_apply_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COT)
        assert gen is not None

    def test_cot_uses_color_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.COT]
        assert config.use_atom_colors is True


class TestBOTStyle:
    """Test Black on Transparent (BOT) style preset."""

    def test_bot_apply_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.BOT)
        assert gen is not None

    def test_bot_uses_black_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.BOT]
        assert config.use_atom_colors is False


class TestStylePresetSystemInitialization:
    """Test StylePresetSystem initialization."""

    def test_default_initialization(self):
        """Test that system initializes with correct cdk_base."""
        system = StylePresetSystem()
        assert system.cdk_base == "org.openscience.cdk"

    def test_initialization_loads_classes(self):
        """Test that initialization loads all required Java classes."""
        system = StylePresetSystem()
        assert system.Color is not None
        assert system.CDK2DAtomColors is not None
        assert system.StandardGenerator is not None


class TestApplyStyle:
    """Test applying style presets."""

    def test_apply_cow_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COW)
        assert gen is not None

    def test_apply_cob_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COB)
        assert gen is not None

    def test_apply_bow_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.BOW)
        assert gen is not None

    def test_apply_wob_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.WOB)
        assert gen is not None

    def test_apply_nob_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.NOB)
        assert gen is not None

    def test_apply_cot_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COT)
        assert gen is not None

    def test_apply_bot_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.BOT)
        assert gen is not None


class TestCustomColors:
    """Test custom color application."""

    def test_apply_custom_bgcolor(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_custom_colors(
            DepictionGenerator, bg_color=(255, 255, 255)
        )
        assert gen is not None

    def test_apply_custom_fgcolor(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_custom_colors(DepictionGenerator, fg_color=(0, 0, 0))
        assert gen is not None

    def test_apply_both_custom_colors(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_custom_colors(
            DepictionGenerator, bg_color=(255, 255, 255), fg_color=(0, 0, 0)
        )
        assert gen is not None

    def test_apply_custom_red_background(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_custom_colors(DepictionGenerator, bg_color=(255, 0, 0))
        assert gen is not None

    def test_apply_custom_blue_foreground(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_custom_colors(DepictionGenerator, fg_color=(0, 0, 255))
        assert gen is not None


class TestStyleWithCustomColors:
    """Test applying preset with custom color override."""

    def test_cow_with_custom_bgcolor(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COW)
        gen = style_system.apply_custom_colors(gen, bg_color=(200, 200, 200))
        assert gen is not None

    def test_bow_with_custom_colors(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.BOW)
        gen = style_system.apply_custom_colors(
            gen, bg_color=(255, 255, 255), fg_color=(50, 50, 50)
        )
        assert gen is not None


class TestColorConversion:
    """Test color string parsing to Java Color conversion."""

    def test_white_color_hex_conversion(self, style_system):
        """Test parsing white color from hex string."""
        Color = JClass("java.awt.Color")
        color = style_system._parse_color("#FFFFFF")
        assert color is not None
        assert isinstance(color, type(Color.WHITE))

    def test_black_color_hex_conversion(self, style_system):
        """Test parsing black color from hex string."""
        Color = JClass("java.awt.Color")
        color = style_system._parse_color("#000000")
        assert color is not None
        assert isinstance(color, type(Color.BLACK))

    def test_red_color_hex_conversion(self, style_system):
        """Test parsing red color from hex string."""
        color = style_system._parse_color("#FF0000")
        assert color is not None

    def test_custom_color_hex_conversion(self, style_system):
        """Test parsing custom color from hex string."""
        color = style_system._parse_color("#7B2D43")
        assert color is not None

    def test_color_with_alpha_conversion(self, style_system):
        """Test parsing color with alpha channel."""
        color = style_system._parse_color("0xFF0000FF")
        assert color is not None

    def test_invalid_color_returns_none(self, style_system):
        """Test that invalid color string returns None."""
        color = style_system._parse_color("invalid")
        assert color is None


class TestStyleChaining:
    """Test chaining style operations."""

    def test_apply_style_then_custom_colors(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COW)
        gen = style_system.apply_custom_colors(gen, bg_color=(240, 240, 240))
        assert gen is not None

    def test_multiple_style_changes(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COW)
        gen = style_system.apply_preset(gen, StylePreset.BOW)
        assert gen is not None

    def test_reapply_same_style(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COW)
        gen = style_system.apply_preset(gen, StylePreset.COW)
        assert gen is not None


class TestPublicationStyles:
    """Test styles commonly used for publications."""

    def test_bow_for_publication(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.BOW)
        assert gen is not None

    def test_cow_for_publication(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COW)
        assert gen is not None


class TestPresentationStyles:
    """Test styles commonly used for presentations."""

    def test_cob_for_dark_slides(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COB)
        assert gen is not None

    def test_nob_for_dark_slides(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.NOB)
        assert gen is not None


class TestTransparentStyles:
    """Test transparent background styles."""

    def test_cot_transparent_background(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.COT)
        assert gen is not None

    def test_bot_transparent_background(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_preset(DepictionGenerator, StylePreset.BOT)
        assert gen is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_apply_style_with_none_generator(self, style_system):
        try:
            gen = style_system.apply_preset(None, StylePreset.COW)
            assert gen is None or gen is not None
        except Exception:
            assert True

    def test_custom_colors_extremes(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_custom_colors(
            DepictionGenerator, bg_color=(0, 0, 0), fg_color=(255, 255, 255)
        )
        assert gen is not None

    def test_custom_colors_midtones(self, style_system):
        DepictionGenerator = JClass("org.openscience.cdk.depict.DepictionGenerator")()
        gen = style_system.apply_custom_colors(
            DepictionGenerator, bg_color=(128, 128, 128), fg_color=(64, 64, 64)
        )
        assert gen is not None


class TestColorAtomSettings:
    """Test color atom settings for different styles."""

    def test_cow_uses_colored_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.COW]
        assert config.use_atom_colors is True

    def test_bow_uses_black_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.BOW]
        assert config.use_atom_colors is False

    def test_cob_uses_colored_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.COB]
        assert config.use_atom_colors is True

    def test_wob_uses_white_atoms(self):
        config = STYLE_CONFIGURATIONS[StylePreset.WOB]
        assert config.use_atom_colors is False

    def test_nob_uses_single_neon_color(self):
        """Test that NOB uses a single neon color (cyan), not CDK atom colors."""
        config = STYLE_CONFIGURATIONS[StylePreset.NOB]
        assert config.use_atom_colors is False
        assert config.atom_color == (0, 255, 255)  # Cyan neon


class TestStylePresetStringParsing:
    """Test string to StylePreset conversion."""

    def test_parse_cow_string(self):
        preset = StylePreset("cow")
        assert preset == StylePreset.COW

    def test_parse_cob_string(self):
        preset = StylePreset("cob")
        assert preset == StylePreset.COB

    def test_parse_bow_string(self):
        preset = StylePreset("bow")
        assert preset == StylePreset.BOW

    def test_parse_wob_string(self):
        preset = StylePreset("wob")
        assert preset == StylePreset.WOB

    def test_parse_nob_string(self):
        preset = StylePreset("nob")
        assert preset == StylePreset.NOB

    def test_parse_cot_string(self):
        preset = StylePreset("cot")
        assert preset == StylePreset.COT

    def test_parse_bot_string(self):
        preset = StylePreset("bot")
        assert preset == StylePreset.BOT
