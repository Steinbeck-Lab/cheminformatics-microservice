"""Style Presets for CDK Molecular Depictions.

This module provides predefined color schemes for molecular structure depictions,
including variations for different backgrounds (white, black, transparent) and
color modes (color, black/white, neon).

Based on CDK DepictController.java style preset functionality.
"""

from __future__ import annotations

import logging
from enum import Enum
from typing import Tuple, Optional, Dict, Any

from jpype import JClass

logger = logging.getLogger(__name__)


class StylePreset(Enum):
    """Enumeration of supported style presets.

    Attributes:
        COW: Color on White (default) - Colored atoms on white background
        COB: Color on Black - Colored atoms on black background
        COT: Color on Transparent - Colored atoms on transparent background
        BOW: Black on White - Black atoms on white background
        BOT: Black on Transparent - Black atoms on transparent background
        WOB: White on Black - White atoms on black background
        NOB: Neon on Black - Neon-colored atoms on black background
    """

    COW = "cow"  # Color on White (default)
    COB = "cob"  # Color on Black
    COT = "cot"  # Color on Transparent
    BOW = "bow"  # Black on White
    BOT = "bot"  # Black on Transparent
    WOB = "wob"  # White on Black
    NOB = "nob"  # Neon on Black


class StyleConfiguration:
    """Configuration for a style preset.

    Attributes:
        use_atom_colors: Whether to use atom-specific colors
        atom_color: Single color for all atoms (if not using atom colors)
        background_color: Background color (None for transparent)
        annotation_color: Color for annotations (atom numbers, etc.)
        selection_color: Color for highlighting selections
    """

    def __init__(
        self,
        use_atom_colors: bool = True,
        atom_color: Optional[Tuple[int, int, int]] = None,
        background_color: Optional[Tuple[int, int, int]] = (255, 255, 255),
        annotation_color: Tuple[int, int, int] = (0, 0, 0),
        selection_color: Tuple[int, int, int] = (170, 255, 170),
    ):
        """Initialize style configuration.

        Args:
            use_atom_colors: Use CDK2DAtomColors (True) or single color (False)
            atom_color: RGB tuple for single atom color (used if use_atom_colors=False)
            background_color: RGB tuple for background (None for transparent)
            annotation_color: RGB tuple for annotations
            selection_color: RGB tuple for highlights/selections
        """
        self.use_atom_colors = use_atom_colors
        self.atom_color = atom_color
        self.background_color = background_color
        self.annotation_color = annotation_color
        self.selection_color = selection_color


# Style preset configurations (matching CDK Java implementation)
STYLE_CONFIGURATIONS: Dict[StylePreset, StyleConfiguration] = {
    StylePreset.COW: StyleConfiguration(
        use_atom_colors=True,
        background_color=(255, 255, 255),  # White
        annotation_color=(0, 0, 0),  # Black
        selection_color=(170, 255, 170),  # Light green
    ),
    StylePreset.COB: StyleConfiguration(
        use_atom_colors=True,
        background_color=(0, 0, 0),  # Black
        annotation_color=(255, 255, 255),  # White
        selection_color=(255, 170, 170),  # Light red
    ),
    StylePreset.COT: StyleConfiguration(
        use_atom_colors=True,
        background_color=None,  # Transparent
        annotation_color=(0, 0, 0),  # Black
        selection_color=(170, 255, 170),  # Light green
    ),
    StylePreset.BOW: StyleConfiguration(
        use_atom_colors=False,
        atom_color=(0, 0, 0),  # Black
        background_color=(255, 255, 255),  # White
        annotation_color=(0, 0, 0),  # Black
        selection_color=(255, 0, 0),  # Red
    ),
    StylePreset.BOT: StyleConfiguration(
        use_atom_colors=False,
        atom_color=(0, 0, 0),  # Black
        background_color=None,  # Transparent
        annotation_color=(0, 0, 0),  # Black
        selection_color=(255, 0, 0),  # Red
    ),
    StylePreset.WOB: StyleConfiguration(
        use_atom_colors=False,
        atom_color=(255, 255, 255),  # White
        background_color=(0, 0, 0),  # Black
        annotation_color=(255, 255, 255),  # White
        selection_color=(255, 0, 0),  # Red
    ),
    StylePreset.NOB: StyleConfiguration(
        use_atom_colors=False,
        atom_color=(0, 255, 255),  # Cyan (neon)
        background_color=(0, 0, 0),  # Black
        annotation_color=(255, 0, 255),  # Magenta
        selection_color=(255, 255, 0),  # Yellow
    ),
}


class StylePresetSystem:
    """Manages style presets for molecular depictions.

    This class provides methods to apply predefined style presets to CDK
    depiction generators, controlling colors, backgrounds, and visual styling.

    Attributes:
        Color: Java Color class
        UniColor: CDK UniColor class
        CDK2DAtomColors: CDK atom coloring class
        StandardGenerator: CDK StandardGenerator class
    """

    def __init__(self):
        """Initialize the style preset system with CDK classes."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.Color = JClass("java.awt.Color")
            self.UniColor = JClass(self.cdk_base + ".renderer.color.UniColor")
            self.CDK2DAtomColors = JClass(
                self.cdk_base + ".renderer.color.CDK2DAtomColors"
            )()
            self.StandardGenerator = JClass(
                self.cdk_base + ".renderer.generators.standard.StandardGenerator"
            )
            logger.debug("StylePresetSystem initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize StylePresetSystem: {e}")
            raise

    def apply_preset(self, depiction_generator: Any, preset: StylePreset) -> Any:
        """Apply a style preset to the depiction generator.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            preset: StylePreset to apply

        Returns:
            Updated DepictionGenerator with style applied

        Example:
            >>> system = StylePresetSystem()
            >>> gen = get_depiction_generator()
            >>> gen = system.apply_preset(gen, StylePreset.COB)
        """
        try:
            config = STYLE_CONFIGURATIONS[preset]
            logger.debug(f"Applying style preset: {preset.value}")

            # Apply atom colors
            if config.use_atom_colors:
                depiction_generator = depiction_generator.withAtomColors(
                    self.CDK2DAtomColors
                )
            else:
                if config.atom_color:
                    r, g, b = config.atom_color
                    color = self.Color(r, g, b)
                    depiction_generator = depiction_generator.withParam(
                        self.StandardGenerator.AtomColor.class_, self.UniColor(color)
                    )

            # Apply background color
            if config.background_color:
                r, g, b = config.background_color
                bg_color = self.Color(r, g, b)
                depiction_generator = depiction_generator.withBackgroundColor(bg_color)
            else:
                # Transparent background
                transparent = self.Color(0, 0, 0, 0)
                depiction_generator = depiction_generator.withBackgroundColor(
                    transparent
                )

            # Apply annotation color
            r, g, b = config.annotation_color
            annot_color = self.Color(r, g, b)
            depiction_generator = depiction_generator.withAnnotationColor(annot_color)

            # Apply selection/highlight color
            r, g, b = config.selection_color
            select_color = self.Color(r, g, b)
            # Selection color is set via RendererModel
            RendererModel = JClass(self.cdk_base + ".renderer.RendererModel")
            depiction_generator = depiction_generator.withParam(
                RendererModel.SelectionColor.class_, select_color
            )

            logger.info(f"Successfully applied style preset: {preset.value}")
            return depiction_generator

        except Exception as e:
            logger.error(f"Failed to apply style preset ({preset.value}): {e}")
            return depiction_generator

    def apply_custom_colors(
        self,
        depiction_generator: Any,
        bg_color: Optional[str] = None,
        fg_color: Optional[str] = None,
    ) -> Any:
        """Apply custom background and foreground colors.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            bg_color: Background color (hex "#RRGGBB" or "0xRRGGBBAA", or "default")
            fg_color: Foreground/atom color (hex "#RRGGBB" or "0xRRGGBBAA", or "default")

        Returns:
            Updated DepictionGenerator with custom colors

        Example:
            >>> system = StylePresetSystem()
            >>> gen = system.apply_custom_colors(gen, bg_color="#000000", fg_color="#FFFFFF")
        """
        try:
            # Apply background color
            if bg_color and bg_color != "default":
                color = self._parse_color(bg_color)
                if color:
                    depiction_generator = depiction_generator.withBackgroundColor(color)
                    logger.debug(f"Applied background color: {bg_color}")

            # Apply foreground color
            if fg_color and fg_color != "default":
                color = self._parse_color(fg_color)
                if color:
                    depiction_generator = depiction_generator.withParam(
                        self.StandardGenerator.AtomColor.class_, self.UniColor(color)
                    )
                    logger.debug(f"Applied foreground color: {fg_color}")

            return depiction_generator

        except Exception as e:
            logger.error(f"Failed to apply custom colors: {e}")
            return depiction_generator

    def _parse_color(self, color_str: str) -> Optional[Any]:
        """Parse color string to Java Color object.

        Supports:
        - Hex format: "#RRGGBB" or "#RRGGBBAA"
        - Integer format: "0xRRGGBB" or "0xRRGGBBAA"

        Args:
            color_str: Color string to parse

        Returns:
            Java Color object or None if parsing fails
        """
        try:
            # Remove prefix and convert to integer
            if color_str.startswith("#"):
                hex_str = color_str[1:]
            elif color_str.startswith("0x") or color_str.startswith("0X"):
                hex_str = color_str[2:]
            else:
                hex_str = color_str

            # Parse as integer
            color_int = int(hex_str, 16)

            # Extract RGBA components
            if len(hex_str) == 6:  # RGB
                r = (color_int >> 16) & 0xFF
                g = (color_int >> 8) & 0xFF
                b = color_int & 0xFF
                return self.Color(r, g, b)
            elif len(hex_str) == 8:  # RGBA
                r = (color_int >> 24) & 0xFF
                g = (color_int >> 16) & 0xFF
                b = (color_int >> 8) & 0xFF
                a = color_int & 0xFF
                return self.Color(r, g, b, a)
            else:
                logger.warning(f"Invalid color format: {color_str}")
                return None

        except Exception as e:
            logger.warning(f"Failed to parse color '{color_str}': {e}")
            return None


def get_style_preset(preset_str: str) -> StylePreset:
    """Convert string to StylePreset enum.

    Args:
        preset_str: String representation of style preset

    Returns:
        Corresponding StylePreset enum value

    Raises:
        ValueError: If preset_str is not a valid style preset

    Example:
        >>> preset = get_style_preset("cow")
        >>> assert preset == StylePreset.COW
    """
    try:
        return StylePreset(preset_str.lower())
    except ValueError:
        valid_presets = [p.value for p in StylePreset]
        raise ValueError(
            f"Invalid style preset: '{preset_str}'. "
            f"Valid presets are: {', '.join(valid_presets)}"
        )


def apply_style_preset(depiction_generator: Any, preset: str = "cow") -> Any:
    """Convenience function to apply a style preset.

    Args:
        depiction_generator: CDK DepictionGenerator instance
        preset: Style preset as string (default: "cow")

    Returns:
        Updated DepictionGenerator with style applied

    Example:
        >>> gen = apply_style_preset(gen, preset="cob")
        >>> gen = apply_style_preset(gen, preset="bow")
    """
    try:
        style_preset = get_style_preset(preset)
        system = StylePresetSystem()
        return system.apply_preset(depiction_generator, style_preset)
    except Exception as e:
        logger.error(f"Failed to apply style preset: {e}")
        return depiction_generator


def apply_custom_colors(
    depiction_generator: Any,
    bg_color: Optional[str] = None,
    fg_color: Optional[str] = None,
) -> Any:
    """Convenience function to apply custom colors.

    Args:
        depiction_generator: CDK DepictionGenerator instance
        bg_color: Background color (hex or default)
        fg_color: Foreground color (hex or default)

    Returns:
        Updated DepictionGenerator with custom colors

    Example:
        >>> gen = apply_custom_colors(gen, bg_color="#000000", fg_color="#FFFFFF")
    """
    try:
        system = StylePresetSystem()
        return system.apply_custom_colors(depiction_generator, bg_color, fg_color)
    except Exception as e:
        logger.error(f"Failed to apply custom colors: {e}")
        return depiction_generator


# Export public API
__all__ = [
    "StylePreset",
    "StyleConfiguration",
    "STYLE_CONFIGURATIONS",
    "StylePresetSystem",
    "get_style_preset",
    "apply_style_preset",
    "apply_custom_colors",
]
