"""Advanced Controls for CDK Molecular Depictions.

This module provides advanced depiction control features including:
- Zoom level control
- Stroke ratio (bond thickness) control
- SVG units control (px, mm, cm, in)
- Structure flipping
- Title display
- Anonymous atom display
- SMARTS hit limiting

Based on CDK DepictController.java advanced features.

"""

from __future__ import annotations

import logging
from enum import Enum
from typing import Any, List

from jpype import JClass

logger = logging.getLogger(__name__)


class SVGUnits(Enum):
    """Enumeration of supported SVG unit types.

    Attributes:
        PX: Pixels (default for web)
        MM: Millimeters
        CM: Centimeters
        IN: Inches
    """

    PX = "px"
    MM = "mm"
    CM = "cm"
    IN = "in"


class AdvancedControls:
    """Manages advanced depiction control parameters.

    This class provides methods to configure advanced rendering parameters
    for CDK molecular depictions.

    Attributes:
        StandardGenerator: CDK StandardGenerator class
        GeometryTools: CDK GeometryTools for transformations
        SymbolVisibility: CDK SymbolVisibility for anonymous display
    """

    def __init__(self):
        """Initialize the advanced controls system with CDK classes."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.StandardGenerator = JClass(
                self.cdk_base + ".renderer.generators.standard.StandardGenerator"
            )
            self.GeometryTools = JClass(self.cdk_base + ".geometry.GeometryTools")
            self.SymbolVisibility = JClass(self.cdk_base + ".renderer.SymbolVisibility")
            self.Math = JClass("java.lang.Math")
            logger.debug("AdvancedControls initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize AdvancedControls: {e}")
            raise

    def set_zoom(self, depiction_generator: Any, zoom: float = 1.3) -> Any:
        """Set zoom level for depiction.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            zoom: Zoom factor (0.1 to 5.0, default: 1.3)

        Returns:
            Updated DepictionGenerator with zoom applied

        Example:
            >>> controls = AdvancedControls()
            >>> gen = controls.set_zoom(gen, zoom=1.5)
        """
        try:
            # Clamp zoom to reasonable range
            zoom = max(0.1, min(5.0, zoom))

            depiction_generator = depiction_generator.withZoom(zoom)
            logger.debug(f"Set zoom level: {zoom}")
            return depiction_generator

        except Exception as e:
            logger.error(f"Failed to set zoom: {e}")
            return depiction_generator

    def set_stroke_ratio(self, depiction_generator: Any, ratio: float = 1.0) -> Any:
        """Set stroke ratio (bond line thickness).

        Args:
            depiction_generator: CDK DepictionGenerator instance
            ratio: Stroke ratio (0.5 to 2.0, default: 1.0)

        Returns:
            Updated DepictionGenerator with stroke ratio applied

        Example:
            >>> controls = AdvancedControls()
            >>> gen = controls.set_stroke_ratio(gen, ratio=1.2)
        """
        try:
            # Clamp ratio to reasonable range
            ratio = max(0.5, min(2.0, ratio))

            depiction_generator = depiction_generator.withParam(
                self.StandardGenerator.StrokeRatio.class_, ratio
            )
            logger.debug(f"Set stroke ratio: {ratio}")
            return depiction_generator

        except Exception as e:
            logger.error(f"Failed to set stroke ratio: {e}")
            return depiction_generator

    def flip_structure(self, molecule: Any) -> None:
        """Flip molecular structure horizontally.

        Performs a horizontal mirror flip and adjusts wedge/hash bonds accordingly.

        Args:
            molecule: CDK IAtomContainer

        Example:
            >>> controls = AdvancedControls()
            >>> controls.flip_structure(mol)
        """
        try:
            # Get 2D center
            center = self.GeometryTools.get2DCenter(molecule)

            # Flip all atoms
            for atom in molecule.atoms():
                point = atom.getPoint2d()
                if point is not None:
                    # Mirror x-coordinate across center
                    new_x = 2 * center.x - point.x
                    point.x = new_x

            # Flip stereo bonds (swap UP and DOWN)
            IBond = JClass(self.cdk_base + ".interfaces.IBond")
            for bond in molecule.bonds():
                stereo = bond.getStereo()
                if stereo == IBond.Stereo.UP:
                    bond.setStereo(IBond.Stereo.DOWN)
                elif stereo == IBond.Stereo.DOWN:
                    bond.setStereo(IBond.Stereo.UP)
                elif stereo == IBond.Stereo.UP_INVERTED:
                    bond.setStereo(IBond.Stereo.DOWN_INVERTED)
                elif stereo == IBond.Stereo.DOWN_INVERTED:
                    bond.setStereo(IBond.Stereo.UP_INVERTED)

            logger.debug("Flipped molecular structure")

        except Exception as e:
            logger.error(f"Failed to flip structure: {e}")

    def rotate_structure(self, molecule: Any, degrees: float) -> None:
        """Rotate molecular structure by specified degrees.

        Args:
            molecule: CDK IAtomContainer
            degrees: Rotation angle in degrees (clockwise)

        Example:
            >>> controls = AdvancedControls()
            >>> controls.rotate_structure(mol, degrees=45)
        """
        try:
            center = self.GeometryTools.get2DCenter(molecule)
            radians = degrees * self.Math.PI / 180.0
            self.GeometryTools.rotate(molecule, center, radians)
            logger.debug(f"Rotated structure by {degrees} degrees")
        except Exception as e:
            logger.error(f"Failed to rotate structure: {e}")

    def set_title_display(
        self,
        depiction_generator: Any,
        show_title: bool = True,
        is_reaction: bool = False,
    ) -> Any:
        """Enable or disable title display in depiction.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            show_title: Whether to show title (default: True)
            is_reaction: Whether depicting a reaction (default: False)

        Returns:
            Updated DepictionGenerator with title display configured

        Example:
            >>> controls = AdvancedControls()
            >>> gen = controls.set_title_display(gen, show_title=True)
        """
        try:
            if show_title:
                if is_reaction:
                    depiction_generator = depiction_generator.withRxnTitle()
                else:
                    depiction_generator = depiction_generator.withMolTitle()
                logger.debug("Enabled title display")

            return depiction_generator

        except Exception as e:
            logger.error(f"Failed to set title display: {e}")
            return depiction_generator

    def set_anonymous_display(
        self, depiction_generator: Any, anonymous: bool = False
    ) -> Any:
        """Set anonymous atom display mode.

        When enabled, only heteroatoms (non-C, non-H) are shown with labels.
        All carbon atoms appear as vertices without labels.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            anonymous: Whether to use anonymous display (default: False)

        Returns:
            Updated DepictionGenerator with anonymous display configured

        Example:
            >>> controls = AdvancedControls()
            >>> gen = controls.set_anonymous_display(gen, anonymous=True)
        """
        try:
            if anonymous:
                # Use iupacRecommendations which shows minimal labels
                visibility = self.SymbolVisibility.iupacRecommendations()
                depiction_generator = depiction_generator.withParam(
                    self.StandardGenerator.Visibility.class_, visibility
                )
                logger.debug("Enabled anonymous atom display")

            return depiction_generator

        except Exception as e:
            logger.error(f"Failed to set anonymous display: {e}")
            return depiction_generator


class SmartsHitLimiter:
    """Manages SMARTS pattern match limiting.

    Provides functionality to limit the number of SMARTS pattern matches
    that are highlighted in a depiction.
    """

    def __init__(self):
        """Initialize the SMARTS hit limiter."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.SmartsPattern = JClass(self.cdk_base + ".smarts.SmartsPattern")
            logger.debug("SmartsHitLimiter initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize SmartsHitLimiter: {e}")
            raise

    def limit_smarts_hits(
        self, molecule: Any, smarts: str, max_hits: int = 100
    ) -> List[List[int]]:
        """Find SMARTS matches with hit limiting.

        Args:
            molecule: CDK IAtomContainer
            smarts: SMARTS pattern string
            max_hits: Maximum number of matches to return (default: 100)

        Returns:
            List of matches, each match is a list of atom indices
            Limited to max_hits matches

        Example:
            >>> limiter = SmartsHitLimiter()
            >>> matches = limiter.limit_smarts_hits(mol, "c1ccccc1", max_hits=50)
        """
        matches = []

        try:
            pattern = self.SmartsPattern.create(smarts)

            # Find all matches
            all_matches = pattern.matchAll(molecule)

            # Limit number of matches
            count = 0
            for match in all_matches:
                if count >= max_hits:
                    logger.warning(f"SMARTS hit limit reached: {max_hits}")
                    break

                # Convert match to list of atom indices
                atom_indices = []
                for atom in match:
                    idx = molecule.indexOf(atom)
                    if idx >= 0:
                        atom_indices.append(idx)

                matches.append(atom_indices)
                count += 1

            logger.debug(f"Found {len(matches)} SMARTS matches (limit: {max_hits})")
            return matches

        except Exception as e:
            logger.error(f"Failed to limit SMARTS hits: {e}")
            return []


def get_svg_units(units_str: str) -> SVGUnits:
    """Convert string to SVGUnits enum.

    Args:
        units_str: String representation of units ("px", "mm", "cm", "in")

    Returns:
        Corresponding SVGUnits enum value

    Raises:
        ValueError: If units_str is not valid

    Example:
        >>> units = get_svg_units("mm")
        >>> assert units == SVGUnits.MM
    """
    try:
        return SVGUnits(units_str.lower())
    except ValueError:
        valid_units = [u.value for u in SVGUnits]
        raise ValueError(
            f"Invalid SVG units: '{units_str}'. "
            f"Valid units are: {', '.join(valid_units)}"
        )


# Convenience functions
def set_zoom(depiction_generator: Any, zoom: float = 1.3) -> Any:
    """Convenience function to set zoom level."""
    try:
        controls = AdvancedControls()
        return controls.set_zoom(depiction_generator, zoom)
    except Exception as e:
        logger.error(f"Failed to set zoom: {e}")
        return depiction_generator


def set_stroke_ratio(depiction_generator: Any, ratio: float = 1.0) -> Any:
    """Convenience function to set stroke ratio."""
    try:
        controls = AdvancedControls()
        return controls.set_stroke_ratio(depiction_generator, ratio)
    except Exception as e:
        logger.error(f"Failed to set stroke ratio: {e}")
        return depiction_generator


def flip_structure(molecule: Any) -> None:
    """Convenience function to flip structure."""
    try:
        controls = AdvancedControls()
        controls.flip_structure(molecule)
    except Exception as e:
        logger.error(f"Failed to flip structure: {e}")


def rotate_structure(molecule: Any, degrees: float) -> None:
    """Convenience function to rotate structure."""
    try:
        controls = AdvancedControls()
        controls.rotate_structure(molecule, degrees)
    except Exception as e:
        logger.error(f"Failed to rotate structure: {e}")


# Export public API
__all__ = [
    "SVGUnits",
    "AdvancedControls",
    "SmartsHitLimiter",
    "get_svg_units",
    "set_zoom",
    "set_stroke_ratio",
    "flip_structure",
    "rotate_structure",
]
