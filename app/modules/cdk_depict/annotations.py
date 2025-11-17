"""Annotation System for CDK Molecular Depictions.

This module provides comprehensive annotation capabilities for molecular structures,
including atom numbering, bond numbering, atom mapping, atom values, and CIP stereochemistry.

Based on CDK DepictController.java annotation functionality.

Author: Kohulan Rajan
License: MIT
Date: 2024-2025
"""

from __future__ import annotations

import logging
from enum import Enum
from typing import Any
from jpype import JClass
from app.modules.toolkits.cdk_wrapper import get_cip_annotation

logger = logging.getLogger(__name__)


class AnnotationMode(Enum):
    """Enumeration of supported annotation display modes.

    Attributes:
        NONE: No annotations (default)
        NUMBER: Display atom numbers (0-based indexing)
        BONDNUMBER: Display bond numbers
        MAPIDX: Display atom mapping numbers (for reactions)
        ATOMVALUE: Display atom values/properties
        COLMAP: Color-code atoms by mapping number
        CIP: Display CIP stereochemistry labels (R/S, E/Z)
    """

    NONE = "none"
    NUMBER = "number"
    BONDNUMBER = "bondnumber"
    MAPIDX = "mapidx"
    ATOMVALUE = "atomvalue"
    COLMAP = "colmap"
    CIP = "cip"


# Color palette for atom mapping visualization (from CDK)
MAPPING_COLORS = [
    None,  # 0 - no mapping
    (0xFF, 0x7F, 0x7F),  # 1 - red
    (0xFF, 0xBF, 0x7F),  # 2 - orange
    (0xFF, 0xFF, 0x7F),  # 3 - yellow
    (0xBF, 0xFF, 0x7F),  # 4 - yellow-green
    (0x7F, 0xFF, 0x7F),  # 5 - green
    (0x7F, 0xFF, 0xBF),  # 6 - cyan-green
    (0x7F, 0xFF, 0xFF),  # 7 - cyan
    (0x7F, 0xBF, 0xFF),  # 8 - cyan-blue
    (0x7F, 0x7F, 0xFF),  # 9 - blue
    (0xBF, 0x7F, 0xFF),  # 10 - purple
    (0xFF, 0x7F, 0xFF),  # 11 - magenta
    (0xFF, 0x7F, 0xBF),  # 12 - pink
]


class AnnotationSystem:
    """Manages annotation modes for molecular depictions.

    This class provides methods to apply various annotation styles to CDK molecules
    and depiction generators, including atom numbering, bond numbering, atom mapping,
    and stereochemistry labels.

    Attributes:
        StandardGenerator: CDK StandardGenerator class
        Color: Java Color class
        CDKConstants: CDK constants
    """

    def __init__(self):
        """Initialize the annotation system with CDK classes."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.StandardGenerator = JClass(
                self.cdk_base + ".renderer.generators.standard.StandardGenerator"
            )
            self.Color = JClass("java.awt.Color")
            self.CDKConstants = JClass(self.cdk_base + ".CDKConstants")
            logger.debug("AnnotationSystem initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize AnnotationSystem: {e}")
            raise

    def apply_annotations(
        self,
        depiction_generator: Any,
        molecule: Any,
        mode: AnnotationMode,
        is_reaction: bool = False,
    ) -> Any:
        """Apply annotation mode to depiction generator and/or molecule.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            molecule: CDK IAtomContainer or IReaction
            mode: Annotation mode to apply
            is_reaction: Whether the molecule is a reaction

        Returns:
            Updated DepictionGenerator with annotations applied

        Example:
            >>> system = AnnotationSystem()
            >>> gen = get_depiction_generator()
            >>> gen = system.apply_annotations(gen, mol, AnnotationMode.NUMBER)
        """
        try:
            if mode == AnnotationMode.NONE:
                logger.debug("No annotations to apply")
                return depiction_generator

            if mode == AnnotationMode.NUMBER:
                return self._apply_atom_numbers(depiction_generator)

            elif mode == AnnotationMode.BONDNUMBER:
                return self._apply_bond_numbers(depiction_generator)

            elif mode == AnnotationMode.MAPIDX:
                return self._apply_atom_mapping_numbers(depiction_generator)

            elif mode == AnnotationMode.ATOMVALUE:
                self._annotate_atom_values(molecule)
                return depiction_generator

            elif mode == AnnotationMode.COLMAP:
                self._annotate_color_mapping(molecule, is_reaction)
                return depiction_generator

            elif mode == AnnotationMode.CIP:
                self._annotate_cip(molecule, is_reaction)
                return depiction_generator

            else:
                logger.warning(f"Unknown annotation mode: {mode}")
                return depiction_generator

        except Exception as e:
            logger.error(f"Failed to apply annotations ({mode}): {e}")
            # Return unmodified generator on error
            return depiction_generator

    def _apply_atom_numbers(self, depiction_generator: Any) -> Any:
        """Apply atom numbering to depiction.

        Args:
            depiction_generator: CDK DepictionGenerator

        Returns:
            Updated DepictionGenerator with atom numbers enabled
        """
        try:
            logger.debug("Applying atom numbers")
            return depiction_generator.withAtomNumbers()
        except Exception as e:
            logger.error(f"Failed to apply atom numbers: {e}")
            return depiction_generator

    def _apply_bond_numbers(self, depiction_generator: Any) -> Any:
        """Apply bond numbering to depiction.

        Args:
            depiction_generator: CDK DepictionGenerator

        Returns:
            Updated DepictionGenerator with bond numbers enabled
        """
        try:
            logger.debug("Applying bond numbers")
            # Bond numbers are set via StandardGenerator visibility
            visibility = self.StandardGenerator.Visibility.BOND_NUMBERS
            return depiction_generator.withParam(
                self.StandardGenerator.Visibility.class_, visibility
            )
        except Exception as e:
            logger.error(f"Failed to apply bond numbers: {e}")
            return depiction_generator

    def _apply_atom_mapping_numbers(self, depiction_generator: Any) -> Any:
        """Apply atom mapping numbers to depiction.

        Args:
            depiction_generator: CDK DepictionGenerator

        Returns:
            Updated DepictionGenerator with mapping numbers enabled
        """
        try:
            logger.debug("Applying atom mapping numbers")
            return depiction_generator.withAtomMapNumbers()
        except Exception as e:
            logger.error(f"Failed to apply atom mapping numbers: {e}")
            return depiction_generator

    def _annotate_atom_values(self, molecule: Any) -> None:
        """Annotate atoms with their values/properties.

        Sets StandardGenerator.ANNOTATION_LABEL property on atoms that have
        values defined in CDK's atom value property.

        Args:
            molecule: CDK IAtomContainer or IReaction
        """
        try:
            logger.debug("Annotating atom values")

            # Handle reactions
            if self._is_reaction(molecule):
                ReactionManipulator = JClass(
                    self.cdk_base + ".reaction.ReactionManipulator"
                )
                containers = ReactionManipulator.getAllAtomContainers(molecule)
                for container in containers:
                    self._annotate_atom_values_container(container)
            else:
                self._annotate_atom_values_container(molecule)

        except Exception as e:
            logger.error(f"Failed to annotate atom values: {e}")

    def _annotate_atom_values_container(self, container: Any) -> None:
        """Annotate atom values for a single container.

        Args:
            container: CDK IAtomContainer
        """
        for atom in container.atoms():
            atom_value = atom.getProperty(self.CDKConstants.COMMENT)
            if atom_value is not None:
                atom.setProperty(
                    self.StandardGenerator.ANNOTATION_LABEL, str(atom_value)
                )
                logger.debug(f"Set atom value annotation: {atom_value}")

    def _annotate_color_mapping(self, molecule: Any, is_reaction: bool) -> None:
        """Color-code atoms by their mapping numbers.

        Args:
            molecule: CDK IAtomContainer or IReaction
            is_reaction: Whether this is a reaction
        """
        try:
            logger.debug("Applying color mapping")

            if is_reaction:
                ReactionManipulator = JClass(
                    self.cdk_base + ".reaction.ReactionManipulator"
                )
                containers = ReactionManipulator.getAllAtomContainers(molecule)
                for container in containers:
                    self._color_map_container(container)
            else:
                self._color_map_container(molecule)

        except Exception as e:
            logger.error(f"Failed to apply color mapping: {e}")

    def _color_map_container(self, container: Any) -> None:
        """Apply color mapping to a single container.

        Args:
            container: CDK IAtomContainer
        """
        for atom in container.atoms():
            mapidx = atom.getProperty(self.CDKConstants.ATOM_ATOM_MAPPING)
            if mapidx is not None and 0 < mapidx < len(MAPPING_COLORS):
                color_tuple = MAPPING_COLORS[mapidx]
                if color_tuple:
                    r, g, b = color_tuple
                    color = self.Color(r, g, b)
                    atom.setProperty(self.StandardGenerator.HIGHLIGHT_COLOR, color)
                    logger.debug(f"Set color for atom mapping {mapidx}")

    def _annotate_cip(self, molecule: Any, is_reaction: bool) -> None:
        """Annotate CIP stereochemistry labels (R/S, E/Z).

        Args:
            molecule: CDK IAtomContainer or IReaction
            is_reaction: Whether this is a reaction
        """
        try:
            logger.debug("Annotating CIP stereochemistry")

            if is_reaction:
                ReactionManipulator = JClass(
                    self.cdk_base + ".reaction.ReactionManipulator"
                )
                containers = ReactionManipulator.getAllAtomContainers(molecule)
                for container in containers:
                    get_cip_annotation(container)
            else:
                get_cip_annotation(molecule)

        except Exception as e:
            logger.error(f"Failed to annotate CIP: {e}")

    def _is_reaction(self, molecule: Any) -> bool:
        """Check if the molecule is a reaction.

        Args:
            molecule: CDK IChemObject

        Returns:
            True if molecule is an IReaction or IReactionSet
        """
        try:
            IReaction = JClass(self.cdk_base + ".interfaces.IReaction")
            IReactionSet = JClass(self.cdk_base + ".interfaces.IReactionSet")
            return isinstance(molecule, (IReaction, IReactionSet))
        except Exception:
            return False


def get_annotation_mode(mode_str: str) -> AnnotationMode:
    """Convert string to AnnotationMode enum.

    Args:
        mode_str: String representation of annotation mode

    Returns:
        Corresponding AnnotationMode enum value

    Raises:
        ValueError: If mode_str is not a valid annotation mode

    Example:
        >>> mode = get_annotation_mode("number")
        >>> assert mode == AnnotationMode.NUMBER
    """
    try:
        return AnnotationMode(mode_str.strip().lower())
    except ValueError:
        valid_modes = [m.value for m in AnnotationMode]
        raise ValueError(
            f"Invalid annotation mode: '{mode_str}'. "
            f"Valid modes are: {', '.join(valid_modes)}"
        )


def apply_annotations(
    depiction_generator: Any,
    molecule: Any,
    mode: str = "none",
    is_reaction: bool = False,
) -> Any:
    """Convenience function to apply annotations to a depiction.

    Args:
        depiction_generator: CDK DepictionGenerator instance
        molecule: CDK IAtomContainer or IReaction
        mode: Annotation mode as string (default: "none")
        is_reaction: Whether the molecule is a reaction

    Returns:
        Updated DepictionGenerator with annotations applied

    Example:
        >>> gen = apply_annotations(gen, mol, mode="number")
        >>> gen = apply_annotations(gen, rxn, mode="mapidx", is_reaction=True)
    """
    try:
        annotation_mode = get_annotation_mode(mode)
        system = AnnotationSystem()
        return system.apply_annotations(
            depiction_generator, molecule, annotation_mode, is_reaction
        )
    except Exception as e:
        logger.error(f"Failed to apply annotations: {e}")
        return depiction_generator


# Export public API
__all__ = [
    "AnnotationMode",
    "AnnotationSystem",
    "MAPPING_COLORS",
    "get_annotation_mode",
    "apply_annotations",
]
