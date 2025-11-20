"""Reaction Arrow Types for CDK Molecular Depictions.

This module provides support for different reaction arrow types including
forward, bidirectional, no-go, retrosynthetic, and resonance arrows.

Based on CDK DepictController.java reaction arrow functionality.

"""

from __future__ import annotations

import logging
from enum import Enum
from typing import Any, Optional

from jpype import JClass

logger = logging.getLogger(__name__)


class ReactionArrowType(Enum):
    """Enumeration of supported reaction arrow types.

    Attributes:
        FORWARD: Standard forward arrow (→)
        BIDIRECTIONAL: Equilibrium/reversible arrow (⇌)
        NO_GO: No-go/blocked arrow (⇏)
        RETRO_SYNTHETIC: Retrosynthetic arrow (⇒)
        RESONANCE: Resonance arrow (↔)
    """

    FORWARD = "forward"
    BIDIRECTIONAL = "bidirectional"
    NO_GO = "no_go"
    RETRO_SYNTHETIC = "retrosynthetic"
    RESONANCE = "resonance"


# Map string abbreviations to arrow types (matching CDK Java)
ARROW_ABBREV_MAP = {
    "": ReactionArrowType.FORWARD,  # Default
    "fwd": ReactionArrowType.FORWARD,
    "forward": ReactionArrowType.FORWARD,
    "equ": ReactionArrowType.BIDIRECTIONAL,
    "bidirectional": ReactionArrowType.BIDIRECTIONAL,
    "equilibrium": ReactionArrowType.BIDIRECTIONAL,
    "ngo": ReactionArrowType.NO_GO,
    "no_go": ReactionArrowType.NO_GO,
    "ret": ReactionArrowType.RETRO_SYNTHETIC,
    "retro": ReactionArrowType.RETRO_SYNTHETIC,
    "retrosynthetic": ReactionArrowType.RETRO_SYNTHETIC,
    "res": ReactionArrowType.RESONANCE,
    "resonance": ReactionArrowType.RESONANCE,
}


class ReactionArrowSystem:
    """Manages reaction arrow types for reaction depictions.

    This class provides methods to set different arrow types for chemical reactions
    in CDK depictions.

    Attributes:
        IReaction: CDK IReaction interface
        Direction: CDK IReaction.Direction enum
    """

    def __init__(self):
        """Initialize the reaction arrow system with CDK classes."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.IReaction = JClass(self.cdk_base + ".interfaces.IReaction")
            # Get Direction enum from IReaction
            self.Direction = self.IReaction.Direction
            logger.debug("ReactionArrowSystem initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize ReactionArrowSystem: {e}")
            raise

    def set_arrow_type(self, reaction: Any, arrow_type: ReactionArrowType) -> None:
        """Set the arrow type for a reaction.

        Args:
            reaction: CDK IReaction or IReactionSet
            arrow_type: Type of arrow to use

        Example:
            >>> system = ReactionArrowSystem()
            >>> system.set_arrow_type(reaction, ReactionArrowType.BIDIRECTIONAL)
        """
        try:
            # Map arrow type to CDK Direction
            direction = self._get_cdk_direction(arrow_type)

            # Handle reaction sets
            if self._is_reaction_set(reaction):
                for i in range(reaction.getReactionCount()):
                    rxn = reaction.getReaction(i)
                    rxn.setDirection(direction)
                logger.debug(
                    f"Set arrow type for {reaction.getReactionCount()} reactions"
                )
            else:
                reaction.setDirection(direction)
                logger.debug(f"Set arrow type: {arrow_type.value}")

        except Exception as e:
            logger.error(f"Failed to set arrow type: {e}")

    def _get_cdk_direction(self, arrow_type: ReactionArrowType) -> Any:
        """Map ReactionArrowType to CDK Direction enum.

        Args:
            arrow_type: ReactionArrowType enum value

        Returns:
            CDK IReaction.Direction enum value
        """
        # Map to CDK Direction enum
        mapping = {
            ReactionArrowType.FORWARD: self.Direction.FORWARD,
            ReactionArrowType.BIDIRECTIONAL: self.Direction.BIDIRECTIONAL,
            ReactionArrowType.NO_GO: self.Direction.NO_GO,
            ReactionArrowType.RETRO_SYNTHETIC: self.Direction.RETRO,
            ReactionArrowType.RESONANCE: self.Direction.BIDIRECTIONAL,  # Use bidirectional for resonance
        }
        return mapping.get(arrow_type, self.Direction.FORWARD)

    def _is_reaction_set(self, obj: Any) -> bool:
        """Check if object is a reaction set.

        Args:
            obj: Object to check

        Returns:
            True if obj is an IReactionSet
        """
        try:
            IReactionSet = JClass(self.cdk_base + ".interfaces.IReactionSet")
            return isinstance(obj, IReactionSet)
        except Exception:
            return False

    def get_arrow_type(self, reaction: Any) -> Optional[ReactionArrowType]:
        """Get the current arrow type of a reaction.

        Args:
            reaction: CDK IReaction

        Returns:
            Current ReactionArrowType or None if not set
        """
        try:
            direction = reaction.getDirection()

            # Map CDK Direction back to ReactionArrowType
            if direction == self.Direction.FORWARD:
                return ReactionArrowType.FORWARD
            elif direction == self.Direction.BIDIRECTIONAL:
                return ReactionArrowType.BIDIRECTIONAL
            elif direction == self.Direction.NO_GO:
                return ReactionArrowType.NO_GO
            elif direction == self.Direction.RETRO:
                return ReactionArrowType.RETRO_SYNTHETIC
            else:
                return None

        except Exception as e:
            logger.error(f"Failed to get arrow type: {e}")
            return None


class ReactionLayoutSystem:
    """Manages reaction layout and alignment settings.

    This class provides methods to control reaction layout parameters such as
    atom map alignment and reaction scheme grid layouts.

    Attributes:
        DepictionGenerator: CDK DepictionGenerator class
    """

    def __init__(self):
        """Initialize the reaction layout system."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.DepictionGenerator = JClass(
                self.cdk_base + ".depict.DepictionGenerator"
            )
            logger.debug("ReactionLayoutSystem initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize ReactionLayoutSystem: {e}")
            raise

    def set_map_alignment(self, depiction_generator: Any, align: bool = True) -> Any:
        """Enable or disable atom map alignment in reactions.

        When enabled, atoms with the same mapping number will be vertically aligned
        across reactants and products for clearer visualization.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            align: Whether to align mapped atoms (default: True)

        Returns:
            Updated DepictionGenerator with alignment configured

        Example:
            >>> system = ReactionLayoutSystem()
            >>> gen = system.set_map_alignment(gen, align=True)
        """
        try:
            depiction_generator = depiction_generator.withMappedRxnAlign(align)
            logger.debug(f"Set reaction map alignment: {align}")
            return depiction_generator
        except Exception as e:
            logger.error(f"Failed to set map alignment: {e}")
            return depiction_generator


def get_arrow_type(arrow_str: str) -> ReactionArrowType:
    """Convert string to ReactionArrowType enum.

    Supports both full names and abbreviations:
    - "forward", "fwd", "" -> FORWARD
    - "bidirectional", "equ", "equilibrium" -> BIDIRECTIONAL
    - "no_go", "ngo" -> NO_GO
    - "retrosynthetic", "ret", "retro" -> RETRO_SYNTHETIC
    - "resonance", "res" -> RESONANCE

    Args:
        arrow_str: String representation of arrow type

    Returns:
        Corresponding ReactionArrowType enum value

    Raises:
        ValueError: If arrow_str is not a valid arrow type

    Example:
        >>> arrow = get_arrow_type("equ")
        >>> assert arrow == ReactionArrowType.BIDIRECTIONAL
    """
    arrow_lower = arrow_str.lower().strip()

    if arrow_lower in ARROW_ABBREV_MAP:
        return ARROW_ABBREV_MAP[arrow_lower]

    # Try direct enum lookup
    try:
        return ReactionArrowType(arrow_lower)
    except ValueError:
        valid_types = list(ARROW_ABBREV_MAP.keys())
        raise ValueError(
            f"Invalid arrow type: '{arrow_str}'. "
            f"Valid types are: {', '.join(valid_types)}"
        )


def set_reaction_arrow(reaction: Any, arrow: str = "forward") -> None:
    """Convenience function to set reaction arrow type.

    Args:
        reaction: CDK IReaction or IReactionSet
        arrow: Arrow type as string (default: "forward")

    Example:
        >>> set_reaction_arrow(reaction, arrow="equ")
        >>> set_reaction_arrow(reaction, arrow="retro")
    """
    try:
        arrow_type = get_arrow_type(arrow)
        system = ReactionArrowSystem()
        system.set_arrow_type(reaction, arrow_type)
    except Exception as e:
        logger.error(f"Failed to set reaction arrow: {e}")


def set_map_alignment(depiction_generator: Any, align_rxnmap: bool = True) -> Any:
    """Convenience function to set reaction map alignment.

    Args:
        depiction_generator: CDK DepictionGenerator instance
        align_rxnmap: Whether to align mapped atoms (default: True)

    Returns:
        Updated DepictionGenerator

    Example:
        >>> gen = set_map_alignment(gen, align_rxnmap=True)
    """
    try:
        system = ReactionLayoutSystem()
        return system.set_map_alignment(depiction_generator, align_rxnmap)
    except Exception as e:
        logger.error(f"Failed to set map alignment: {e}")
        return depiction_generator


# Export public API
__all__ = [
    "ReactionArrowType",
    "ReactionArrowSystem",
    "ReactionLayoutSystem",
    "ARROW_ABBREV_MAP",
    "get_arrow_type",
    "set_reaction_arrow",
    "set_map_alignment",
]
