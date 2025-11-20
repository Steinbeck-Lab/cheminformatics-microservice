"""Aromatic Display (Donuts) for CDK Molecular Depictions.

This module provides functionality to display aromatic rings with a circle inside
the ring (donut representation) instead of alternating double bonds.

Based on CDK DepictController.java aromatic display functionality.
"""

from __future__ import annotations

import logging
from typing import Any, List, Set

from jpype import JClass

logger = logging.getLogger(__name__)


class AromaticDisplaySystem:
    """Manages aromatic ring display modes for molecular depictions."""

    def __init__(self):
        """Initialize the aromatic display system with CDK classes."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.Aromaticity = JClass(self.cdk_base + ".aromaticity.Aromaticity")
            self.ElectronDonation = JClass(
                self.cdk_base + ".aromaticity.ElectronDonation"
            )
            self.Cycles = JClass(self.cdk_base + ".graph.Cycles")
            logger.debug("AromaticDisplaySystem initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize AromaticDisplaySystem: {e}")
            raise

    def apply_donut_display(
        self, depiction_generator: Any, molecule: Any, enable: bool = True
    ) -> Any:
        """Apply or disable donut (circle-in-ring) display for aromatic rings.

        Args:
            depiction_generator: CDK DepictionGenerator instance
            molecule: CDK IAtomContainer
            enable: Whether to enable (True) or disable (False) donut display

        Returns:
            Updated DepictionGenerator with aromatic display configured
        """
        try:
            if enable:
                # Perceive aromaticity first
                self._perceive_aromaticity(molecule)

                # Enable donut display - simple method call
                depiction_generator = depiction_generator.withAromaticDisplay()
                logger.debug("Enabled donut (circle-in-ring) aromatic display")
            else:
                # Donut display is off by default, no action needed
                logger.debug("Donut display disabled (default Kekule structure)")

            return depiction_generator

        except Exception as e:
            logger.error(f"Failed to apply donut display: {e}")
            return depiction_generator

    def _perceive_aromaticity(self, molecule: Any) -> None:
        """Perceive aromaticity in the molecule."""
        try:
            aromaticity = self.Aromaticity(
                self.ElectronDonation.daylight(),
                self.Cycles.or_(self.Cycles.all(), self.Cycles.relevant()),
            )
            aromaticity.apply(molecule)
            logger.debug("Aromaticity perception completed")
        except Exception as e:
            logger.warning(f"Failed to perceive aromaticity: {e}")

    def get_aromatic_rings(self, molecule: Any) -> List[Set[int]]:
        """Get list of aromatic ring systems."""
        aromatic_rings = []

        try:
            self._perceive_aromaticity(molecule)

            AllRingsFinder = JClass(self.cdk_base + ".ringsearch.AllRingsFinder")
            ringFinder = AllRingsFinder()
            ringSet = ringFinder.findAllRings(molecule)

            for i in range(ringSet.getAtomContainerCount()):
                ring = ringSet.getAtomContainer(i)
                if self._is_aromatic_ring(ring):
                    atom_indices = set()
                    for atom in ring.atoms():
                        idx = molecule.indexOf(atom)
                        if idx >= 0:
                            atom_indices.add(idx)

                    if atom_indices:
                        aromatic_rings.append(atom_indices)

            logger.debug(f"Found {len(aromatic_rings)} aromatic rings")
            return aromatic_rings

        except Exception as e:
            logger.error(f"Failed to get aromatic rings: {e}")
            return []

    def _is_aromatic_ring(self, ring: Any) -> bool:
        """Check if a ring is aromatic."""
        try:
            for atom in ring.atoms():
                if not atom.isAromatic():
                    return False

            for bond in ring.bonds():
                if not bond.isAromatic():
                    return False

            return True
        except Exception:
            return False

    def count_aromatic_atoms(self, molecule: Any) -> int:
        """Count number of aromatic atoms in molecule."""
        try:
            self._perceive_aromaticity(molecule)

            count = 0
            for atom in molecule.atoms():
                if atom.isAromatic():
                    count += 1

            logger.debug(f"Found {count} aromatic atoms")
            return count

        except Exception as e:
            logger.error(f"Failed to count aromatic atoms: {e}")
            return 0


def apply_aromatic_display(
    depiction_generator: Any, molecule: Any, donuts: bool = True
) -> Any:
    """Convenience function to apply aromatic display mode.

    Args:
        depiction_generator: CDK DepictionGenerator instance
        molecule: CDK IAtomContainer
        donuts: Whether to use donut (circle-in-ring) display (default: True)

    Returns:
        Updated DepictionGenerator with aromatic display configured
    """
    try:
        system = AromaticDisplaySystem()
        return system.apply_donut_display(depiction_generator, molecule, donuts)
    except Exception as e:
        logger.error(f"Failed to apply aromatic display: {e}")
        return depiction_generator


def perceive_aromaticity(molecule: Any) -> None:
    """Convenience function to perceive aromaticity."""
    try:
        system = AromaticDisplaySystem()
        system._perceive_aromaticity(molecule)
    except Exception as e:
        logger.error(f"Failed to perceive aromaticity: {e}")


__all__ = [
    "AromaticDisplaySystem",
    "apply_aromatic_display",
    "perceive_aromaticity",
]
