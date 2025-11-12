"""MDL HILITE Support for CDK Molecular Depictions.

This module provides functionality to parse and apply MDL/SDF HILITE properties
for atom and bond highlighting in V3000 format molfiles.

Based on CDK MDL HILITE parsing functionality.

"""

from __future__ import annotations

import logging
import re
from typing import Any, Set, Tuple, Optional, List

from jpype import JClass

logger = logging.getLogger(__name__)


class MDLHiliteParser:
    """Parses MDL HILITE properties for atom/bond highlighting.

    This class provides methods to extract highlight information from MDL V3000
    format molfiles and apply highlighting to CDK depictions.

    HILITE format in V3000:
        M  V30 HIGHLIGHT ATOMS=(count index1 index2 ...)
        M  V30 HIGHLIGHT BONDS=(count index1 index2 ...)

    Attributes:
        Color: Java Color class
        StandardGenerator: CDK StandardGenerator class
    """

    def __init__(self):
        """Initialize the MDL HILITE parser with CDK classes."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.Color = JClass("java.awt.Color")
            self.StandardGenerator = JClass(
                self.cdk_base + ".renderer.generators.standard.StandardGenerator"
            )
            self.CDKConstants = JClass(self.cdk_base + ".CDKConstants")
            logger.debug("MDLHiliteParser initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize MDLHiliteParser: {e}")
            raise

    def parse_hilite_from_mol(self, molecule: Any) -> Tuple[Set[int], Set[int]]:
        """Extract HILITE information from molecule properties.

        Args:
            molecule: CDK IAtomContainer with HILITE properties

        Returns:
            Tuple of (highlighted_atoms, highlighted_bonds)
            where each is a set of 0-based indices

        Example:
            >>> parser = MDLHiliteParser()
            >>> atoms, bonds = parser.parse_hilite_from_mol(mol)
            >>> print(f"Highlighted {len(atoms)} atoms and {len(bonds)} bonds")
        """
        highlighted_atoms = set()
        highlighted_bonds = set()

        try:
            # Check for HILITE property in molecule
            hilite_atoms = molecule.getProperty("HILITE_ATOMS")
            hilite_bonds = molecule.getProperty("HILITE_BONDS")

            if hilite_atoms:
                highlighted_atoms = self._parse_hilite_indices(hilite_atoms)
                logger.debug(
                    f"Parsed {len(highlighted_atoms)} highlighted atoms from property"
                )

            if hilite_bonds:
                highlighted_bonds = self._parse_hilite_indices(hilite_bonds)
                logger.debug(
                    f"Parsed {len(highlighted_bonds)} highlighted bonds from property"
                )

            return highlighted_atoms, highlighted_bonds

        except Exception as e:
            logger.error(f"Failed to parse HILITE from molecule: {e}")
            return set(), set()

    def parse_hilite_from_mdl(self, mdl_text: str) -> Tuple[Set[int], Set[int]]:
        """Parse HILITE information from MDL V3000 text.

        Args:
            mdl_text: MDL V3000 format molfile text

        Returns:
            Tuple of (highlighted_atoms, highlighted_bonds)

        Example:
            >>> parser = MDLHiliteParser()
            >>> atoms, bonds = parser.parse_hilite_from_mdl(mdl_content)
        """
        highlighted_atoms = set()
        highlighted_bonds = set()

        try:
            # Look for V30 HIGHLIGHT lines
            for line in mdl_text.split("\n"):
                if "V30 HIGHLIGHT" in line:
                    if "ATOMS" in line:
                        atoms = self._extract_indices_from_line(line)
                        highlighted_atoms.update(atoms)
                        logger.debug(f"Parsed {len(atoms)} atoms from HILITE line")
                    elif "BONDS" in line:
                        bonds = self._extract_indices_from_line(line)
                        highlighted_bonds.update(bonds)
                        logger.debug(f"Parsed {len(bonds)} bonds from HILITE line")

            return highlighted_atoms, highlighted_bonds

        except Exception as e:
            logger.error(f"Failed to parse HILITE from MDL: {e}")
            return set(), set()

    def _parse_hilite_indices(self, hilite_str: Any) -> Set[int]:
        """Parse HILITE indices from property string.

        Args:
            hilite_str: HILITE property value (string or array)

        Returns:
            Set of 0-based atom/bond indices
        """
        indices = set()

        try:
            # Convert to string if needed
            if not isinstance(hilite_str, str):
                hilite_str = str(hilite_str)

            # Parse comma-separated or space-separated indices
            parts = re.split(r"[,\s]+", hilite_str.strip())
            for part in parts:
                try:
                    idx = int(part)
                    # Convert from 1-based to 0-based indexing
                    indices.add(idx - 1)
                except ValueError:
                    continue

            return indices

        except Exception as e:
            logger.warning(f"Failed to parse HILITE indices: {e}")
            return set()

    def _extract_indices_from_line(self, line: str) -> Set[int]:
        """Extract indices from V30 HIGHLIGHT line.

        Format: M  V30 HIGHLIGHT ATOMS=(count index1 index2 ...)

        Args:
            line: V30 HIGHLIGHT line from molfile

        Returns:
            Set of 0-based indices
        """
        indices = set()

        try:
            # Find content within parentheses
            match = re.search(r"\(([^)]+)\)", line)
            if match:
                content = match.group(1)
                parts = content.split()

                # First number is count, rest are indices
                if len(parts) > 1:
                    for idx_str in parts[1:]:
                        try:
                            idx = int(idx_str)
                            # Convert from 1-based to 0-based
                            indices.add(idx - 1)
                        except ValueError:
                            continue

            return indices

        except Exception as e:
            logger.warning(f"Failed to extract indices from line: {e}")
            return set()

    def apply_hilite_to_molecule(
        self,
        molecule: Any,
        highlighted_atoms: Set[int],
        highlighted_bonds: Set[int],
        highlight_color: Optional[Tuple[int, int, int]] = None,
    ) -> None:
        """Apply highlighting to atoms and bonds in molecule.

        Sets HIGHLIGHT_COLOR property on specified atoms and bonds.

        Args:
            molecule: CDK IAtomContainer
            highlighted_atoms: Set of atom indices to highlight
            highlighted_bonds: Set of bond indices to highlight
            highlight_color: RGB tuple for highlight color (default: light green)

        Example:
            >>> parser = MDLHiliteParser()
            >>> parser.apply_hilite_to_molecule(mol, {0, 1, 2}, {0, 1})
        """
        try:
            # Default highlight color (light green)
            if highlight_color is None:
                highlight_color = (170, 255, 170)

            r, g, b = highlight_color
            color = self.Color(r, g, b)

            # Highlight atoms
            for atom_idx in highlighted_atoms:
                if 0 <= atom_idx < molecule.getAtomCount():
                    atom = molecule.getAtom(atom_idx)
                    atom.setProperty(self.StandardGenerator.HIGHLIGHT_COLOR, color)
                    logger.debug(f"Highlighted atom {atom_idx}")

            # Highlight bonds
            for bond_idx in highlighted_bonds:
                if 0 <= bond_idx < molecule.getBondCount():
                    bond = molecule.getBond(bond_idx)
                    bond.setProperty(self.StandardGenerator.HIGHLIGHT_COLOR, color)
                    logger.debug(f"Highlighted bond {bond_idx}")

        except Exception as e:
            logger.error(f"Failed to apply highlighting: {e}")

    def has_hilite_properties(self, molecule: Any) -> bool:
        """Check if molecule has HILITE properties.

        Args:
            molecule: CDK IAtomContainer

        Returns:
            True if HILITE properties are present
        """
        try:
            return (
                molecule.getProperty("HILITE_ATOMS") is not None
                or molecule.getProperty("HILITE_BONDS") is not None
            )
        except:
            return False


def parse_mdl_hilite(molecule: Any) -> Tuple[Set[int], Set[int]]:
    """Convenience function to parse HILITE from molecule.

    Args:
        molecule: CDK IAtomContainer

    Returns:
        Tuple of (highlighted_atoms, highlighted_bonds)

    Example:
        >>> atoms, bonds = parse_mdl_hilite(mol)
    """
    try:
        parser = MDLHiliteParser()
        return parser.parse_hilite_from_mol(molecule)
    except Exception as e:
        logger.error(f"Failed to parse MDL HILITE: {e}")
        return set(), set()


def apply_mdl_hilite(
    molecule: Any, highlight_color: Optional[Tuple[int, int, int]] = None
) -> None:
    """Convenience function to parse and apply HILITE from molecule.

    Args:
        molecule: CDK IAtomContainer with HILITE properties
        highlight_color: Optional RGB tuple for highlight color

    Example:
        >>> apply_mdl_hilite(mol)
        >>> apply_mdl_hilite(mol, highlight_color=(255, 200, 200))
    """
    try:
        parser = MDLHiliteParser()
        atoms, bonds = parser.parse_hilite_from_mol(molecule)
        if atoms or bonds:
            parser.apply_hilite_to_molecule(molecule, atoms, bonds, highlight_color)
    except Exception as e:
        logger.error(f"Failed to apply MDL HILITE: {e}")


# Export public API
__all__ = [
    "MDLHiliteParser",
    "parse_mdl_hilite",
    "apply_mdl_hilite",
]
