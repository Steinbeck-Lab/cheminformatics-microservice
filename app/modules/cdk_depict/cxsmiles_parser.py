"""CXSMILES Handler.

CXSMILES extends SMILES with additional information like:
- Highlighted atoms (ha:)
- Highlighted bonds (hb:)  
- 2D/3D coordinates
- Atom labels, values, etc.

This module provides utilities to parse and extract that information.
"""

from __future__ import annotations

import re
from typing import Any
from jpype import JClass


def parse_cxsmiles_highlighting_from_string(smiles_string: str) -> tuple[set, set]:
    """Parse CXSMILES highlighting directly from SMILES string.

    Since CDK's SmilesParser doesn't preserve CXSMILES highlighting in atom properties,
    we need to manually extract it from the input string.

    Args:
        smiles_string: CXSMILES string (e.g., "c1ccccc1 |ha:0,1,2|")

    Returns:
        Tuple of (highlighted_atom_indices, highlighted_bond_indices) as Python sets

    Examples:
        >>> parse_cxsmiles_highlighting_from_string("c1ccccc1 |ha:0,1,2|")
        ({0, 1, 2}, set())
        >>> parse_cxsmiles_highlighting_from_string("CCO |ha:0,1,hb:0|")
        ({0, 1}, {0})
    """
    highlighted_atoms = set()
    highlighted_bonds = set()

    if not smiles_string or "|" not in smiles_string:
        return highlighted_atoms, highlighted_bonds

    try:
        # Extract CXSMILES extension (everything after |)
        parts = smiles_string.split("|", 1)
        if len(parts) < 2:
            return highlighted_atoms, highlighted_bonds

        cxsmiles_ext = parts[1].rstrip("|")

        # Parse highlighted atoms: ha:0,1,2
        ha_match = re.search(r"ha:([0-9,]+)", cxsmiles_ext)
        if ha_match:
            atom_indices = ha_match.group(1).split(",")
            highlighted_atoms = {
                int(idx.strip()) for idx in atom_indices if idx.strip()
            }

        # Parse highlighted bonds: hb:0,1,2
        hb_match = re.search(r"hb:([0-9,]+)", cxsmiles_ext)
        if hb_match:
            bond_indices = hb_match.group(1).split(",")
            highlighted_bonds = {
                int(idx.strip()) for idx in bond_indices if idx.strip()
            }

    except Exception as e:
        # If parsing fails, return empty sets
        pass

    return highlighted_atoms, highlighted_bonds


def extract_cxsmiles_highlighting(molecule: Any) -> tuple[set, set]:
    """Extract CXSMILES highlighting from molecule.

    NOTE: CDK's SmilesParser does NOT preserve CXSMILES highlighting in atom/bond properties.
    Use parse_cxsmiles_highlighting_from_string() instead on the original SMILES string.

    This function is kept for backward compatibility but will return empty sets.

    Args:
        molecule: CDK IAtomContainer (already parsed from CXSMILES by SmilesParser)

    Returns:
        Tuple of (highlighted_atom_indices, highlighted_bond_indices) as Python sets
        (always empty - CDK doesn't preserve this information)
    """
    # CDK doesn't preserve CXSMILES highlighting in atom/bond properties
    # Use parse_cxsmiles_highlighting_from_string() instead
    return set(), set()


def apply_cxsmiles_highlighting_to_depiction(
    depiction_generator: Any, molecule: Any, highlight_color: Any = None
) -> Any:
    """Apply CXSMILES highlighting to depiction generator.

    Args:
        depiction_generator: CDK DepictionGenerator
        molecule: CDK IAtomContainer (parsed from CXSMILES)
        highlight_color: Java Color object (default: light green)

    Returns:
        Updated DepictionGenerator with highlighting applied
    """
    Color = JClass("java.awt.Color")

    if highlight_color is None:
        highlight_color = Color(0xAA, 0xFF, 0xAA)  # Light green

    try:
        # Extract highlighted atoms/bonds
        highlighted_atoms, highlighted_bonds = extract_cxsmiles_highlighting(molecule)

        # If we have any highlighting, apply it
        if not highlighted_atoms.isEmpty() or not highlighted_bonds.isEmpty():
            # Combine atoms and bonds into one set
            HashSet = JClass("java.util.HashSet")
            highlight_set = HashSet()

            # Add all highlighted atoms
            for atom in highlighted_atoms:
                highlight_set.add(atom)

            # Add all highlighted bonds
            for bond in highlighted_bonds:
                highlight_set.add(bond)

            # Apply highlighting to depiction generator
            depiction_generator = depiction_generator.withHighlight(
                highlight_set, highlight_color
            )

    except Exception as e:
        pass

    return depiction_generator


def parse_cxsmiles(cxsmiles: str, cdk_base: str = "org.openscience.cdk") -> Any:
    """Parse CXSMILES string and return molecule with highlighting.

    CDK's SmilesParser handles all CXSMILES extensions automatically.

    Args:
        cxsmiles: CXSMILES string (e.g., "CCO |ha:0,1|")
        cdk_base: Base CDK package path

    Returns:
        CDK IAtomContainer with highlighting information stored
    """
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())

    try:
        # CDK automatically parses CXSMILES extensions!
        molecule = SmilesParser.parseSmiles(cxsmiles)
        return molecule

    except Exception as e:
        raise


# Export public API
__all__ = [
    "parse_cxsmiles",
    "parse_cxsmiles_highlighting_from_string",
    "extract_cxsmiles_highlighting",
    "apply_cxsmiles_highlighting_to_depiction",
]
