"""Radical Perception for CDK Molecular Depictions.

This module provides functionality to detect and display radical electrons
(unpaired electrons) on atoms in molecular structures.

Based on CDK radical handling functionality.

"""

from __future__ import annotations

from typing import Any, List, Dict, Optional

from jpype import JClass


class RadicalPerception:
    """Manages radical detection and display for molecular depictions.

    This class provides methods to detect unpaired electrons on atoms and
    configure their display in CDK depictions.

    Attributes:
        SPIN_MULTIPLICITY: CDK constant for spin multiplicity property
        IAtom: CDK IAtom interface
        Integer: Java Integer class
    """

    def __init__(self):
        """Initialize the radical perception system with CDK classes."""
        self.cdk_base = "org.openscience.cdk"
        try:
            self.CDKConstants = JClass(self.cdk_base + ".CDKConstants")
            self.SPIN_MULTIPLICITY = self.CDKConstants.SPIN_MULTIPLICITY
            self.IAtom = JClass(self.cdk_base + ".interfaces.IAtom")
            self.Integer = JClass("java.lang.Integer")
        except Exception:
            raise

    def perceive_radicals(self, molecule: Any) -> Dict[int, int]:
        """Detect radical electrons in molecule.

        Analyzes valence, formal charge, and connectivity to identify atoms
        with unpaired electrons.

        Args:
            molecule: CDK IAtomContainer

        Returns:
            Dictionary mapping atom index to number of unpaired electrons

        Example:
            >>> perceiver = RadicalPerception()
            >>> radicals = perceiver.perceive_radicals(mol)
            >>> for atom_idx, radical_count in radicals.items():
            ...     print(f"Atom {atom_idx} has {radical_count} radical(s)")
        """
        radicals = {}

        try:
            for i, atom in enumerate(molecule.atoms()):
                radical_count = self._get_radical_count(atom, molecule)
                if radical_count > 0:
                    radicals[i] = radical_count

            return radicals

        except Exception:
            return {}

    def _get_radical_count(self, atom: Any, molecule: Any) -> int:
        """Calculate number of unpaired electrons on an atom.

        Args:
            atom: CDK IAtom
            molecule: CDK IAtomContainer (parent molecule)

        Returns:
            Number of unpaired electrons (0 if none)
        """
        try:
            # Check if spin multiplicity is already set
            spin_mult = atom.getProperty(self.SPIN_MULTIPLICITY)
            if spin_mult is not None:
                # Spin multiplicity = 2S + 1, where S is total spin
                # For radicals: doublet (S=1/2, mult=2), triplet (S=1, mult=3), etc.
                return max(0, spin_mult - 1)

            # Calculate from valence
            element = atom.getSymbol()
            formal_charge = (
                atom.getFormalCharge() if atom.getFormalCharge() is not None else 0
            )

            # Get number of bonds
            bond_count = 0
            bond_order_sum = 0.0
            for bond in molecule.getConnectedBondsList(atom):
                bond_count += 1
                order = bond.getOrder()
                if order is not None:
                    # Convert CDK IBond.Order to numeric value
                    if order.numeric() is not None:
                        bond_order_sum += order.numeric()
                    else:
                        bond_order_sum += 1.0  # Assume single bond if unknown

            # Get implicit hydrogen count
            implicit_h = atom.getImplicitHydrogenCount()
            if implicit_h is None:
                implicit_h = 0

            # Calculate expected valence for common elements
            expected_valence = self._get_expected_valence(element, formal_charge)
            if expected_valence is None:
                return 0  # Unknown element

            # Calculate actual valence
            actual_valence = bond_order_sum + implicit_h

            # Radical electrons = expected - actual
            radical_count = int(expected_valence - actual_valence)

            return max(0, radical_count)

        except Exception:
            return 0

    def _get_expected_valence(self, element: str, formal_charge: int) -> Optional[int]:
        """Get expected valence for an element considering formal charge.

        Args:
            element: Element symbol (e.g., "C", "N", "O")
            formal_charge: Formal charge on atom

        Returns:
            Expected valence or None if unknown
        """
        # Common valences (neutral atoms)
        base_valences = {
            "H": 1,
            "He": 0,
            "B": 3,
            "C": 4,
            "N": 3,
            "O": 2,
            "F": 1,
            "Ne": 0,
            "Al": 3,
            "Si": 4,
            "P": 3,
            "S": 2,
            "Cl": 1,
            "Ar": 0,
            "Br": 1,
            "I": 1,
        }

        if element not in base_valences:
            return None

        base_valence = base_valences[element]

        # Adjust for formal charge
        # Positive charge generally reduces available electrons
        # Negative charge generally adds electrons
        adjusted_valence = base_valence - formal_charge

        return max(0, adjusted_valence)

    def mark_radicals(self, molecule: Any) -> None:
        """Mark radical atoms with spin multiplicity property.

        Sets the SPIN_MULTIPLICITY property on atoms that have unpaired electrons.

        Args:
            molecule: CDK IAtomContainer

        Example:
            >>> perceiver = RadicalPerception()
            >>> perceiver.mark_radicals(mol)
        """
        try:
            radicals = self.perceive_radicals(molecule)

            for atom_idx, radical_count in radicals.items():
                atom = molecule.getAtom(atom_idx)
                # Spin multiplicity = 2S + 1
                spin_mult = radical_count + 1
                atom.setProperty(self.SPIN_MULTIPLICITY, self.Integer(spin_mult))

        except Exception:
            pass

    def count_radicals(self, molecule: Any) -> int:
        """Count total number of radical electrons in molecule.

        Args:
            molecule: CDK IAtomContainer

        Returns:
            Total number of unpaired electrons
        """
        try:
            radicals = self.perceive_radicals(molecule)
            total = sum(radicals.values())
            return total
        except Exception:
            return 0

    def get_radical_atoms(self, molecule: Any) -> List[int]:
        """Get list of atom indices that have radicals.

        Args:
            molecule: CDK IAtomContainer

        Returns:
            List of atom indices with unpaired electrons
        """
        try:
            radicals = self.perceive_radicals(molecule)
            return list(radicals.keys())
        except Exception:
            return []


def perceive_radicals(molecule: Any) -> Dict[int, int]:
    """Convenience function to detect radicals in a molecule.

    Args:
        molecule: CDK IAtomContainer

    Returns:
        Dictionary mapping atom index to number of unpaired electrons

    Example:
        >>> radicals = perceive_radicals(mol)
        >>> print(f"Found {len(radicals)} atoms with radicals")
    """
    try:
        perceiver = RadicalPerception()
        return perceiver.perceive_radicals(molecule)
    except Exception:
        return {}


def mark_radicals(molecule: Any) -> None:
    """Convenience function to mark radical atoms.

    Sets spin multiplicity properties on atoms with unpaired electrons.

    Args:
        molecule: CDK IAtomContainer

    Example:
        >>> mark_radicals(mol)
    """
    try:
        perceiver = RadicalPerception()
        perceiver.mark_radicals(molecule)
    except Exception:
        pass


# Export public API
__all__ = [
    "RadicalPerception",
    "perceive_radicals",
    "mark_radicals",
]
