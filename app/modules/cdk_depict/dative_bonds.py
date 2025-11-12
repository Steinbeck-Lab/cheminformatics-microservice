"""Dative Bond Perception Module.

This module provides comprehensive support for perceiving and displaying
dative (coordinate) bonds in chemical structures. Dative bonds occur when
both electrons in a bond come from one atom (donor) to another (acceptor).

Common cases:
- N, P, O → Metal (coordinate complexes)
- N, P → B (borane complexes)
- O → S, P (sulfoxides, phosphine oxides)

Based on CDK's dative bond perception with three modes:
- Always: Perceive all dative bonds
- Metals: Only metal-ligand dative bonds (default)
- Never: No dative bond perception
"""

from __future__ import annotations

from enum import Enum

from jpype import JClass


class DativeBondMode(Enum):
    """Modes for dative bond perception."""

    ALWAYS = "always"  # Perceive all dative bonds
    METALS = "metals"  # Only metal-ligand bonds (default)
    NEVER = "never"  # No dative bond perception


class DativeBondPerception:
    """Perception and display of dative (coordinate) bonds.

    Dative bonds are coordinate covalent bonds where both electrons come
    from one atom. This class detects and marks such bonds for proper
    depiction with arrow notation.

    Example:
        >>> perceiver = DativeBondPerception()
        >>> mol = get_CDK_IAtomContainer("[NH3]B(F)(F)F")
        >>> perceiver.perceive(mol, mode=DativeBondMode.ALWAYS)
        >>> # N→B bond will be marked as dative
    """

    def __init__(self, cdk_base: str = "org.openscience.cdk"):
        """Initialize dative bond perception.

        Args:
            cdk_base: Base package path for CDK classes
        """
        self.cdk_base = cdk_base
        # Cache Java Integer class for charge setting
        self._Integer = JClass("java.lang.Integer")

    def perceive(
        self, molecule: any, mode: DativeBondMode = DativeBondMode.METALS
    ) -> int:
        """Perceive and mark dative bonds in a molecule.

        Args:
            molecule: CDK IAtomContainer
            mode: Dative bond perception mode

        Returns:
            Number of dative bonds detected and marked
        """
        if mode == DativeBondMode.NEVER:
            return 0

        IBond = JClass(self.cdk_base + ".interfaces.IBond")

        count = 0

        try:
            # First pass: Handle charged dative bonds (positive donor, negative acceptor)
            for bond in molecule.bonds():
                begin = bond.getBegin()
                end = bond.getEnd()

                # Check for positive donor → negative acceptor
                if self._is_pos_dative_donor(
                    end, mode
                ) and self._is_neg_dative_acceptor(begin, mode):
                    bond.setDisplay(IBond.Display.ArrowBeg)
                    # Adjust formal charges
                    begin_charge = (
                        begin.getFormalCharge() if begin.getFormalCharge() else 0
                    )
                    end_charge = end.getFormalCharge() if end.getFormalCharge() else 0
                    # FIXED: Convert to Java Integer objects
                    begin.setFormalCharge(self._Integer(begin_charge + 1))
                    end.setFormalCharge(self._Integer(end_charge - 1))
                    count += 1

                elif self._is_pos_dative_donor(
                    begin, mode
                ) and self._is_neg_dative_acceptor(end, mode):
                    bond.setDisplay(IBond.Display.ArrowEnd)
                    # Adjust formal charges
                    begin_charge = (
                        begin.getFormalCharge() if begin.getFormalCharge() else 0
                    )
                    end_charge = end.getFormalCharge() if end.getFormalCharge() else 0
                    # FIXED: Convert to Java Integer objects
                    begin.setFormalCharge(self._Integer(begin_charge - 1))
                    end.setFormalCharge(self._Integer(end_charge + 1))
                    count += 1

            # Second pass: Handle neutral dative bonds
            for bond in molecule.bonds():
                begin = bond.getBegin()
                end = bond.getEnd()

                # Check for neutral donor → acceptor
                if self._is_dative_donor(end, mode) and self._is_dative_acceptor(
                    begin, mode
                ):
                    if (
                        bond.getDisplay() != IBond.Display.ArrowBeg
                    ):  # Not already marked
                        bond.setDisplay(IBond.Display.ArrowBeg)
                        count += 1

                elif self._is_dative_donor(begin, mode) and self._is_dative_acceptor(
                    end, mode
                ):
                    if (
                        bond.getDisplay() != IBond.Display.ArrowEnd
                    ):  # Not already marked
                        bond.setDisplay(IBond.Display.ArrowEnd)
                        count += 1

        except Exception as e:
            pass

        return count

    def _calc_valence(self, atom: any) -> int:
        """Calculate valence of an atom.

        Args:
            atom: CDK IAtom

        Returns:
            Total valence (implicit H + bond orders)
        """
        IBond = JClass(self.cdk_base + ".interfaces.IBond")

        h_count = atom.getImplicitHydrogenCount()
        valence = h_count if h_count is not None else 0

        for bond in atom.bonds():
            order = bond.getOrder()
            if order is not None and order != IBond.Order.UNSET:
                valence += order.numeric()

        return valence

    def _is_dative_donor(self, atom: any, mode: DativeBondMode) -> bool:
        """Check if atom can be a dative bond donor (electron pair donor).

        Typical donors:
        - N with 4 bonds, formal charge 0 (e.g., NH3→BF3)
        - P with 4 bonds, formal charge 0
        - O with 3 bonds, formal charge 0 (e.g., R3N→O)

        Args:
            atom: CDK IAtom
            mode: Dative bond mode

        Returns:
            True if atom can donate electron pair
        """
        atomic_num = atom.getAtomicNumber()
        charge = atom.getFormalCharge()
        if charge is None:
            charge = 0

        if charge != 0:
            return False

        valence = self._calc_valence(atom)

        # Nitrogen (7) or Phosphorus (15)
        if atomic_num in [7, 15]:
            return valence == 4

        # Oxygen (8)
        if atomic_num == 8:
            return valence == 3

        return False

    def _is_dative_acceptor(self, atom: any, mode: DativeBondMode) -> bool:
        """Check if atom can be a dative bond acceptor (electron pair acceptor).

        Typical acceptors:
        - Metals (any mode)
        - B with 4 bonds (only in ALWAYS mode)
        - O with 1 bond (only in ALWAYS mode)

        Args:
            atom: CDK IAtom
            mode: Dative bond mode

        Returns:
            True if atom can accept electron pair
        """
        Elements = JClass(self.cdk_base + ".config.Elements")

        # Metals are always acceptors
        if Elements.isMetal(atom):
            return True

        # Non-metals only in ALWAYS mode
        if mode == DativeBondMode.METALS:
            return False

        atomic_num = atom.getAtomicNumber()
        charge = atom.getFormalCharge()
        if charge is None:
            charge = 0

        if charge != 0:
            return False

        valence = self._calc_valence(atom)

        # Boron (5) - electron deficient
        if atomic_num == 5:
            return valence == 4

        # Oxygen (8) - can accept in some cases
        if atomic_num == 8:
            return valence == 1

        return False

    def _is_pos_dative_donor(self, atom: any, mode: DativeBondMode) -> bool:
        """Check if positively charged atom can be a dative donor.

        Examples:
        - [NH4]+ can donate to form [NH3]→X
        - [OH3]+ can donate to form H2O→X

        Args:
            atom: CDK IAtom
            mode: Dative bond mode

        Returns:
            True if positive atom can donate
        """
        atomic_num = atom.getAtomicNumber()
        charge = atom.getFormalCharge()
        if charge is None:
            charge = 0

        if charge != 1:
            return False

        valence = self._calc_valence(atom)

        # Nitrogen (7) or Phosphorus (15)
        if atomic_num in [7, 15]:
            return valence == 4

        # Oxygen (8)
        if atomic_num == 8:
            return valence == 3

        return False

    def _is_neg_dative_acceptor(self, atom: any, mode: DativeBondMode) -> bool:
        """Check if negatively charged atom can be a dative acceptor.

        Args:
            atom: CDK IAtom
            mode: Dative bond mode

        Returns:
            True if negative atom can accept
        """
        Elements = JClass(self.cdk_base + ".config.Elements")

        charge = atom.getFormalCharge()
        if charge is None:
            charge = 0

        if charge != -1:
            return False

        # Metals can accept even when negatively charged
        if Elements.isMetal(atom):
            return True

        # Non-metals only in ALWAYS mode
        if mode == DativeBondMode.METALS:
            return False

        atomic_num = atom.getAtomicNumber()
        valence = self._calc_valence(atom)

        # Boron (5)
        if atomic_num == 5:
            return valence == 4

        # Oxygen (8)
        if atomic_num == 8:
            return valence == 1

        return False


def perceive_dative_bonds(
    molecule: any, mode: str = "metals", cdk_base: str = "org.openscience.cdk"
) -> int:
    """Convenience function to perceive dative bonds.

    Args:
        molecule: CDK IAtomContainer
        mode: Mode string ("always", "metals", "never")
        cdk_base: Base CDK package path

    Returns:
        Number of dative bonds perceived
    """
    mode_enum = get_dative_bond_mode(mode)
    perceiver = DativeBondPerception(cdk_base)
    return perceiver.perceive(molecule, mode_enum)


def get_dative_bond_mode(mode_str: str) -> DativeBondMode:
    """Convert string to DativeBondMode enum.

    Args:
        mode_str: Mode string

    Returns:
        DativeBondMode enum value
    """
    mode_str = mode_str.lower().strip()

    if mode_str in ["y", "yes", "true", "always", "all"]:
        return DativeBondMode.ALWAYS
    elif mode_str in ["m", "metals", "metal"]:
        return DativeBondMode.METALS
    elif mode_str in ["n", "no", "false", "never", "none"]:
        return DativeBondMode.NEVER
    else:
        return DativeBondMode.METALS
