"""Multicenter Bond Module.

This module provides comprehensive support for multicenter bonding in
organometallic and coordination compounds. Multicenter bonds occur when
a single metal atom bonds to multiple atoms simultaneously, such as:

- η5-Cyclopentadienyl (Cp) complexes: Fe(η5-C5H5)2 (ferrocene)
- η6-Benzene complexes: Cr(η6-C6H6)2
- π-Allyl complexes: Pd(η3-C3H5)
- Agostic interactions
- Bridging ligands

Display styles:
- Provided: Show as-is from input
- Dative: Arrow notation from ring to metal
- Dashed: Dashed lines (preserving charges)
- DashedNeutral: Dashed lines (neutralize charges)
- Hidden: Hide bonds (preserving charges)
- HiddenNeutral: Hide bonds (neutralize charges)

"""

from __future__ import annotations

from enum import Enum
from typing import Set

from jpype import JClass


class MulticenterStyle(Enum):
    """Display styles for multicenter bonds."""

    PROVIDED = "provided"  # As-is from input
    DATIVE = "dative"  # Arrow notation
    DASHED = "dashed"  # Dashed lines with charges
    DASHED_NEUTRAL = "dashed_neutral"  # Dashed lines, neutralize charges
    HIDDEN = "hidden"  # Hide bonds with charges
    HIDDEN_NEUTRAL = "hidden_neutral"  # Hide bonds, neutralize charges


class MulticenterBonds:
    """Handler for multicenter bonding in organometallic compounds.

    Multicenter bonds (η-bonds, hapto bonds) are common in organometallic
    chemistry where a metal center bonds to multiple atoms of a ligand
    simultaneously.

    Example:
        >>> handler = MulticenterBonds()
        >>> mol = get_CDK_IAtomContainer("[Fe]c1ccccc1")  # Ferrocene fragment
        >>> handler.set_style(mol, MulticenterStyle.DATIVE)
    """

    def __init__(self, cdk_base: str = "org.openscience.cdk"):
        """Initialize multicenter bond handler.

        Args:
            cdk_base: Base package path for CDK classes
        """
        self.cdk_base = cdk_base

    def set_style(
        self, molecule: any, style: MulticenterStyle = MulticenterStyle.PROVIDED
    ) -> int:
        """Set display style for multicenter bonds in a molecule.

        Args:
            molecule: CDK IAtomContainer
            style: Display style to apply

        Returns:
            Number of multicenter bonds processed
        """
        if style == MulticenterStyle.PROVIDED:
            return 0

        try:
            CDKConstants = JClass(self.cdk_base + ".CDKConstants")
            SgroupType = JClass(self.cdk_base + ".sgroup.SgroupType")
            Elements = JClass(self.cdk_base + ".config.Elements")

            # Get Sgroups from molecule
            sgroups = molecule.getProperty(CDKConstants.CTAB_SGROUPS)
            if sgroups is None:
                return 0

            count = 0

            # Process each Sgroup
            for sgroup in sgroups:
                # Only process ExtMulticenter Sgroups
                if sgroup.getType() != SgroupType.ExtMulticenter:
                    continue

                atoms = sgroup.getAtoms()
                bonds = sgroup.getBonds()

                # Must have exactly one bond
                if bonds.size() != 1:
                    continue

                bond = bonds.iterator().next()
                begin = bond.getBegin()
                end = bond.getEnd()

                # Determine which end is the metal and which is the attachment point
                if Elements.isMetal(begin) and atoms.contains(end):
                    # Metal is at begin
                    self._apply_bond_style(bond, style)

                    if style in [
                        MulticenterStyle.DATIVE,
                        MulticenterStyle.DASHED_NEUTRAL,
                        MulticenterStyle.HIDDEN_NEUTRAL,
                    ]:
                        self._neutralize_charges(atoms, begin, end)

                    count += 1

                elif Elements.isMetal(end) and atoms.contains(begin):
                    # Metal is at end
                    self._apply_bond_style(bond, style)

                    if style in [
                        MulticenterStyle.DATIVE,
                        MulticenterStyle.DASHED_NEUTRAL,
                        MulticenterStyle.HIDDEN_NEUTRAL,
                    ]:
                        self._neutralize_charges(atoms, end, begin)

                    count += 1

        except Exception as e:
            pass

        return count

    def _apply_bond_style(self, bond: any, style: MulticenterStyle) -> None:
        """Apply display style to a bond.

        Args:
            bond: CDK IBond
            style: Display style to apply
        """
        IBond = JClass(self.cdk_base + ".interfaces.IBond")
        StandardGenerator = JClass(
            self.cdk_base + ".renderer.generators.standard.StandardGenerator"
        )

        if style in [MulticenterStyle.HIDDEN, MulticenterStyle.HIDDEN_NEUTRAL]:
            # Hide the bond
            bond.setProperty(StandardGenerator.HIDDEN, True)

        elif style in [MulticenterStyle.DASHED, MulticenterStyle.DASHED_NEUTRAL]:
            # Make bond dashed
            bond.setDisplay(IBond.Display.Dash)

        elif style == MulticenterStyle.DATIVE:
            # Use arrow notation (direction based on which atom has fewer bonds)
            begin = bond.getBegin()
            end = bond.getEnd()

            if begin.getBondCount() == 1:
                bond.setDisplay(IBond.Display.ArrowEnd)
            else:
                bond.setDisplay(IBond.Display.ArrowBeg)

    def _neutralize_charges(self, ring_atoms: any, metal: any, attachment: any) -> None:
        """Neutralize formal charges in a multicenter bond system.

        This adjusts the formal charges on the ring system and metal
        to make the depiction cleaner.

        Args:
            ring_atoms: Java set of atoms in the ring system
            metal: Metal atom
            attachment: Attachment point atom
        """
        try:
            # Calculate total charge on ring atoms
            charge_on_ring = 0

            for atom in ring_atoms:
                # Skip the attachment point
                if atom.equals(attachment):
                    continue

                charge = atom.getFormalCharge()
                if charge is not None:
                    charge_on_ring += charge

            # If metal can compensate the ring charge, neutralize
            metal_charge = metal.getFormalCharge()
            if metal_charge is None:
                metal_charge = 0

            if metal_charge >= -charge_on_ring:
                # Neutralize ring atoms
                for atom in ring_atoms:
                    atom.setFormalCharge(0)

                # Adjust metal charge
                metal.setFormalCharge(metal_charge - (-charge_on_ring))

        except Exception as e:
            pass

    def detect_multicenter_bonds(self, molecule: any) -> int:
        """Detect potential multicenter bonding patterns.

        This is a heuristic method to identify multicenter bonds that
        may not be explicitly marked in the input structure.

        Args:
            molecule: CDK IAtomContainer

        Returns:
            Number of potential multicenter bonds detected
        """
        try:
            Elements = JClass(self.cdk_base + ".config.Elements")

            count = 0

            # Look for metals bonded to aromatic systems
            for atom in molecule.atoms():
                if not Elements.isMetal(atom):
                    continue

                # Check each neighbor
                for bond in atom.bonds():
                    neighbor = bond.getOther(atom)

                    # If neighbor is aromatic and in a ring
                    if neighbor.isAromatic():
                        # This could be a multicenter bond
                        ring_atoms = self._get_aromatic_ring(molecule, neighbor)

                        if len(ring_atoms) >= 5:  # Likely η5 or η6
                            count += 1

            return count

        except Exception as e:
            return 0

    def _get_aromatic_ring(self, molecule: any, start_atom: any) -> Set[any]:
        """Get all atoms in an aromatic ring containing the start atom.

        Args:
            molecule: CDK IAtomContainer
            start_atom: Starting atom in the ring

        Returns:
            Set of atoms in the aromatic ring
        """
        visited = set()
        to_visit = [start_atom]
        ring_atoms = set()

        while to_visit:
            atom = to_visit.pop()

            if atom in visited:
                continue

            visited.add(atom)

            if atom.isAromatic():
                ring_atoms.add(atom)

                # Add aromatic neighbors
                for bond in atom.bonds():
                    if bond.isAromatic():
                        neighbor = bond.getOther(atom)
                        if neighbor not in visited:
                            to_visit.append(neighbor)

        return ring_atoms

    def mark_multicenter_bond(
        self, molecule: any, metal_index: int, ring_atom_indices: list[int]
    ) -> bool:
        """Manually mark a multicenter bond using Sgroups.

        Args:
            molecule: CDK IAtomContainer
            metal_index: Index of metal atom
            ring_atom_indices: List of ring atom indices

        Returns:
            True if successfully marked
        """
        try:
            CDKConstants = JClass(self.cdk_base + ".CDKConstants")
            Sgroup = JClass(self.cdk_base + ".sgroup.Sgroup")
            SgroupType = JClass(self.cdk_base + ".sgroup.SgroupType")
            ArrayList = JClass("java.util.ArrayList")
            HashSet = JClass("java.util.HashSet")

            # Get or create Sgroup list
            sgroups = molecule.getProperty(CDKConstants.CTAB_SGROUPS)
            if sgroups is None:
                sgroups = ArrayList()
                molecule.setProperty(CDKConstants.CTAB_SGROUPS, sgroups)

            # Create new Sgroup
            sgroup = Sgroup()
            sgroup.setType(SgroupType.ExtMulticenter)

            # Add ring atoms to Sgroup
            for idx in ring_atom_indices:
                if idx < molecule.getAtomCount():
                    sgroup.addAtom(molecule.getAtom(idx))

            # Find bond between metal and one of the ring atoms
            metal = molecule.getAtom(metal_index)
            bond_found = False

            for bond in metal.bonds():
                neighbor = bond.getOther(metal)
                neighbor_idx = molecule.indexOf(neighbor)

                if neighbor_idx in ring_atom_indices:
                    sgroup.addBond(bond)
                    bond_found = True
                    break

            if not bond_found:
                return False

            # Add Sgroup to molecule
            sgroups.add(sgroup)

            return True

        except Exception as e:
            return False


def set_multicenter_style(
    molecule: any, style: str = "provided", cdk_base: str = "org.openscience.cdk"
) -> int:
    """Convenience function to set multicenter bond style.

    Args:
        molecule: CDK IAtomContainer
        style: Style string
        cdk_base: Base CDK package path

    Returns:
        Number of multicenter bonds processed
    """
    style_enum = get_multicenter_style(style)
    handler = MulticenterBonds(cdk_base)
    return handler.set_style(molecule, style_enum)


def get_multicenter_style(style_str: str) -> MulticenterStyle:
    """Convert string to MulticenterStyle enum.

    Args:
        style_str: Style string

    Returns:
        MulticenterStyle enum value
    """
    style_str = style_str.lower().strip()

    if style_str in ["p", "provided", "default"]:
        return MulticenterStyle.PROVIDED
    elif style_str in ["d", "dative", "arrow"]:
        return MulticenterStyle.DATIVE
    elif style_str in ["a", "dashed", "dash"]:
        return MulticenterStyle.DASHED
    elif style_str in ["an", "dashed_neutral", "dashneutral"]:
        return MulticenterStyle.DASHED_NEUTRAL
    elif style_str in ["h", "hidden", "hide"]:
        return MulticenterStyle.HIDDEN
    elif style_str in ["hn", "hidden_neutral", "hideneutral"]:
        return MulticenterStyle.HIDDEN_NEUTRAL
    else:
        return MulticenterStyle.PROVIDED
