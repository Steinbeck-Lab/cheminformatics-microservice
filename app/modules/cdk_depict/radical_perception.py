"""Radical Perception for CDK Molecular Depictions.

This module provides functionality to detect and display radical electrons
(unpaired electrons) on atoms in molecular structures.

Direct port of CDK Java MolOp.perceiveRadicals() implementation.
"""

from __future__ import annotations

from typing import Any

from jpype import JClass


class RadicalPerception:
    """Manages radical detection and display for molecular depictions.

    This is a direct port of the CDK Java MolOp.perceiveRadicals() method
    which uses ISingleElectron objects for radical representation.

    Attributes:
        cdk_base: Base package path for CDK classes
    """

    def __init__(self):
        """Initialize the radical perception system."""
        self.cdk_base = "org.openscience.cdk"

    def _calc_valence(self, atom: Any, molecule: Any) -> int:
        """Calculate valence from implicit hydrogens and bond orders.

        This replicates the Java calcValence() method exactly.

        Args:
            atom: CDK IAtom
            molecule: CDK IAtomContainer

        Returns:
            Total valence (implicit H + sum of bond orders)
        """
        v = atom.getImplicitHydrogenCount()
        if v is None:
            v = 0

        for bond in molecule.getConnectedBondsList(atom):
            order = bond.getOrder()
            if order is not None:
                numeric = order.numeric()
                if numeric is not None:
                    v += int(numeric)

        return v

    def perceive_radicals(self, molecule: Any) -> None:
        """Detect and mark radical electrons using CDK's ISingleElectron.

        This is a direct port of the Java MolOp.perceiveRadicals() method.
        Radicals are represented by adding ISingleElectron objects to the molecule.

        The logic handles:
        - Carbon (atomic number 6): Can have 1-2 radicals depending on valence
        - Nitrogen (atomic number 7): Can have 1 radical when under-coordinated
        - Oxygen (atomic number 8): Can have 1-2 radicals depending on valence

        Only processes neutral (formal charge = 0) and non-aromatic atoms.

        Args:
            molecule: CDK IAtomContainer to process

        Example:
            >>> perceiver = RadicalPerception()
            >>> perceiver.perceive_radicals(mol)
        """
        try:
            for atom in molecule.atoms():
                q = atom.getFormalCharge()
                if q is None:
                    q = 0

                # Skip aromatic atoms
                if atom.isAromatic():
                    continue

                atomic_num = atom.getAtomicNumber()

                # Carbon (atomic number 6)
                if atomic_num == 6 and q == 0:
                    v = self._calc_valence(atom, molecule)
                    # Note: These are separate if statements, both can execute
                    if v == 2:
                        molecule.addSingleElectron(molecule.indexOf(atom))
                    if v < 4:
                        molecule.addSingleElectron(molecule.indexOf(atom))

                # Nitrogen (atomic number 7)
                elif atomic_num == 7 and q == 0:
                    v = self._calc_valence(atom, molecule)
                    if v < 3:
                        molecule.addSingleElectron(molecule.indexOf(atom))

                # Oxygen (atomic number 8)
                elif atomic_num == 8 and q == 0:
                    v = self._calc_valence(atom, molecule)
                    # Note: These are separate if statements, both can execute
                    if v < 2:
                        molecule.addSingleElectron(molecule.indexOf(atom))
                    if v < 1:
                        molecule.addSingleElectron(molecule.indexOf(atom))

        except Exception:
            # Silent fail like Java implementation
            pass

    def perceive_radicals_reaction(self, reaction: Any) -> None:
        """Perceive radicals in all components of a reaction.

        Applies radical perception to:
        - All reactant molecules
        - All product molecules
        - All agent/catalyst molecules

        Args:
            reaction: CDK IReaction

        Example:
            >>> perceiver = RadicalPerception()
            >>> perceiver.perceive_radicals_reaction(reaction)
        """
        try:
            # Process reactants
            reactants = reaction.getReactants()
            for mol in reactants.atomContainers():
                self.perceive_radicals(mol)

            # Process products
            products = reaction.getProducts()
            for mol in products.atomContainers():
                self.perceive_radicals(mol)

            # Process agents/catalysts
            agents = reaction.getAgents()
            for mol in agents.atomContainers():
                self.perceive_radicals(mol)

        except Exception:
            pass

    def perceive_radicals_reaction_set(self, reaction_set: Any) -> None:
        """Perceive radicals in all reactions in a reaction set.

        Args:
            reaction_set: CDK IReactionSet

        Example:
            >>> perceiver = RadicalPerception()
            >>> perceiver.perceive_radicals_reaction_set(rxn_set)
        """
        try:
            for rxn in reaction_set.reactions():
                self.perceive_radicals_reaction(rxn)
        except Exception:
            pass


def perceive_radicals(molecule_or_reaction: Any) -> None:
    """Convenience function to detect and mark radicals.

    Automatically detects whether input is a molecule, reaction, or reaction set
    and applies appropriate radical perception.

    Args:
        molecule_or_reaction: CDK IAtomContainer, IReaction, or IReactionSet

    Example:
        >>> from app.modules.cdk_depict.radical_perception import perceive_radicals
        >>> perceive_radicals(mol)
        >>> perceive_radicals(reaction)
        >>> perceive_radicals(reaction_set)
    """
    try:
        perceiver = RadicalPerception()
        cdk_base = "org.openscience.cdk"

        # Check if it's a reaction set
        try:
            IReactionSet = JClass(cdk_base + ".interfaces.IReactionSet")
            if isinstance(molecule_or_reaction, IReactionSet):
                perceiver.perceive_radicals_reaction_set(molecule_or_reaction)
                return
        except Exception:
            pass

        # Check if it's a reaction
        try:
            IReaction = JClass(cdk_base + ".interfaces.IReaction")
            if isinstance(molecule_or_reaction, IReaction):
                perceiver.perceive_radicals_reaction(molecule_or_reaction)
                return
        except Exception:
            pass

        # Otherwise treat as molecule
        perceiver.perceive_radicals(molecule_or_reaction)

    except Exception:
        pass


# Export public API
__all__ = [
    "RadicalPerception",
    "perceive_radicals",
]
