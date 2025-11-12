"""Hydrogen Display Control for CDK Molecular Depictions.

This module provides comprehensive control over hydrogen atom display modes
in molecular structures, including:
- Minimal: Suppress all hydrogens (implicit only)
- Explicit: Show all hydrogens as explicit atoms
- Stereo: Show only stereo-relevant hydrogens (chiral centers, E/Z bonds)
- Smart: Intelligent hydrogen display considering ring systems
- Provided: Keep hydrogens as provided (no changes)

Based on CDK hydrogen display functionality.

"""

from __future__ import annotations

import logging
from typing import Literal, Any

from jpype import JClass

logger = logging.getLogger(__name__)


HydrogenDisplayType = Literal["Minimal", "Explicit", "Stereo", "Smart", "Provided"]


def set_hydrogen_display(
    molecule: any,
    display_type: HydrogenDisplayType = "Provided",
    cdk_base: str = "org.openscience.cdk",
) -> None:
    """Set the hydrogen display mode for a CDK IAtomContainer molecule.

    This function controls how hydrogen atoms are displayed in molecular structures
    by modifying the molecule in-place. It supports five different display modes:
    - Minimal: Suppress all hydrogens (implicit only)
    - Explicit: Show all hydrogens as explicit atoms
    - Stereo: Show only stereo-relevant hydrogens (chiral centers, E/Z bonds)
    - Smart: Intelligent hydrogen display considering ring systems
    - Provided: Keep hydrogens as provided (no changes)

    The function directly modifies the input molecule and does not return a value.

    Args:
        molecule (any): CDK IAtomContainer molecule object to modify. This should be
                       a parsed CDK molecule obtained from functions like
                       get_CDK_IAtomContainer(smiles).
        display_type (str, optional): The hydrogen display mode to apply. Must be one of:
                                     "Minimal", "Explicit", "Stereo", "Smart", or "Provided".
                                     Defaults to "Provided" (no changes).
        cdk_base (str, optional): Base package path for CDK classes.
                                 Defaults to "org.openscience.cdk". Only change this if
                                 using a custom CDK installation.

    Returns:
        None: The molecule is modified in-place.

    Raises:
        ValueError: If display_type is not one of the valid options.

    Example:
        >>> from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
        >>>
        >>> # Create a chiral molecule
        >>> molecule = get_CDK_IAtomContainer("C[C@H](N)C(=O)O")  # (S)-Alanine
        >>>
        >>> # Show only stereo-relevant hydrogens
        >>> set_hydrogen_display(molecule, "Stereo")
        >>>
        >>> # Now the hydrogen at the chiral center is explicit
        >>> print(f"Atoms: {molecule.getAtomCount()}")

    Notes:
        - Apply hydrogen display BEFORE generating 2D coordinates for best results
        - For chiral molecules, use "Stereo" or "Smart" mode
        - For publication-quality images, use "Minimal" mode
        - The molecule object is modified in-place and cannot be undone
    """

    AtomContainerManipulator = JClass(
        cdk_base + ".tools.manipulator.AtomContainerManipulator"
    )

    if display_type == "Minimal":
        AtomContainerManipulator.suppressHydrogens(molecule)

    elif display_type == "Explicit":
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule)

    elif display_type == "Stereo":
        _set_hydrogen_display_stereo(molecule, cdk_base)

    elif display_type == "Smart":
        _set_hydrogen_display_smart(molecule, cdk_base)

    elif display_type == "Provided":
        pass

    else:
        raise ValueError(
            f"Invalid display_type: {display_type}. "
            f"Must be one of: 'Minimal', 'Explicit', 'Stereo', 'Smart', 'Provided'"
        )


def _set_hydrogen_display_stereo(molecule: any, cdk_base: str) -> None:
    """Internal function to handle 'Stereo' hydrogen display mode.

    This function suppresses all hydrogens and then adds them back explicitly
    only for stereogenic centers (chiral atoms, E/Z double bonds, and allenes).
    The stereo elements are rebuilt with explicit hydrogens where needed.

    Args:
        molecule (any): CDK IAtomContainer molecule to modify.
        cdk_base (str): Base package path for CDK classes.

    Returns:
        None: The molecule is modified in-place.

    Notes:
        - Handles tetrahedral chirality (R/S configuration)
        - Handles cis-trans double bonds (E/Z configuration)
        - Handles allenes (extended tetrahedral chirality)
        - Preserves stereo element group information (racemic, relative)
    """
    AtomContainerManipulator = JClass(
        cdk_base + ".tools.manipulator.AtomContainerManipulator"
    )
    TetrahedralChirality = JClass(cdk_base + ".stereo.TetrahedralChirality")
    IStereoElement = JClass(cdk_base + ".interfaces.IStereoElement")
    ExtendedTetrahedral = JClass(cdk_base + ".stereo.ExtendedTetrahedral")
    ArrayList = JClass("java.util.ArrayList")
    HashMap = JClass("java.util.HashMap")
    Integer = JClass("java.lang.Integer")

    AtomContainerManipulator.suppressHydrogens(molecule)

    new_stereo_elements = ArrayList()

    for stereo_element in molecule.stereoElements():
        config_class = stereo_element.getConfigClass()

        if config_class == IStereoElement.Tetrahedral:
            focus = stereo_element.getFocus()

            if focus.getImplicitHydrogenCount() == 1:
                focus.setImplicitHydrogenCount(Integer(0))
                hydrogen = _sprout_hydrogen(molecule, focus)

                # Create Java HashMap for mapping
                mapping = HashMap()
                mapping.put(focus, hydrogen)
                temp_stereo = stereo_element.map(mapping)

                # Get carriers as Java array
                carriers_list = []
                for carrier in temp_stereo.getCarriers():
                    carriers_list.append(carrier)

                new_stereo = TetrahedralChirality(
                    focus, carriers_list, temp_stereo.getConfig()
                )
                new_stereo.setGroupInfo(stereo_element.getGroupInfo())
                new_stereo_elements.add(new_stereo)
            else:
                new_stereo_elements.add(stereo_element)

        elif config_class == IStereoElement.CisTrans:
            focus_bond = stereo_element.getFocus()
            begin_atom = focus_bond.getBegin()
            end_atom = focus_bond.getEnd()

            if begin_atom.getImplicitHydrogenCount() == 1:
                begin_atom.setImplicitHydrogenCount(Integer(0))
                _sprout_hydrogen(molecule, begin_atom)

            if end_atom.getImplicitHydrogenCount() == 1:
                end_atom.setImplicitHydrogenCount(Integer(0))
                _sprout_hydrogen(molecule, end_atom)

            new_stereo_elements.add(stereo_element)

        elif config_class == IStereoElement.Allenal:
            focus = stereo_element.getFocus()
            terminals = ExtendedTetrahedral.findTerminalAtoms(molecule, focus)

            mapping = HashMap()
            if terminals[0].getImplicitHydrogenCount() == 1:
                terminals[0].setImplicitHydrogenCount(Integer(0))
                hydrogen = _sprout_hydrogen(molecule, terminals[0])
                mapping.put(terminals[0], hydrogen)

            if terminals[1].getImplicitHydrogenCount() == 1:
                terminals[1].setImplicitHydrogenCount(Integer(0))
                hydrogen = _sprout_hydrogen(molecule, terminals[1])
                mapping.put(terminals[1], hydrogen)

            if not mapping.isEmpty():
                new_stereo_elements.add(stereo_element.map(mapping))
            else:
                new_stereo_elements.add(stereo_element)

        else:
            new_stereo_elements.add(stereo_element)

    molecule.setStereoElements(new_stereo_elements)


def _set_hydrogen_display_smart(molecule: any, cdk_base: str) -> None:
    """Internal function to handle 'Smart' hydrogen display mode.

    This function is similar to Stereo mode but uses intelligent heuristics to decide
    when hydrogens should be displayed. It considers ring systems, neighboring stereo
    centers, and coordination geometries to make smarter decisions about hydrogen
    placement. This mode is particularly useful for complex molecules with multiple
    stereocenters or unusual coordination geometries.

    Args:
        molecule (any): CDK IAtomContainer molecule to modify.
        cdk_base (str): Base package path for CDK classes.

    Returns:
        None: The molecule is modified in-place.

    Notes:
        - Uses _should_add_hydrogen() to make intelligent decisions
        - Handles standard tetrahedral chirality with smart logic
        - Supports advanced coordination geometries:
          * Square planar (e.g., Pt(II) complexes)
          * Trigonal bipyramidal (e.g., Fe(CO)5)
          * Octahedral (e.g., Fe(CN)6)
        - Marks ring atoms and bonds for context-aware decisions
        - Considers hydrogen isotopes (D, T) which should always be shown
    """
    AtomContainerManipulator = JClass(
        cdk_base + ".tools.manipulator.AtomContainerManipulator"
    )
    Cycles = JClass(cdk_base + ".graph.Cycles")
    TetrahedralChirality = JClass(cdk_base + ".stereo.TetrahedralChirality")
    SquarePlanar = JClass(cdk_base + ".stereo.SquarePlanar")
    TrigonalBipyramidal = JClass(cdk_base + ".stereo.TrigonalBipyramidal")
    Octahedral = JClass(cdk_base + ".stereo.Octahedral")
    IStereoElement = JClass(cdk_base + ".interfaces.IStereoElement")
    ExtendedTetrahedral = JClass(cdk_base + ".stereo.ExtendedTetrahedral")
    ArrayList = JClass("java.util.ArrayList")
    HashMap = JClass("java.util.HashMap")
    Integer = JClass("java.lang.Integer")

    AtomContainerManipulator.suppressHydrogens(molecule)
    Cycles.markRingAtomsAndBonds(molecule)

    new_stereo_elements = ArrayList()

    for stereo_element in molecule.stereoElements():
        config_class = stereo_element.getConfigClass()

        if config_class == IStereoElement.Tetrahedral:
            focus = stereo_element.getFocus()

            if focus.getImplicitHydrogenCount() == 1:
                bonds = molecule.getConnectedBondsList(focus)

                if _should_add_hydrogen(molecule, focus, bonds):
                    carriers = _get_explicit_hydrogen_carriers(
                        molecule, stereo_element, focus
                    )
                    new_stereo = TetrahedralChirality(
                        focus, carriers, stereo_element.getConfig()
                    )
                    new_stereo.setGroupInfo(stereo_element.getGroupInfo())
                    new_stereo_elements.add(new_stereo)
                else:
                    new_stereo_elements.add(stereo_element)
            else:
                new_stereo_elements.add(stereo_element)

        elif config_class == IStereoElement.SquarePlanar:
            focus = stereo_element.getFocus()
            if focus.getImplicitHydrogenCount() > 0:
                carriers = _get_explicit_hydrogen_carriers(
                    molecule, stereo_element, focus
                )
                new_stereo_elements.add(
                    SquarePlanar(focus, carriers, stereo_element.getConfig())
                )
            else:
                new_stereo_elements.add(stereo_element)

        elif config_class == IStereoElement.TrigonalBipyramidal:
            focus = stereo_element.getFocus()
            if focus.getImplicitHydrogenCount() > 0:
                carriers = _get_explicit_hydrogen_carriers(
                    molecule, stereo_element, focus
                )
                new_stereo_elements.add(
                    TrigonalBipyramidal(focus, carriers, stereo_element.getConfig())
                )
            else:
                new_stereo_elements.add(stereo_element)

        elif config_class == IStereoElement.Octahedral:
            focus = stereo_element.getFocus()
            if focus.getImplicitHydrogenCount() > 0:
                carriers = _get_explicit_hydrogen_carriers(
                    molecule, stereo_element, focus
                )
                new_stereo_elements.add(
                    Octahedral(focus, carriers, stereo_element.getConfig())
                )
            else:
                new_stereo_elements.add(stereo_element)

        elif config_class == IStereoElement.CisTrans:
            focus_bond = stereo_element.getFocus()
            begin_atom = focus_bond.getBegin()
            end_atom = focus_bond.getEnd()

            mapping = HashMap()

            if begin_atom.getImplicitHydrogenCount() == 1:
                bonds = molecule.getConnectedBondsList(begin_atom)
                if _should_add_hydrogen(molecule, begin_atom, bonds):
                    begin_atom.setImplicitHydrogenCount(Integer(0))
                    hydrogen = _sprout_hydrogen(molecule, begin_atom)
                    mapping.put(begin_atom, hydrogen)

            if end_atom.getImplicitHydrogenCount() == 1:
                bonds = molecule.getConnectedBondsList(end_atom)
                if _should_add_hydrogen(molecule, end_atom, bonds):
                    end_atom.setImplicitHydrogenCount(Integer(0))
                    hydrogen = _sprout_hydrogen(molecule, end_atom)
                    mapping.put(end_atom, hydrogen)

            if not mapping.isEmpty():
                new_stereo_elements.add(stereo_element.map(mapping))
            else:
                new_stereo_elements.add(stereo_element)

        elif config_class == IStereoElement.Allenal:
            focus = stereo_element.getFocus()
            terminals = ExtendedTetrahedral.findTerminalAtoms(molecule, focus)

            mapping = HashMap()
            if terminals[0].getImplicitHydrogenCount() == 1:
                terminals[0].setImplicitHydrogenCount(Integer(0))
                hydrogen = _sprout_hydrogen(molecule, terminals[0])
                mapping.put(terminals[0], hydrogen)

            if terminals[1].getImplicitHydrogenCount() == 1:
                terminals[1].setImplicitHydrogenCount(Integer(0))
                hydrogen = _sprout_hydrogen(molecule, terminals[1])
                mapping.put(terminals[1], hydrogen)

            if not mapping.isEmpty():
                new_stereo_elements.add(stereo_element.map(mapping))
            else:
                new_stereo_elements.add(stereo_element)

        else:
            new_stereo_elements.add(stereo_element)

    molecule.setStereoElements(new_stereo_elements)


def _sprout_hydrogen(molecule: any, focus_atom: any) -> any:
    """Add an explicit hydrogen atom connected to the focus atom.

    This utility function creates a new hydrogen atom, adds it to the molecule,
    and creates a single bond between the hydrogen and the specified focus atom.
    The hydrogen is given atomic number 1, symbol "H", and zero implicit hydrogens.

    Args:
        molecule (any): CDK IAtomContainer molecule to add the hydrogen to.
        focus_atom (any): The CDK IAtom to which the hydrogen will be bonded.

    Returns:
        any: The newly created hydrogen atom (CDK IAtom object).

    Notes:
        - The hydrogen is added to the end of the atom list
        - A single bond is created automatically
        - The bond uses the CDK IBond.Order.SINGLE constant
        - Proper Java Integer wrapper is used for setAtomicNumber
    """
    IBond = JClass("org.openscience.cdk.interfaces.IBond")
    Integer = JClass("java.lang.Integer")

    hydrogen = molecule.getBuilder().newAtom()
    hydrogen.setAtomicNumber(Integer(1))
    hydrogen.setSymbol("H")
    hydrogen.setImplicitHydrogenCount(Integer(0))

    molecule.addAtom(hydrogen)

    focus_idx = molecule.indexOf(focus_atom)
    hydrogen_idx = molecule.getAtomCount() - 1
    molecule.addBond(focus_idx, hydrogen_idx, IBond.Order.SINGLE)

    return molecule.getAtom(hydrogen_idx)


def _should_add_hydrogen(molecule: any, atom: any, bonds) -> bool:
    """Determine if a hydrogen should be added based on smart heuristics.

    This function implements intelligent decision-making for hydrogen display by
    considering the chemical environment of the atom. It counts "important" neighbors
    (ring bonds and neighboring stereo centers) to decide if showing the hydrogen
    would improve structure clarity. The heuristic is: if an atom has 3 important
    neighbors, the 4th position (with the hydrogen) is stereochemically relevant.

    Args:
        molecule (any): CDK IAtomContainer molecule being processed.
        atom (any): The CDK IAtom under consideration.
        bonds (iterable): Iterable of CDK IBond objects connected to the atom.

    Returns:
        bool: True if hydrogen should be added explicitly, False otherwise.

    Notes:
        - Counts ring bonds as "important" (contributes to stereochemistry)
        - Counts bonds to neighboring stereo centers as "important"
        - Always returns True for hydrogen isotopes (deuterium, tritium)
        - Returns True if count == 3 (meaning the H position is stereochemically relevant)

    Example:
        A chiral carbon in a ring with 3 ring bonds: the H is important
        A chiral carbon with 3 substituents: the H is stereochemically relevant
    """
    IStereoElement = JClass("org.openscience.cdk.interfaces.IStereoElement")

    count = 0

    for bond in bonds:
        neighbor = bond.getOther(atom)

        if bond.isInRing():
            count += 1
        else:
            for stereo_element in molecule.stereoElements():
                if (
                    stereo_element.getConfigClass() == IStereoElement.TH
                    and stereo_element.getFocus().equals(neighbor)
                ):
                    count += 1

        # Hydrogen isotopes (D, T) should always be shown explicitly
        if neighbor.getAtomicNumber() == 1 and neighbor.getMassNumber() is not None:
            return True

    return count == 3


def _get_explicit_hydrogen_carriers(
    molecule: any, stereo_element: any, focus_atom: any
) -> list:
    """Get carrier atoms with explicit hydrogens for stereo elements.

    This function creates explicit hydrogen atoms for all implicit hydrogens on
    the focus atom and returns a carrier list where focus atom entries in the
    stereo element's carrier list are replaced with the explicit hydrogens.
    This is necessary for properly defining stereochemistry when hydrogens need
    to be shown explicitly.

    Args:
        molecule (any): CDK IAtomContainer molecule being processed.
        stereo_element (any): The CDK IStereoElement being updated.
        focus_atom (any): The CDK IAtom that is the focus of the stereo element.

    Returns:
        list: Python list of carrier atoms with explicit hydrogens replacing focus
              atom entries. This list can be passed to CDK stereo element constructors.

    Notes:
        - Creates one hydrogen for each implicit hydrogen on focus atom
        - Sets focus atom's implicit hydrogen count to 0
        - Replaces focus atom entries in carriers with the new hydrogens
        - Other carrier atoms are kept unchanged
        - Uses proper Java Integer wrapper for setImplicitHydrogenCount

    Example:
        For a tetrahedral center [C@H](A)(B)C, if H is implicit:
        - Creates explicit H atom
        - Carrier list [C, A, B, C] becomes [H, A, B, C]
    """
    Integer = JClass("java.lang.Integer")

    hydrogens = []
    h_count = focus_atom.getImplicitHydrogenCount()
    for i in range(h_count):
        hydrogens.append(_sprout_hydrogen(molecule, focus_atom))

    focus_atom.setImplicitHydrogenCount(Integer(0))

    carriers = []
    hydrogen_idx = 0

    for carrier in stereo_element.getCarriers():
        if carrier.equals(focus_atom) and hydrogen_idx < len(hydrogens):
            carriers.append(hydrogens[hydrogen_idx])
            hydrogen_idx += 1
        else:
            carriers.append(carrier)

    return carriers


def setHydrogenDisplay(
    molecule: any, display_type: HydrogenDisplayType = "Provided"
) -> None:
    """Convenience alias for set_hydrogen_display using camelCase naming.

    This function is provided for backwards compatibility with Java naming conventions
    and to match the naming style of the original CDK Java implementation. It simply
    calls set_hydrogen_display() with the same parameters.

    Args:
        molecule (any): CDK IAtomContainer molecule to modify.
        display_type (str, optional): The hydrogen display mode. Defaults to "Provided".

    Returns:
        None: The molecule is modified in-place.

    Example:
        >>> # Both of these are equivalent:
        >>> set_hydrogen_display(molecule, "Stereo")
        >>> setHydrogenDisplay(molecule, "Stereo")
    """
    set_hydrogen_display(molecule, display_type)
