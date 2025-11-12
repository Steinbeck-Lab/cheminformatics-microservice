"""Chemical Abbreviations Module

This module provides comprehensive support for chemical abbreviations including:
- Common reagent abbreviations (THF, DMF, NaOH, etc.)
- Functional group abbreviations (Ph, Me, Et, Boc, Fmoc, etc.)
- Auto-contraction of terminal groups
- Singleton abbreviation handling
- Hydrate contraction (·nH2O notation)

Based on CDK Abbreviations system with 198 reagent and 60+ group abbreviations.
"""

from __future__ import annotations

from typing import Optional, Set
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from jpype import JClass


class AbbreviationMode(Enum):
    """Modes for abbreviation application."""

    OFF = "off"
    GROUPS = "groups"
    REAGENTS = "reagents"
    ALL = "on"  # Both groups and reagents


@dataclass
class AbbreviationOptions:
    """Configuration options for abbreviation system."""

    mode: AbbreviationMode = AbbreviationMode.REAGENTS
    allow_singleton: bool = False
    auto_contract_terminal: bool = True
    auto_contract_hetero: bool = False
    custom_group_file: Optional[str] = None
    custom_reagent_file: Optional[str] = None


class ChemicalAbbreviations:
    """Chemical abbreviation system for molecular depiction.

    This class manages two types of abbreviations:
    1. Group abbreviations (Ph, Me, Et, Boc, Fmoc, etc.)
    2. Reagent abbreviations (THF, DMF, DCM, NaOH, etc.)

    Example:
        >>> abbr = ChemicalAbbreviations()
        >>> abbr.load_default_abbreviations()
        >>> mol = get_CDK_IAtomContainer("c1ccccc1C")
        >>> abbr.apply(mol, mode=AbbreviationMode.GROUPS)
    """

    # Default CDK internal resource paths
    DEFAULT_GROUP_ABBR_PATH = "app/modules/cdk_depict/data/group_abbr.smi"
    DEFAULT_REAGENT_ABBR_PATH = "app/modules/cdk_depict/data/reagent_abbr.smi"

    def __init__(self, cdk_base: str = "org.openscience.cdk"):
        """Initialize abbreviation system.

        Args:
            cdk_base: Base package path for CDK classes
        """
        self.cdk_base = cdk_base
        self.group_abbr = None
        self.reagent_abbr = None
        self._initialized = False

    def initialize(
        self,
        custom_group_file: Optional[str] = None,
        custom_reagent_file: Optional[str] = None,
    ):
        """Initialize CDK Abbreviations objects.

        Args:
            custom_group_file: Optional path to custom group abbreviations file (local filesystem)
            custom_reagent_file: Optional path to custom reagent abbreviations file (local filesystem)
        """
        if self._initialized:
            return

        try:
            Abbreviations = JClass(self.cdk_base + ".depict.Abbreviations")

            # Initialize group abbreviations
            self.group_abbr = Abbreviations()

            # Try to load from CDK internal resources first
            try:
                self.group_abbr.loadFromFile(self.DEFAULT_GROUP_ABBR_PATH)
            except Exception as e:
                if custom_group_file and Path(custom_group_file).exists():
                    self.group_abbr.loadFromFile(custom_group_file)

            # Load additional custom group abbreviations if provided
            if custom_group_file and Path(custom_group_file).exists():
                try:
                    self.group_abbr.loadFromFile(custom_group_file)
                except Exception as e:
                    pass

            self.group_abbr.without(Abbreviations.Option.ALLOW_SINGLETON).with_(
                Abbreviations.Option.AUTO_CONTRACT_TERMINAL
            ).without(Abbreviations.Option.AUTO_CONTRACT_HETERO)

            # Initialize reagent abbreviations
            self.reagent_abbr = Abbreviations()

            # Try to load from CDK internal resources
            # NOTE: Only load reagent file here, NOT the group file
            # This ensures "reagents" mode only abbreviates standalone reagents (THF, DMF, etc.)
            # and NOT functional groups (Ph, Me, Et, etc.)
            try:
                self.reagent_abbr.loadFromFile(self.DEFAULT_REAGENT_ABBR_PATH)
            except Exception as e:
                if custom_reagent_file and Path(custom_reagent_file).exists():
                    self.reagent_abbr.loadFromFile(custom_reagent_file)

            # Load additional custom reagent abbreviations if provided
            if custom_reagent_file and Path(custom_reagent_file).exists():
                try:
                    self.reagent_abbr.loadFromFile(custom_reagent_file)
                except Exception as e:
                    pass

            self.reagent_abbr.with_(Abbreviations.Option.ALLOW_SINGLETON).with_(
                Abbreviations.Option.AUTO_CONTRACT_TERMINAL
            ).without(Abbreviations.Option.AUTO_CONTRACT_HETERO)

            self._initialized = True

        except Exception as e:
            raise RuntimeError(f"Abbreviation initialization failed: {e}") from e

    def apply(
        self,
        molecule: any,
        mode: AbbreviationMode = AbbreviationMode.REAGENTS,
        highlighted_atoms: Optional[Set[int]] = None,
    ) -> None:
        """Apply abbreviations to a molecule.

        Args:
            molecule: CDK IAtomContainer
            mode: Abbreviation mode (off, groups, reagents, all)
            highlighted_atoms: Set of atom indices that should not be abbreviated
        """
        if not self._initialized:
            self.initialize()

        if mode == AbbreviationMode.OFF:
            return

        # Create atom set for highlighting (atoms to preserve)
        HashMap = JClass("java.util.HashMap")
        atom_set = HashMap()

        if highlighted_atoms:
            for atom_idx in highlighted_atoms:
                if atom_idx < molecule.getAtomCount():
                    atom = molecule.getAtom(atom_idx)
                    atom_set.put(atom, 2)  # Mark as highlighted

        try:
            if mode == AbbreviationMode.ALL:
                # Apply both groups and reagents
                self._contract_hydrates(molecule)
                self.group_abbr.apply(molecule, atom_set)

            elif mode == AbbreviationMode.GROUPS:
                # Only group abbreviations
                self._contract_hydrates(molecule)
                self.group_abbr.apply(molecule, atom_set)

            elif mode == AbbreviationMode.REAGENTS:
                # Only reagent abbreviations
                self._contract_hydrates(molecule)
                self.reagent_abbr.apply(molecule, atom_set)

        except Exception as e:
            pass

    def apply_to_reaction(
        self,
        reaction: any,
        mode: AbbreviationMode = AbbreviationMode.REAGENTS,
        highlighted_atoms: Optional[Set[int]] = None,
    ) -> None:
        """Apply abbreviations to a reaction.

        Args:
            reaction: CDK IReaction
            mode: Abbreviation mode
            highlighted_atoms: Set of atom indices to preserve
        """
        if not self._initialized:
            self.initialize()

        if mode == AbbreviationMode.OFF:
            return

        try:
            # Apply to reactants
            for mol in reaction.getReactants().atomContainers():
                self.apply(mol, mode, highlighted_atoms)

            # Apply to products
            for mol in reaction.getProducts().atomContainers():
                self.apply(mol, mode, highlighted_atoms)

            # Apply to agents (typically more aggressive abbreviation)
            for mol in reaction.getAgents().atomContainers():
                if mode in [AbbreviationMode.ALL, AbbreviationMode.REAGENTS]:
                    self.apply(mol, AbbreviationMode.REAGENTS, highlighted_atoms)
                else:
                    self.apply(mol, mode, highlighted_atoms)

        except Exception as e:
            pass

    def _contract_hydrates(self, molecule: any) -> None:
        """Contract water molecules into hydrate notation (·nH2O).

        This finds isolated water molecules and groups them with Sgroups.

        Args:
            molecule: CDK IAtomContainer
        """
        try:
            CDKConstants = JClass(self.cdk_base + ".CDKConstants")
            Sgroup = JClass(self.cdk_base + ".sgroup.Sgroup")
            SgroupKey = JClass(self.cdk_base + ".sgroup.SgroupKey")
            SgroupType = JClass(self.cdk_base + ".sgroup.SgroupType")
            ArrayList = JClass("java.util.ArrayList")
            HashSet = JClass("java.util.HashSet")

            # Get or create Sgroup list
            sgroups = molecule.getProperty(CDKConstants.CTAB_SGROUPS)
            if sgroups is None:
                sgroups = ArrayList()
                molecule.setProperty(CDKConstants.CTAB_SGROUPS, sgroups)

            # Find water molecules
            hydrate = ArrayList()

            for atom in molecule.atoms():
                # Check if atom is oxygen
                if atom.getAtomicNumber() != 8:  # Not oxygen
                    continue

                # Check if it has exactly 2 hydrogens
                h_count = atom.getImplicitHydrogenCount()
                if h_count is None:
                    h_count = 0

                # Count explicit hydrogens
                explicit_h = 0
                for neighbor in molecule.getConnectedAtomsList(atom):
                    if neighbor.getAtomicNumber() == 1:  # Hydrogen
                        explicit_h += 1

                total_h = h_count + explicit_h

                # Check if isolated water (O with 2H and no other bonds)
                if total_h == 2 and atom.getBondCount() == explicit_h:
                    hydrate.add(atom)
                    # Add explicit hydrogens
                    for neighbor in molecule.getConnectedAtomsList(atom):
                        if neighbor.getAtomicNumber() == 1:
                            hydrate.add(neighbor)

            if hydrate.isEmpty():
                return

            # Check if we should add to existing Sgroup or create new one
            existing_sgroup = None

            for sgrp in sgroups:
                if sgrp.getType() == SgroupType.CtabMultipleGroup:
                    atoms = sgrp.getAtoms()
                    # Check if no overlap
                    okay = True
                    for atom in hydrate:
                        if atoms.contains(atom):
                            okay = False
                            break

                    if okay and (
                        atoms.size() + hydrate.size() == molecule.getAtomCount()
                    ):
                        existing_sgroup = sgrp
                        break

            if existing_sgroup is not None:
                # Add to existing group
                for atom in hydrate:
                    existing_sgroup.addAtom(atom)

                subscript = existing_sgroup.getSubscript()
                new_subscript = subscript + "·" + str(hydrate.size() // 3) + "H2O"
                existing_sgroup.setSubscript(new_subscript)

            else:
                # Create new Sgroup
                sgrp = Sgroup()
                for atom in hydrate:
                    sgrp.addAtom(atom)

                # Set parent atom list
                parent_set = HashSet()
                parent_set.add(hydrate.iterator().next())
                sgrp.putValue(SgroupKey.CtabParentAtomList, parent_set)

                sgrp.setType(SgroupType.CtabMultipleGroup)
                sgrp.setSubscript(str(hydrate.size() // 3))
                sgroups.add(sgrp)

        except Exception as e:
            pass


def get_abbreviation_mode(mode_str: str) -> AbbreviationMode:
    """Convert string to AbbreviationMode enum.

    Args:
        mode_str: Mode string ("off", "groups", "reagents", "on")

    Returns:
        AbbreviationMode enum value
    """
    mode_str = mode_str.lower().strip()

    if mode_str in ["off", "false", "no"]:
        return AbbreviationMode.OFF
    elif mode_str in ["groups", "group"]:
        return AbbreviationMode.GROUPS
    elif mode_str in ["reagents", "reagent"]:
        return AbbreviationMode.REAGENTS
    elif mode_str in ["on", "true", "yes", "all", "both"]:
        return AbbreviationMode.ALL
    else:
        return AbbreviationMode.REAGENTS
