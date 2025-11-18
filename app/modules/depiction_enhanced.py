"""Enhanced Molecular Depiction Module using CDK.

This module provides advanced 2D molecular depiction functionality using the
Chemistry Development Kit (CDK) exclusively. It supports:

- CXSMILES highlighting (atoms and bonds)
- Chemical abbreviations (functional groups and reagents)
- Dative bond perception for coordination chemistry
- Multicenter bond display for η-complexes
- Multiple style presets and color schemes
- Comprehensive annotation modes
- Aromatic donut display
- Radical perception
- MDL HILITE support
- Advanced rendering controls (zoom, ratio, flip, etc.)

Author: Kohulan Rajan
License: MIT
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

from jpype import JClass

from app.modules.toolkits.cdk_wrapper import get_CDK_SDG, get_cip_annotation

from app.modules.cdk_depict import (
    set_hydrogen_display,
    ChemicalAbbreviations,
    AbbreviationMode,
    DativeBondPerception,
    DativeBondMode,
    MulticenterBonds,
    MulticenterStyle,
    AnnotationSystem,
    AnnotationMode,
    AromaticDisplaySystem,
    StylePresetSystem,
    StylePreset,
    AdvancedControls,
    SVGUnits,
    ReactionArrowSystem,
    ReactionArrowType,
    RadicalPerception,
    MDLHiliteParser,
)

from app.modules.cdk_depict.cxsmiles_parser import (
    parse_cxsmiles_highlighting_from_string,
)


def get_cdk_depiction(
    molecule: any,
    molSize=(512, 512),
    rotate=0,
    kekulize=True,
    CIP=True,
    unicolor=False,
    highlight="",
    highlight_atoms=None,
    showAtomNumbers=False,
    hydrogen_display="Smart",
    # Chemical structure enhancements
    abbreviate="reagents",
    dative="metals",
    multicenter="provided",
    # Visual style and annotations
    annotate="none",
    style="cow",
    donuts=False,
    arrow="",
    alignrxnmap=True,
    showtitle=False,
    bgcolor=None,
    fgcolor=None,
    # Advanced rendering controls
    zoom=1.3,
    ratio=1.0,
    flip=False,
    anon=False,
    smalim=100,
    svgunits="px",
    perceive_radicals=False,
    apply_mdl_highlighting=True,
    smiles_string=None,
):
    """Generate enhanced 2D molecular depiction using CDK with all features.

    Args:
        molecule: CDK IAtomContainer (pre-parsed from SMILES/CXSMILES)
        smiles_string: Optional original SMILES/CXSMILES string for parsing highlighting
        molSize: Output image size
        rotate: Rotation angle in degrees
        kekulize: Whether to kekulize the molecule
        CIP: Whether to annotate CIP stereochemistry
        unicolor: Whether to use black and white colors
        highlight: SMARTS pattern to highlight
        highlight_atoms: List of atom indices to highlight
        showAtomNumbers: Whether to display atom numbers
        hydrogen_display: Hydrogen display mode (Provided/Minimal/Explicit/Stereo/Smart)

        # Chemical structure enhancements
        abbreviate: Abbreviation mode (off/groups/reagents/on)
        dative: Dative bond mode (never/metals/always)
        multicenter: Multicenter bond style

        # Visual style and annotations
        annotate: Annotation mode
        style: Style preset (cow/cob/bow/nob/etc)
        donuts: Whether to use aromatic circles
        arrow: Reaction arrow type
        alignrxnmap: Whether to align reaction mapping
        showtitle: Whether to display title
        bgcolor: Background color
        fgcolor: Foreground color

        # Advanced rendering controls
        zoom: Zoom level (0.1-5.0)
        ratio: Stroke ratio (0.5-2.0)
        flip: Whether to flip structure
        anon: Whether to use anonymous display
        smalim: SMARTS hit limit
        svgunits: SVG units (px/mm/cm/in)
        perceive_radicals: Whether to detect radicals
        apply_mdl_highlighting: Whether to apply MDL HILITE

    Returns:
        str: SVG string of molecular depiction
    """
    cdk_base = "org.openscience.cdk"

    try:
        # Initialize CDK classes
        StandardGenerator = JClass(
            cdk_base + ".renderer.generators.standard.StandardGenerator"
        )
        Color = JClass("java.awt.Color")
        UniColor = JClass(cdk_base + ".renderer.color.UniColor")
        CDK2DAtomColors = JClass(cdk_base + ".renderer.color.CDK2DAtomColors")()
        Kekulization = JClass(cdk_base + ".aromaticity.Kekulization")
        SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")

        # Determine if this is a reaction
        IReaction = JClass(cdk_base + ".interfaces.IReaction")
        IReactionSet = JClass(cdk_base + ".interfaces.IReactionSet")
        is_reaction = isinstance(molecule, (IReaction, IReactionSet))

        # ====== Extract CXSMILES highlighting from original SMILES string ======
        cx_highlight_atoms, cx_highlight_bonds = set(), set()
        if smiles_string:
            (
                cx_highlight_atoms,
                cx_highlight_bonds,
            ) = parse_cxsmiles_highlighting_from_string(smiles_string)

        # ====== Radical Perception ======
        if perceive_radicals:
            if is_reaction:
                # For reactions, apply to all components (reactants, products, agents)
                if isinstance(molecule, IReactionSet):
                    for rxn in molecule.reactions():
                        # Process reactants
                        for mol in rxn.getReactants().atomContainers():
                            RadicalPerception.perceive_radicals(mol)
                        # Process products
                        for mol in rxn.getProducts().atomContainers():
                            RadicalPerception.perceive_radicals(mol)
                        # Process agents
                        for mol in rxn.getAgents().atomContainers():
                            RadicalPerception.perceive_radicals(mol)
                else:  # Single IReaction
                    # Process reactants
                    for mol in molecule.getReactants().atomContainers():
                        RadicalPerception.perceive_radicals(mol)
                    # Process products
                    for mol in molecule.getProducts().atomContainers():
                        RadicalPerception.perceive_radicals(mol)
                    # Process agents
                    for mol in molecule.getAgents().atomContainers():
                        RadicalPerception.perceive_radicals(mol)
            else:
                # For molecules, apply directly
                RadicalPerception.perceive_radicals(molecule)

        # ====== MDL HILITE Support ======
        mdl_highlight_atoms = set()
        mdl_highlight_bonds = set()
        if apply_mdl_highlighting and not is_reaction:
            try:
                mdl_parser = MDLHiliteParser()
                mdl_atoms, mdl_bonds = mdl_parser.parse_hilite_from_mol(molecule)
                mdl_highlight_atoms = mdl_atoms
                mdl_highlight_bonds = mdl_bonds
                if mdl_atoms or mdl_bonds:
                    pass
            except Exception:
                pass

        # ====== Merge all highlighting sources ======
        # Priority: explicit highlight_atoms > CXSMILES > MDL HILITE
        if highlight_atoms is None:
            highlight_atoms = list(cx_highlight_atoms | mdl_highlight_atoms)
        else:
            # User provided explicit atoms - merge with CXSMILES and MDL
            highlight_atoms = list(
                set(highlight_atoms) | cx_highlight_atoms | mdl_highlight_atoms
            )

        # Merge bond highlighting from all sources
        highlight_bonds = list(cx_highlight_bonds | mdl_highlight_bonds)

        # ====== Configure Depiction Generator ======
        # Initialize base generator
        DepictionGenerator = JClass(cdk_base + ".depict.DepictionGenerator")()

        # Configure unicolor BEFORE style preset (takes priority over style presets)
        if unicolor:
            pass
            DepictionGenerator = (
                DepictionGenerator.withSize(molSize[0], molSize[1])
                .withParam(StandardGenerator.StrokeRatio.class_, 1.0)
                .withAnnotationColor(Color.BLACK)
                .withParam(StandardGenerator.AtomColor.class_, UniColor(Color.BLACK))
                .withBackgroundColor(Color.WHITE)
                .withFillToFit()
            )
        else:
            # ====== Style Preset Configuration ======

            # Apply style preset
            style_system = StylePresetSystem()
            try:
                style_enum = StylePreset(style.lower())
                DepictionGenerator = style_system.apply_preset(
                    DepictionGenerator, style_enum
                )
            except ValueError:
                # Apply default colors if style preset fails
                DepictionGenerator = (
                    DepictionGenerator.withAtomColors(CDK2DAtomColors)
                    .withSize(molSize[0], molSize[1])
                    .withParam(StandardGenerator.StrokeRatio.class_, 1.0)
                    .withFillToFit()
                    .withBackgroundColor(Color.WHITE)
                )

            # Override with custom colors if provided (only when not in unicolor mode)
            if bgcolor or fgcolor:
                if bgcolor:
                    try:
                        bg_color = (
                            Color.decode(bgcolor)
                            if bgcolor != "default"
                            else Color.WHITE
                        )
                        DepictionGenerator = DepictionGenerator.withBackgroundColor(
                            bg_color
                        )
                    except Exception:
                        pass
                if fgcolor:
                    try:
                        fg_color = (
                            Color.decode(fgcolor)
                            if fgcolor != "default"
                            else Color.BLACK
                        )
                        DepictionGenerator = DepictionGenerator.withAnnotationColor(
                            fg_color
                        )
                    except Exception:
                        pass

            # Apply size
            DepictionGenerator = DepictionGenerator.withSize(molSize[0], molSize[1])

        # ====== Zoom and Stroke Ratio ======
        controls = AdvancedControls()
        DepictionGenerator = controls.set_zoom(DepictionGenerator, zoom)
        DepictionGenerator = controls.set_stroke_ratio(DepictionGenerator, ratio)

        # ====== Map Alignment (for reactions) ======
        if is_reaction and alignrxnmap:
            DepictionGenerator = DepictionGenerator.withMappedRxnAlign()

        # ====== Anonymous Display ======
        if anon:
            SymbolVisibility = JClass(cdk_base + ".renderer.SymbolVisibility")
            DepictionGenerator = DepictionGenerator.withParam(
                StandardGenerator.Visibility.class_,
                SymbolVisibility.iupacRecommendations(),
            )

        DepictionGenerator = DepictionGenerator.withFillToFit()

        # ====== Abbreviations ======
        if abbreviate and abbreviate.lower() != "off":
            try:
                abbr_system = ChemicalAbbreviations(cdk_base)
                abbr_system.initialize()

                abbr_mode_map = {
                    "off": AbbreviationMode.OFF,
                    "groups": AbbreviationMode.GROUPS,
                    "reagents": AbbreviationMode.REAGENTS,
                    "on": AbbreviationMode.ALL,
                }
                abbr_mode = abbr_mode_map.get(
                    abbreviate.lower(), AbbreviationMode.REAGENTS
                )

                highlight_set = set(highlight_atoms) if highlight_atoms else set()
                abbr_system.apply(
                    molecule, mode=abbr_mode, highlighted_atoms=highlight_set
                )
            except Exception:
                pass

        # ====== Dative Bonds ======
        if dative and dative.lower() != "never":
            try:
                dative_perceiver = DativeBondPerception(cdk_base)

                dative_mode_map = {
                    "never": DativeBondMode.NEVER,
                    "metals": DativeBondMode.METALS,
                    "always": DativeBondMode.ALWAYS,
                }
                dative_mode = dative_mode_map.get(dative.lower(), DativeBondMode.METALS)

                count = dative_perceiver.perceive(molecule, mode=dative_mode)
                if count > 0:
                    pass
            except Exception:
                pass

        # ====== Hydrogen Display ======
        try:
            set_hydrogen_display(
                molecule if not is_reaction else None, hydrogen_display
            )
        except ValueError:
            pass

        # ====== Generate 2D Coordinates with CIP if requested ======
        if CIP:
            SDGMol = get_cip_annotation(molecule)
        else:
            SDGMol = get_CDK_SDG(molecule)

        if not SDGMol:
            return "<svg><text>Error reading SMILES string</text></svg>"

        # ====== Reaction Arrow Type ======
        if is_reaction and arrow:
            try:
                arrow_system = ReactionArrowSystem()
                arrow_map = {
                    "forward": ReactionArrowType.FORWARD,
                    "equ": ReactionArrowType.BIDIRECTIONAL,
                    "ngo": ReactionArrowType.NO_GO,
                    "ret": ReactionArrowType.RETRO_SYNTHETIC,
                    "res": ReactionArrowType.RESONANCE,
                }
                arrow_type = arrow_map.get(arrow.lower(), ReactionArrowType.FORWARD)
                arrow_system.set_arrow_type(SDGMol, arrow_type)
            except Exception:
                pass

        # ====== Multicenter Bonds ======
        if multicenter and multicenter.lower() != "provided":
            try:
                mc_handler = MulticenterBonds(cdk_base)

                mc_style_map = {
                    "provided": MulticenterStyle.PROVIDED,
                    "dative": MulticenterStyle.DATIVE,
                    "dashed": MulticenterStyle.DASHED,
                    "dashed_neutral": MulticenterStyle.DASHED_NEUTRAL,
                    "hidden": MulticenterStyle.HIDDEN,
                    "hidden_neutral": MulticenterStyle.HIDDEN_NEUTRAL,
                }
                mc_style = mc_style_map.get(
                    multicenter.lower(), MulticenterStyle.PROVIDED
                )

                count = mc_handler.set_style(SDGMol, style=mc_style)
                if count > 0:
                    pass
            except Exception:
                pass

        # ====== Kekulize ======
        if kekulize and not is_reaction:
            try:
                Kekulization.kekulize(SDGMol)
            except Exception:
                pass

        # ====== Aromatic Display (Donuts) ======
        if not is_reaction:
            aromatic_system = AromaticDisplaySystem()
            DepictionGenerator = aromatic_system.apply_donut_display(
                DepictionGenerator, SDGMol, enable=donuts
            )

        # ====== Flip Structure ======
        if flip:
            try:
                controls.flip_structure(SDGMol)
            except Exception:
                pass

        # ====== Rotate ======
        if rotate != 0:
            point = JClass(cdk_base + ".geometry.GeometryTools").get2DCenter(SDGMol)
            JClass(cdk_base + ".geometry.GeometryTools").rotate(
                SDGMol,
                point,
                (rotate * JClass("java.lang.Math").PI / 180.0),
            )

        # ====== Annotations ======
        if annotate != "none" or showAtomNumbers:
            annotation_mode = "number" if showAtomNumbers else annotate

            annotation_map = {
                "none": AnnotationMode.NONE,
                "number": AnnotationMode.NUMBER,
                "bondnumber": AnnotationMode.BONDNUMBER,
                "mapidx": AnnotationMode.MAPIDX,
                "atomvalue": AnnotationMode.ATOMVALUE,
                "colmap": AnnotationMode.COLMAP,
                "cip": AnnotationMode.CIP,
            }
            annot_enum = annotation_map.get(
                annotation_mode.lower(), AnnotationMode.NONE
            )

            annotation_system = AnnotationSystem()
            DepictionGenerator = annotation_system.apply_annotations(
                DepictionGenerator, SDGMol, annot_enum, is_reaction
            )

        # ====== Title Display ======
        if showtitle:
            if is_reaction:
                DepictionGenerator = DepictionGenerator.withRxnTitle()
            else:
                DepictionGenerator = DepictionGenerator.withMolTitle()

        # ====== Highlighting ======
        DepictionGenerator = _apply_highlighting(
            DepictionGenerator,
            SDGMol,
            highlight,
            highlight_atoms,
            highlight_bonds,
            cdk_base,
            SCOB,
            smalim,
            style,
        )

        # ====== Generate Depiction ======
        units_map = {
            "px": SVGUnits.PX,
            "mm": SVGUnits.MM,
            "cm": SVGUnits.CM,
            "in": SVGUnits.IN,
        }
        units_enum = units_map.get(svgunits.lower(), SVGUnits.PX)
        units_str = units_enum.value

        mol_imageSVG = DepictionGenerator.depict(SDGMol).toSvgStr(units_str).getBytes()
        encoded_image = ET.tostring(
            ET.fromstring(mol_imageSVG),
            encoding="unicode",
        )

        return encoded_image

    except Exception as e:
        return f"<svg><text>Error: {str(e)}</text></svg>"


def _apply_highlighting(
    depiction_generator,
    molecule,
    highlight_pattern,
    highlight_atoms,
    highlight_bonds,
    cdk_base,
    scob,
    smalim=100,
    style="cow",
):
    """Apply highlighting to the depiction generator.

    Args:
        depiction_generator: CDK DepictionGenerator instance
        molecule: CDK IAtomContainer with 2D coordinates
        highlight_pattern: SMARTS pattern for highlighting
        highlight_atoms: List of atom indices to highlight
        highlight_bonds: List of bond indices to highlight (from CXSMILES/MDL)
        cdk_base: Base CDK package path
        scob: SilentChemObjectBuilder instance
        smalim: SMARTS hit limit
        style: Style preset for choosing highlight color

    Returns:
        Updated DepictionGenerator with highlighting applied
    """
    try:
        Color = JClass("java.awt.Color")

        # Choose highlight color based on style
        if style in ["nob"]:
            hgCol = Color(255, 170, 170)
        elif style in ["bow", "wob", "bot"]:
            hgCol = Color.RED
        else:
            hgCol = Color(170, 255, 170)

        # Priority 1: Explicit atom indices
        if highlight_atoms and len(highlight_atoms) > 0:
            AtomContainer = JClass(cdk_base + ".AtomContainer")
            AtomContainerSet = JClass(cdk_base + ".AtomContainerSet")
            tmpSubstructures = AtomContainerSet()

            is_multiple = isinstance(highlight_atoms[0], (list, tuple))
            atom_index_groups = highlight_atoms if is_multiple else [highlight_atoms]

            for atom_indices in atom_index_groups:
                if len(atom_indices) > 0:
                    subset = AtomContainer()

                    for idx in atom_indices:
                        if idx < molecule.getAtomCount():
                            subset.addAtom(molecule.getAtom(idx))

                    # Add bonds between highlighted atoms
                    for i, idx1 in enumerate(atom_indices):
                        for idx2 in atom_indices[i + 1:]:
                            if (
                                idx1 < molecule.getAtomCount()
                                and idx2 < molecule.getAtomCount()
                            ):
                                bond = molecule.getBond(
                                    molecule.getAtom(idx1), molecule.getAtom(idx2)
                                )
                                if bond is not None:
                                    subset.addBond(bond)

                    # Add explicitly highlighted bonds (from CXSMILES/MDL)
                    if highlight_bonds and len(highlight_bonds) > 0:
                        for bond_idx in highlight_bonds:
                            if bond_idx < molecule.getBondCount():
                                bond = molecule.getBond(bond_idx)
                                # Only add if both atoms are in the subset
                                if (
                                    bond.getBegin() in subset.atoms()
                                    and bond.getEnd() in subset.atoms()
                                ):
                                    if bond not in subset.bonds():
                                        subset.addBond(bond)

                    tmpSubstructures.addAtomContainer(subset)

            return depiction_generator.withHighlight(
                tmpSubstructures, hgCol
            ).withOuterGlowHighlight()

        # Priority 2: Only bond highlighting (no atoms) - from CXSMILES/MDL
        elif highlight_bonds and len(highlight_bonds) > 0:
            AtomContainer = JClass(cdk_base + ".AtomContainer")
            AtomContainerSet = JClass(cdk_base + ".AtomContainerSet")
            tmpSubstructures = AtomContainerSet()
            subset = AtomContainer()

            # Add bonds and their atoms
            for bond_idx in highlight_bonds:
                if bond_idx < molecule.getBondCount():
                    bond = molecule.getBond(bond_idx)
                    if bond.getBegin() not in subset.atoms():
                        subset.addAtom(bond.getBegin())
                    if bond.getEnd() not in subset.atoms():
                        subset.addAtom(bond.getEnd())
                    if bond not in subset.bonds():
                        subset.addBond(bond)

            if subset.getAtomCount() > 0:
                tmpSubstructures.addAtomContainer(subset)
                return depiction_generator.withHighlight(
                    tmpSubstructures, hgCol
                ).withOuterGlowHighlight()

        # Priority 3: SMARTS pattern matching
        elif highlight_pattern and highlight_pattern.strip():
            try:
                SmartsPattern = JClass(cdk_base + ".smarts.SmartsPattern")
                tmpPattern = SmartsPattern.create(highlight_pattern, scob.getInstance())
                SmartsPattern.prepare(molecule)
                tmpMappings = tmpPattern.matchAll(molecule)

                # Apply hit limiting
                ArrayList = JClass("java.util.ArrayList")
                limited_mappings = ArrayList()

                count = 0
                iterator = tmpMappings.iterator()
                while iterator.hasNext() and count < smalim:
                    limited_mappings.add(iterator.next())
                    count += 1

                if count > 0:
                    AtomContainerSet = JClass(cdk_base + ".AtomContainerSet")
                    tmpSubstructures = AtomContainerSet()

                    for mapping in limited_mappings:
                        AtomContainer = JClass(cdk_base + ".AtomContainer")
                        subset = AtomContainer()

                        for idx in mapping:
                            subset.addAtom(molecule.getAtom(idx))

                        for i, idx1 in enumerate(mapping):
                            for idx2 in mapping[i + 1:]:
                                bond = molecule.getBond(
                                    molecule.getAtom(idx1), molecule.getAtom(idx2)
                                )
                                if bond is not None:
                                    subset.addBond(bond)

                        tmpSubstructures.addAtomContainer(subset)

                    return depiction_generator.withHighlight(
                        tmpSubstructures, hgCol
                    ).withOuterGlowHighlight()

            except Exception:
                pass

        return depiction_generator

    except Exception:
        return depiction_generator
