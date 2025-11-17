from __future__ import annotations

import xml.etree.ElementTree as ET

from jpype import JClass
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from app.modules.cdk_depict.hydrogen_display import set_hydrogen_display

from app.modules.toolkits.cdk_wrapper import get_CDK_SDG
from app.modules.toolkits.cdk_wrapper import get_cip_annotation


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
):
    """This function takes the user input SMILES and Depicts it.

    using the CDK Depiction Generator.

    Args:
        molecule (any): CDK IAtomContainer parsed from SMILES string given by the user.
        molSize (tuple, optional): Size of the output image. Defaults to (512, 512).
        rotate (int, optional): Rotation angle in degrees. Defaults to 0.
        kekulize (bool, optional): Whether to kekulize the molecule. Defaults to True.
        CIP (bool, optional): Whether to annotate CIP stereochemistry. Defaults to True.
        unicolor (bool, optional): Whether to use black and white colors. Defaults to False.
        highlight (str, optional): SMARTS pattern to highlight. Defaults to empty.
        highlight_atoms (list, optional): List of atom indices to highlight. Defaults to None.
        showAtomNumbers (bool, optional): Whether to display atom numbers. Defaults to False.
        hydrogen_display (str, optional): Hydrogen display mode: "Provided", "Minimal",
                                         "Explicit", "Stereo", "Smart". Defaults to "Smart".

    Returns:
        tuple: (encoded_image (str), mol_imageSVG (bytes)) - CDK Structure Depiction as SVG.
               Returns (error_message (str), None) on error.
    """
    cdk_base = "org.openscience.cdk"
    StandardGenerator = JClass(
        cdk_base + ".renderer.generators.standard.StandardGenerator",
    )
    Color = JClass("java.awt.Color")
    UniColor = JClass(cdk_base + ".renderer.color.UniColor")
    CDK2DAtomColors = JClass(cdk_base + ".renderer.color.CDK2DAtomColors")()
    Kekulization = JClass(cdk_base + ".aromaticity.Kekulization")
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")

    # Configure depiction generator based on color scheme
    if unicolor:
        DepictionGenerator = (
            JClass(cdk_base + ".depict.DepictionGenerator")()
            .withSize(molSize[0], molSize[1])
            .withParam(StandardGenerator.StrokeRatio.class_, 1.0)
            .withAnnotationColor(Color.BLACK)
            .withParam(StandardGenerator.AtomColor.class_, UniColor(Color.BLACK))
            .withBackgroundColor(Color.WHITE)
            .withFillToFit()
        )
    else:
        DepictionGenerator = (
            JClass(cdk_base + ".depict.DepictionGenerator")()
            .withAtomColors(CDK2DAtomColors)
            .withSize(molSize[0], molSize[1])
            .withParam(StandardGenerator.StrokeRatio.class_, 1.0)
            .withFillToFit()
            .withBackgroundColor(Color.WHITE)
        )

    # This ensures the hydrogen settings are preserved in the coordinate generation
    try:
        set_hydrogen_display(molecule, hydrogen_display)
    except ValueError as e:
        return f"Invalid hydrogen display mode: {e}", None

    # Generate 2D coordinates with CIP annotation if requested
    if CIP:
        SDGMol = get_cip_annotation(molecule)
    else:
        SDGMol = get_CDK_SDG(molecule)

    if not SDGMol:
        return "Error reading SMILES string, check again.", None

    # Kekulize the molecule if requested
    if kekulize:
        try:
            Kekulization.kekulize(SDGMol)
        except Exception as e:
            print(str(e) + " Can't Kekulize molecule")

    # Rotate molecule if requested
    if rotate != 0:
        point = JClass(cdk_base + ".geometry.GeometryTools").get2DCenter(SDGMol)
        JClass(cdk_base + ".geometry.GeometryTools").rotate(
            SDGMol,
            point,
            (rotate * JClass("java.lang.Math").PI / 180.0),
        )

    # Add atom numbers if requested
    if showAtomNumbers:
        DepictionGenerator = DepictionGenerator.withAtomNumbers()

    # Handle highlighting
    DepictionGenerator = _apply_highlighting(
        DepictionGenerator, SDGMol, highlight, highlight_atoms, cdk_base, SCOB
    )

    # Generate the depiction
    try:
        mol_imageSVG = DepictionGenerator.depict(SDGMol).toSvgStr("px").getBytes()
        encoded_image = ET.tostring(
            ET.fromstring(mol_imageSVG),
            encoding="unicode",
        )
        return encoded_image
    except Exception as e:
        return f"Error generating depiction: {str(e)}", None


def _apply_highlighting(
    depiction_generator, molecule, highlight_pattern, highlight_atoms, cdk_base, scob
):
    """Apply highlighting to the depiction generator.

    Helper function to handle atom/bond highlighting. Prioritizes explicit atom indices
    over SMARTS pattern matching.

    Args:
        depiction_generator: CDK DepictionGenerator instance
        molecule: CDK IAtomContainer with 2D coordinates
        highlight_pattern (str): SMARTS pattern for highlighting
        highlight_atoms (list): List of atom indices or list of lists for multiple substructures
        cdk_base (str): Base CDK package path
        scob: SilentChemObjectBuilder instance

    Returns:
        DepictionGenerator: Updated generator with highlighting configured
    """
    Color = JClass("java.awt.Color")
    lightBlue = Color(173, 216, 230)

    # Priority 1: Explicit atom indices
    if highlight_atoms and len(highlight_atoms) > 0:
        AtomContainer = JClass(cdk_base + ".AtomContainer")
        AtomContainerSet = JClass(cdk_base + ".AtomContainerSet")
        tmpSubstructures = AtomContainerSet()

        # Check if we have multiple substructures (list of lists) or single substructure (flat list)
        is_multiple = isinstance(highlight_atoms[0], (list, tuple))
        atom_index_groups = highlight_atoms if is_multiple else [highlight_atoms]

        for atom_indices in atom_index_groups:
            if len(atom_indices) > 0:
                subset = AtomContainer()

                # Add atoms to subset
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

                tmpSubstructures.addAtomContainer(subset)

        return depiction_generator.withHighlight(
            tmpSubstructures, lightBlue
        ).withOuterGlowHighlight()

    # Priority 2: SMARTS pattern matching
    elif highlight_pattern and highlight_pattern.strip():
        try:
            SmartsPattern = JClass(cdk_base + ".smarts.SmartsPattern")
            tmpPattern = SmartsPattern.create(highlight_pattern, scob.getInstance())
            SmartsPattern.prepare(molecule)
            tmpMappings = tmpPattern.matchAll(molecule)
            tmpSubstructures = tmpMappings.toSubstructures()

            return depiction_generator.withHighlight(
                tmpSubstructures, lightBlue
            ).withOuterGlowHighlight()
        except Exception as e:
            print(f"Warning: Invalid SMARTS pattern '{highlight_pattern}': {e}")
            return depiction_generator

    # No highlighting
    return depiction_generator


def get_rdkit_depiction(
    molecule: Chem.Mol,
    mol_size=(512, 512),
    rotate=0,
    kekulize=True,
    CIP=False,
    unicolor=False,
    highlight: str = "",
    highlight_atoms=None,
    showAtomNumbers=False,
) -> str:
    """
    Generate a 2D depiction of the input molecule using RDKit.

    Args:
        molecule (Chem.Mol): RDKit molecule object.
        mol_size (tuple, optional): Size of the output image. Defaults to (512, 512).
        rotate (int, optional): Rotation angle of the molecule. Defaults to 0.
        kekulize (bool, optional): Whether to kekulize the molecule. Defaults to True.
        CIP (bool, optional): Whether to assign CIP stereochemistry. Defaults to False.
        unicolor (bool, optional): Whether to use a unicolor palette. Defaults to False.
        highlight (str, optional): SMARTS pattern to highlight atoms/bonds. Defaults to empty.
        highlight_atoms (list, optional): List of atom indices to highlight. Defaults to None.
        showAtomNumbers (bool, optional): Whether to display atom numbers. Defaults to False.

    Returns:
        str: RDKit Structure Depiction as an SVG image.
    """

    mc = Chem.Mol(molecule.ToBinary())

    if kekulize:
        try:
            Chem.Kekulize(mc)
        except Chem.KekulizeException:
            mc = Chem.Mol(molecule.ToBinary())

    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)

    if CIP:
        Chem.AssignStereochemistry(mc, force=True, cleanIt=True)

    drawer = rdMolDraw2D.MolDraw2DSVG(mol_size[0], mol_size[1])
    drawer.drawOptions().rotate = rotate
    drawer.drawOptions().addStereoAnnotation = CIP

    if unicolor:
        drawer.drawOptions().useBWAtomPalette()

    # Add atom numbers if requested
    if showAtomNumbers:
        # Set atom numbers as notes on each atom
        for atom in mc.GetAtoms():
            atom.SetProp("atomNote", str(atom.GetIdx()))

    # Handle highlighting based on priority: anchor atoms + SMARTS pattern, then atom indices, then SMARTS pattern alone
    if highlight_atoms and len(highlight_atoms) > 0 and highlight:
        # Combined approach: Use SMARTS pattern but only highlight the match that contains the anchor atoms
        patt = Chem.MolFromSmarts(highlight)
        if patt:
            # Find all matches of the SMARTS pattern
            all_matches = mc.GetSubstructMatches(patt)
            target_match = None

            # Find the match that contains the anchor atoms
            for match in all_matches:
                if any(anchor_atom in match for anchor_atom in highlight_atoms):
                    target_match = match
                    break

            if target_match:
                hit_ats = target_match
                # Find all bonds within this specific match
                hit_bonds = []
                # Check all pairs of atoms in the match to find bonds between them
                for i in range(len(hit_ats)):
                    for j in range(i + 1, len(hit_ats)):
                        bond = mc.GetBondBetweenAtoms(hit_ats[i], hit_ats[j])
                        if bond:
                            hit_bonds.append(bond.GetIdx())

                rdMolDraw2D.PrepareAndDrawMolecule(
                    drawer, mc, highlightAtoms=hit_ats, highlightBonds=hit_bonds
                )
            else:
                # Fallback to just highlighting the anchor atoms if no pattern match contains them
                hit_ats = tuple(highlight_atoms)
                rdMolDraw2D.PrepareAndDrawMolecule(
                    drawer, mc, highlightAtoms=hit_ats, highlightBonds=[]
                )
        else:
            # Invalid SMARTS pattern, fallback to anchor atoms only
            hit_ats = tuple(highlight_atoms)
            rdMolDraw2D.PrepareAndDrawMolecule(
                drawer, mc, highlightAtoms=hit_ats, highlightBonds=[]
            )
    elif highlight_atoms and len(highlight_atoms) > 0:
        # Use specific atom indices for precise highlighting
        hit_ats = tuple(highlight_atoms)
        # Find ALL bonds that connect atoms within the functional group
        hit_bonds = []
        for i in range(len(hit_ats)):
            for j in range(i + 1, len(hit_ats)):
                bond = mc.GetBondBetweenAtoms(hit_ats[i], hit_ats[j])
                if bond:
                    hit_bonds.append(bond.GetIdx())

        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer, mc, highlightAtoms=hit_ats, highlightBonds=hit_bonds
        )
    elif highlight:
        # Fallback to SMARTS pattern highlighting
        patt = Chem.MolFromSmarts(highlight)
        if patt:
            hit_ats = mc.GetSubstructMatch(patt)
            hit_bonds = []
            # Find all bonds within the matched substructure
            for i in range(len(hit_ats)):
                for j in range(i + 1, len(hit_ats)):
                    bond = mc.GetBondBetweenAtoms(hit_ats[i], hit_ats[j])
                    if bond:
                        hit_bonds.append(bond.GetIdx())
            rdMolDraw2D.PrepareAndDrawMolecule(
                drawer, mc, highlightAtoms=hit_ats, highlightBonds=hit_bonds
            )
        else:
            drawer.DrawMolecule(mc)
    else:
        drawer.DrawMolecule(mc)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace("svg:", "")
