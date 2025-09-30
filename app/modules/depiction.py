from __future__ import annotations

import xml.etree.ElementTree as ET

from jpype import JClass
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

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

    Returns:
        image (SVG): CDK Structure Depiction as an SVG image.
    """
    print(unicolor)

    cdk_base = "org.openscience.cdk"
    StandardGenerator = JClass(
        cdk_base + ".renderer.generators.standard.StandardGenerator",
    )
    Color = JClass("java.awt.Color")
    UniColor = JClass(cdk_base + ".renderer.color.UniColor")
    CDK2DAtomColors = JClass(cdk_base + ".renderer.color.CDK2DAtomColors")()
    Kekulization = JClass(cdk_base + ".aromaticity.Kekulization")
    SmartsPattern = JClass(cdk_base + ".smarts.SmartsPattern")
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")

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

    if CIP:
        SDGMol = get_cip_annotation(molecule)
    else:
        SDGMol = get_CDK_SDG(molecule)

    if SDGMol:
        # Rotate molecule
        if kekulize:
            try:
                Kekulization.kekulize(SDGMol)
            except Exception as e:
                print(str(e) + " Can't Kekulize molecule")

        point = JClass(
            cdk_base + ".geometry.GeometryTools",
        ).get2DCenter(SDGMol)
        JClass(cdk_base + ".geometry.GeometryTools").rotate(
            SDGMol,
            point,
            (rotate * JClass("java.lang.Math").PI / 180.0),
        )

        # Add atom numbers if requested
        if showAtomNumbers:
            DepictionGenerator = DepictionGenerator.withAtomNumbers()

        # Handle highlighting: prioritize atom indices over SMARTS patterns
        if highlight_atoms and len(highlight_atoms) > 0:
            # For CDK, we need to create substructures from atom indices
            # This is more complex and would require additional CDK classes
            # For now, fall back to SMARTS pattern if available
            if highlight and highlight.strip():
                tmpPattern = SmartsPattern.create(highlight, SCOB.getInstance())
                SmartsPattern.prepare(SDGMol)
                tmpMappings = tmpPattern.matchAll(SDGMol)
                tmpSubstructures = tmpMappings.toSubstructures()
                lightBlue = Color(173, 216, 230)
                DepictionGenerator = DepictionGenerator.withHighlight(
                    tmpSubstructures, lightBlue
                ).withOuterGlowHighlight()
            # Note: Direct atom index highlighting in CDK requires more complex implementation
            # This would need creating IAtomContainerSet from specific atoms
        elif highlight and highlight.strip():
            tmpPattern = SmartsPattern.create(highlight, SCOB.getInstance())
            SmartsPattern.prepare(SDGMol)
            tmpMappings = tmpPattern.matchAll(SDGMol)
            tmpSubstructures = tmpMappings.toSubstructures()
            lightBlue = Color(173, 216, 230)
            DepictionGenerator = DepictionGenerator.withHighlight(
                tmpSubstructures, lightBlue
            ).withOuterGlowHighlight()

        mol_imageSVG = (
            DepictionGenerator.depict(
                SDGMol,
            )
            .toSvgStr("px")
            .getBytes()
        )
        encoded_image = ET.tostring(
            ET.fromstring(mol_imageSVG),
            encoding="unicode",
        )

        return encoded_image
    else:
        return "Error reading SMILES string, check again."


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
        for i, atom1_idx in enumerate(hit_ats):
            for j, atom2_idx in enumerate(hit_ats):
                if i < j:  # Avoid duplicate bonds
                    bond = mc.GetBondBetweenAtoms(atom1_idx, atom2_idx)
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
