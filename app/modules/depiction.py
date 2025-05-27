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
    format="svg",
):
    """This function takes the user input SMILES and Depicts it using the CDK Depiction Generator.

    Args:
        molecule (any): CDK IAtomContainer parsed from SMILES string given by the user.
        molSize (tuple, optional): Size of the output image. Defaults to (512, 512).
        rotate (int, optional): Rotation angle in degrees. Defaults to 0.
        kekulize (bool, optional): Whether to kekulize the molecule. Defaults to True.
        CIP (bool, optional): Whether to include CIP stereochemistry. Defaults to True.
        unicolor (bool, optional): Whether to use a single color. Defaults to False.
        highlight (str, optional): SMARTS pattern to highlight. Defaults to "".
        format (str, optional): Output format ("svg" or "png"). Defaults to "svg".

    Returns:
        str: CDK Structure Depiction as an SVG or PNG image.
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
                print(e + "Can't Kekulize molecule")

        point = JClass(
            cdk_base + ".geometry.GeometryTools",
        ).get2DCenter(SDGMol)
        JClass(cdk_base + ".geometry.GeometryTools").rotate(
            SDGMol,
            point,
            (rotate * JClass("java.lang.Math").PI / 180.0),
        )

        if highlight and highlight.strip():
            tmpPattern = SmartsPattern.create(highlight, SCOB.getInstance())
            SmartsPattern.prepare(SDGMol)
            tmpMappings = tmpPattern.matchAll(SDGMol)
            tmpSubstructures = tmpMappings.toSubstructures()
            lightBlue = Color(173, 216, 230)
            DepictionGenerator = DepictionGenerator.withHighlight(
                tmpSubstructures, lightBlue
            ).withOuterGlowHighlight()

        if format == "svg":
            mol_image = (
                DepictionGenerator.depict(
                    SDGMol,
                )
                .toSvgStr("px")
                .getBytes()
            )
            encoded_image = ET.tostring(
                ET.fromstring(mol_image),
                encoding="unicode",
            )
        else:  # PNG format
            mol_image = (
                DepictionGenerator.depict(
                    SDGMol,
                )
                .toPngStr()
                .getBytes()
            )
            encoded_image = mol_image

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
    format: str = "svg",
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
        format (str, optional): Output format ("svg" or "png"). Defaults to "svg".

    Returns:
        str: RDKit Structure Depiction as an SVG or PNG image.
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

    if format == "svg":
        drawer = rdMolDraw2D.MolDraw2DSVG(mol_size[0], mol_size[1])
    else:  # PNG format
        drawer = rdMolDraw2D.MolDraw2DCairo(mol_size[0], mol_size[1])

    drawer.drawOptions().rotate = rotate
    drawer.drawOptions().addStereoAnnotation = CIP

    if unicolor:
        drawer.drawOptions().useBWAtomPalette()

    if highlight:
        patt = Chem.MolFromSmarts(highlight)
        if patt:
            hit_ats = mc.GetSubstructMatch(patt)
            hit_bonds = [
                mc.GetBondBetweenAtoms(at1, at2).GetIdx()
                for at1, at2 in zip(hit_ats[:-1], hit_ats[1:])
            ]
            rdMolDraw2D.PrepareAndDrawMolecule(
                drawer, mc, highlightAtoms=hit_ats, highlightBonds=hit_bonds
            )
        else:
            drawer.DrawMolecule(mc)
    else:
        drawer.DrawMolecule(mc)

    drawer.FinishDrawing()

    if format == "svg":
        svg = drawer.GetDrawingText()
        return svg.replace("svg:", "")
    else:  # PNG format
        return drawer.GetDrawingText()
