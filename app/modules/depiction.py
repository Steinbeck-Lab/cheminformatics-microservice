import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from app.modules.toolkits.cdk_wrapper import getCDKSDG, getCIPAnnotation
from jpype import JClass


def getCDKDepiction(
    molecule: any, molSize=(512, 512), rotate=0, CIP=True, unicolor=False
):
    """
    This function takes the user input SMILES and Depicts it
    using the CDK Depiction Generator.

    Args:
        molecule (any): CDK IAtomContainer parsed from SMILES string given by the user.

    Returns:
        image (SVG): CDK Structure Depiction as an SVG image.
    """

    cdk_base = "org.openscience.cdk"
    StandardGenerator = JClass(
        cdk_base + ".renderer.generators.standard.StandardGenerator"
    )
    Color = JClass("java.awt.Color")
    UniColor = JClass(cdk_base + ".renderer.color.UniColor")
    CDK2DAtomColors = JClass(cdk_base + ".renderer.color.CDK2DAtomColors")()
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
        SDGMol = getCIPAnnotation(molecule)
    else:
        SDGMol = getCDKSDG(molecule)

    if SDGMol:
        # Rotate molecule
        point = JClass(cdk_base + ".geometry.GeometryTools").get2DCenter(SDGMol)
        JClass(cdk_base + ".geometry.GeometryTools").rotate(
            SDGMol, point, (rotate * JClass("java.lang.Math").PI / 180.0)
        )

        mol_imageSVG = DepictionGenerator.depict(SDGMol).toSvgStr("px").getBytes()
        encoded_image = ET.tostring(ET.fromstring(mol_imageSVG), encoding="unicode")

        return encoded_image
    else:
        return "Error reading SMILES string, check again."


def getRDKitDepiction(molecule: any, molSize=(512, 512), rotate=0, kekulize=True):
    """
    This function takes the user input SMILES and Canonicalize it using the RDKit.

    Args:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        image (SVG): RDKit Structure Depiction as an SVG image.
    """
    if molecule:
        mc = Chem.Mol(molecule.ToBinary())
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except Exception as e:
                print(e, flush=True)
                mc = Chem.Mol(molecule.ToBinary())
        if not mc.GetNumConformers():
            rdDepictor.Compute2DCoords(mc)
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
        drawer.drawOptions().rotate = rotate
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace("svg:", "")
    else:
        return "Error reading SMILES string, check again."
