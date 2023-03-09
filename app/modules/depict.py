from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import app.modules.cdkmodules as cdkmodules
from IPython.display import SVG


def getCDKDepiction(smiles: str, size=512.0):
    """This function takes the user input SMILES and Depicts it
       using the CDK Depiction Generator.
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            imag (PIL): CDK Structure Depiction as a pillow image.
            image (png): CDK Structure Depiction as a PNG image.
    """
    cdk_base = "org.openscience.cdk"
    StandardGenerator = cdkmodules.JClass(
        cdk_base + ".renderer.generators.standard.StandardGenerator"
    )
    DepictionGenerator = cdkmodules.JClass(cdk_base + ".depict.DepictionGenerator")()
    Color = cdkmodules.JClass("java.awt.Color")
    UniColor = cdkmodules.JClass(cdk_base + ".renderer.color.UniColor")

    DepictionGenerator.withSize(size, size).withAtomValues().withParam(
        StandardGenerator.StrokeRatio.class_, 1.5
    ).withAnnotationColor(Color.BLACK).withParam(
        StandardGenerator.AtomColor.class_, UniColor(Color.BLACK)
    ).withBackgroundColor(
        Color.WHITE
    ).withZoom(
        2.0
    )

    moleculeSDG = cdkmodules.getCDKSDG(smiles)
    mol_image = DepictionGenerator.depict(moleculeSDG).toSvgStr("px")

    return mol_image.getBytes()


def getRDKitDepiction(smiles, molSize=(512, 512), kekulize=True):
    """This function takes the user input SMILES and Canonicalize it
       using the RDKit.
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            imag (PIL): CDK Structure Depiction as a pillow image.
            image (png): CDK Structure Depiction as a PNG image.
    """
    mol = Chem.MolFromSmiles(smiles)
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return SVG(svg.replace("svg:", ""))
