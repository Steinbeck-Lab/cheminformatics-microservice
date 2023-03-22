import os
import pystow
from jpype import startJVM, getDefaultJVMPath
from jpype import JClass, JVMNotFoundException, isJVMStarted

# Start JVM to use CDK in python
try:
    jvmPath = getDefaultJVMPath()
except JVMNotFoundException:
    print(
        "If you see this message, for some reason JPype",
        "cannot find jvm.dll.",
        "This indicates that the environment varibale JAVA_HOME",
        "is not set properly.",
        "You can set it or set it manually in the code",
    )
    jvmPath = "Define/path/or/set/JAVA_HOME/variable/properly"
if not isJVMStarted():
    cdk_path = "https://github.com/cdk/cdk/releases/download/cdk-2.8/cdk-2.8.jar"
    sru_path = "https://github.com/JonasSchaub/SugarRemoval/releases/download/v1.3.2/SugarRemovalUtility-jar-with-dependencies.jar"
    cdkjar_path = str(pystow.join("STOUT-V2")) + "/cdk-2.8.jar"
    srujar_path = (
        str(pystow.join("STOUT-V2")) + "/SugarRemovalUtility-jar-with-dependencies.jar"
    )

    if not os.path.exists(cdkjar_path):
        jar_path = pystow.ensure("STOUT-V2", url=cdk_path)

    if not os.path.exists(srujar_path):
        jar_path = pystow.ensure("STOUT-V2", url=sru_path)

    startJVM("-ea", classpath=[cdkjar_path, srujar_path])


def getCDKSDG(smiles: str):
    """This function takes the user input SMILES and Creates a
       Structure Diagram Layout using the CDK.
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            mol object : mol object with CDK SDG.
    """
    cdk_base = "org.openscience.cdk"
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    molecule = SmilesParser.parseSmiles(smiles)
    StructureDiagramGenerator = JClass(cdk_base + ".layout.StructureDiagramGenerator")()
    StructureDiagramGenerator.generateCoordinates(molecule)
    molecule_ = StructureDiagramGenerator.getMolecule()

    return molecule_


def getSugarInfo(smiles: str):
    """This function uses the sugar removal utility and checks
    whether a molecule has ring or linear sugars
    Args:
        smiles (string): SMILES string given by the user.
    Returns:
        (boolean): True or false values whtehr or not molecule has sugar.
    """
    cdk_base = "org.openscience.cdk"
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    molecule = SmilesParser.parseSmiles(smiles)

    sru_base = "de.unijena.cheminf.deglycosylation"

    SugarRemovalUtility = JClass(sru_base + ".SugarRemovalUtility")(SCOB.getInstance())
    hasCircularOrLinearSugars = SugarRemovalUtility.hasCircularOrLinearSugars(molecule)

    if hasCircularOrLinearSugars:
        hasLinearSugar = SugarRemovalUtility.hasLinearSugars(molecule)
        hasCircularSugars = SugarRemovalUtility.hasCircularSugars(molecule)
        return hasLinearSugar, hasCircularSugars
    else:
        return (False, False)


def getMurkoFramework(smiles: str):
    """This function takes the user input SMILES and returns
    the Murko framework
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            smiles (string) : Murko Framework as SMILES.
    """
    cdk_base = "org.openscience.cdk"
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    MurkoFragmenter = JClass(cdk_base + ".fragment.MurckoFragmenter")(True, 3)
    SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    molecule = SmilesParser.parseSmiles(smiles)
    MurkoFragmenter.generateFragments(molecule)
    return str(MurkoFragmenter.getFrameworks()[0])


def getCDKSDGMol(smiles: str):
    """This function takes the user input SMILES and returns a mol
       block as a string with Structure Diagram Layout.
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            mol object (string): CDK Structure Diagram Layout mol block.
    """
    cdk_base = "org.openscience.cdk"
    StringW = JClass("java.io.StringWriter")()

    moleculeSDG = getCDKSDG(smiles)
    SDFW = JClass(cdk_base + ".io.SDFWriter")(StringW)
    SDFW.write(moleculeSDG)
    SDFW.flush()
    mol_str = str(StringW.toString())
    return mol_str
