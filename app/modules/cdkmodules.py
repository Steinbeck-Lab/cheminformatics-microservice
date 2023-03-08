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
    jar_path = str(pystow.join("STOUT-V2")) + "/cdk-2.8.jar"

    if not os.path.exists(jar_path):
        jar_path = pystow.ensure("STOUT-V2", url=cdk_path)

    startJVM(jvmPath, "-ea", "-Djava.class.path=" + str(jar_path))


def getCDKSDG(smiles: str):
    """This function takes the user input SMILES and Canonicalize it
       using the CDK Canonicalisation algorthim.
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            mol object (string): CDK Structure Diagram Layout mol block.
    """
    cdk_base = "org.openscience.cdk"
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    StringW = JClass("java.io.StringWriter")()
    SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    molecule = SmilesParser.parseSmiles(smiles)
    StructureDiagramGenerator = JClass(cdk_base + ".layout.StructureDiagramGenerator")()
    StructureDiagramGenerator.generateCoordinates(molecule)
    molecule_ = StructureDiagramGenerator.getMolecule()

    SDFW = JClass(cdk_base + ".io.SDFWriter")(StringW)
    SDFW.write(molecule_)
    SDFW.close()
    mol_str = StringW.toString()
    return mol_str
