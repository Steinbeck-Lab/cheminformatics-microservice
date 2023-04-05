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
    cdk_base = "org.openscience.cdk"


def getCDKSDG(smiles: str):
    """This function takes the user input SMILES and Creates a
       Structure Diagram Layout using the CDK.
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            mol object : mol object with CDK SDG.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
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
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
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
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
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
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    StringW = JClass("java.io.StringWriter")()

    moleculeSDG = getCDKSDG(smiles)
    SDFW = JClass(cdk_base + ".io.SDFWriter")(StringW)
    SDFW.write(moleculeSDG)
    SDFW.flush()
    mol_str = str(StringW.toString())
    return mol_str


def getAromaticRingCount(mol):
    """This function is adapted from CDK to
    calculate the number of Aromatic Rings
    present in a given molecule.
    Args (mol): CDK mol object as input.
    Returns (int): Number if aromatic rings present
    """
    Cycles = JClass(cdk_base + ".graph.Cycles")
    ElectronDonation = JClass(cdk_base + ".aromaticity.ElectronDonation")

    Aromaticity = JClass(cdk_base + ".aromaticity.Aromaticity")(
        ElectronDonation.daylight(), Cycles.cdkAromaticSet()
    )
    Aromaticity.apply(mol)
    MCBRings = Cycles.mcb(mol).toRingSet()
    NumberOfAromaticRings = 0
    for RingContainer in MCBRings.atomContainers():
        AreAllRingBondsAromatic = True
        for Bond in RingContainer.bonds():
            if not Bond.isAromatic():
                AreAllRingBondsAromatic = False
                break
        if AreAllRingBondsAromatic:
            NumberOfAromaticRings += 1
    return NumberOfAromaticRings


def getCDKDescriptors(smiles: str):
    """Take an input SMILES and generate a selected set of molecular
    descriptors generated using CDK as a list.
    Args (str): SMILES string
    Returns (list): a list of calculated descriptors
    """
    Mol = getCDKSDG(smiles)
    if Mol:
        AtomCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.AtomCountDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        HeavyAtomsC = Mol.getAtomCount()
        WeightDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.WeightDescriptor")()
            .calculate(Mol)
            .getValue()
            .toString()
        )
        TotalExactMass = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator"
        ).getTotalExactMass(Mol)
        ALogP = (
            JClass(cdk_base + ".qsar.descriptors.molecular.ALOGPDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        NumRotatableBonds = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.RotatableBondsCountDescriptor"
            )()
            .calculate(Mol)
            .getValue()
        )
        TPSADescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.TPSADescriptor")()
            .calculate(Mol)
            .getValue()
            .toString()
        )
        HBondAcceptorCountDescriptor = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.HBondAcceptorCountDescriptor"
            )()
            .calculate(Mol)
            .getValue()
        )
        HBondDonorCountDescriptor = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.HBondAcceptorCountDescriptor"
            )()
            .calculate(Mol)
            .getValue()
        )
        RuleOfFiveDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.RuleOfFiveDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        AromaticRings = getAromaticRingCount(Mol)
        QEDWeighted = None
        FormalCharge = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator"
        ).getTotalFormalCharge(Mol)
        FractionalCSP3Descriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.FractionalCSP3Descriptor")()
            .calculate(Mol)
            .getValue()
            .toString()
        )
        NumRings = JClass(cdk_base + ".graph.Cycles").mcb(Mol).numberOfCycles()

        return (
            int(str(AtomCountDescriptor)),
            int(HeavyAtomsC),
            float("{:.2f}".format(float(str(WeightDescriptor)))),
            float("{:.5f}".format(float(str(TotalExactMass)))),
            float("{:.2f}".format(float(str(ALogP).split(",")[0]))),
            int(str(NumRotatableBonds)),
            float("{:.2f}".format(float(str(TPSADescriptor)))),
            int(str(HBondAcceptorCountDescriptor)),
            int(str(HBondDonorCountDescriptor)),
            int(str(HBondAcceptorCountDescriptor)),
            int(str(HBondDonorCountDescriptor)),
            int(str(RuleOfFiveDescriptor)),
            int(AromaticRings),
            str(QEDWeighted),
            int(FormalCharge),
            float("{:.2f}".format(float(str(FractionalCSP3Descriptor)))),
            int(NumRings),
        )
