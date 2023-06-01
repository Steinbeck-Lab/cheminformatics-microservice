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
    centres_path = (
        "https://github.com/SiMolecule/centres/releases/download/1.0/centres.jar"
    )
    cdkjar_path = str(pystow.join("STOUT-V2")) + "/cdk-2.8.jar"
    srujar_path = (
        str(pystow.join("STOUT-V2")) + "/SugarRemovalUtility-jar-with-dependencies.jar"
    )
    centresjar_path = str(pystow.join("STOUT-V2")) + "/centres.jar"

    if not os.path.exists(cdkjar_path):
        jar_path = pystow.ensure("STOUT-V2", url=cdk_path)

    if not os.path.exists(srujar_path):
        jar_path = pystow.ensure("STOUT-V2", url=sru_path)

    if not os.path.exists(centresjar_path):
        jar_path = pystow.ensure("STOUT-V2", url=centres_path)

    startJVM("-ea", classpath=[cdkjar_path, srujar_path, centresjar_path])
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


def getTanimotoSimilarityCDK(smiles1: str, smiles2: str):
    """
    Take two SMILES strings and calculate
    Tanimoto similarity index using Pubchem
    Fingerprints.
    Args (str,str): SMILES strings.
    Returns (float): Tanimoto similarity.
    """
    if any(char.isspace() for char in smiles1):
        smiles1 = smiles1.replace(" ", "+")
    if any(char.isspace() for char in smiles2):
        smiles2 = smiles2.replace(" ", "+")

    Tanimoto = JClass(cdk_base + ".similarity.Tanimoto")
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    PubchemFingerprinter = JClass(cdk_base + ".fingerprint.PubchemFingerprinter")(
        SCOB.getInstance()
    )

    # parse molecules to get IAtomContainers
    mol1 = SmilesParser.parseSmiles(smiles1)
    mol2 = SmilesParser.parseSmiles(smiles2)

    # Generate BitSets using PubChemFingerprinter
    fingerprint1 = PubchemFingerprinter.getBitFingerprint(mol1).asBitSet()
    fingerprint2 = PubchemFingerprinter.getBitFingerprint(mol2).asBitSet()

    # Calculate Tanimoto similarity
    Similarity = Tanimoto.calculate(fingerprint1, fingerprint2)

    return "{:.5f}".format(float(str(Similarity)))


def getCIPAnnotation(smiles: str):
    """
    The function return the CIP annotations using the CDK
    CIP toolkit.
    Args: mol block
    Returns: CIP annotated mol block
    """
    mol = getCDKSDG(smiles)
    centres_base = "com.simolecule.centres"
    Cycles = JClass(cdk_base + ".graph.Cycles")
    IBond = JClass(cdk_base + ".interfaces.IBond")
    IStereoElement = JClass(cdk_base + ".interfaces.IStereoElement")
    Stereocenters = JClass(cdk_base + ".stereo.Stereocenters")
    StandardGenerator = JClass(
        cdk_base + ".renderer.generators.standard.StandardGenerator"
    )

    BaseMol = JClass(centres_base + ".BaseMol")
    CdkLabeller = JClass(centres_base + ".CdkLabeller")
    Descriptor = JClass(centres_base + ".Descriptor")

    stereocenters = Stereocenters.of(mol)
    for atom in mol.atoms():
        if (
            stereocenters.isStereocenter(atom.getIndex())
            and stereocenters.elementType(atom.getIndex())
            == Stereocenters.Type.Tetracoordinate
        ):
            atom.setProperty(StandardGenerator.ANNOTATION_LABEL, "(?)")

    for bond in mol.bonds():
        if bond.getOrder() != IBond.Order.DOUBLE:
            continue
        begIdx = bond.getBegin().getIndex()
        endIdx = bond.getEnd().getIndex()
        if (
            stereocenters.elementType(begIdx) == Stereocenters.Type.Tricoordinate
            and stereocenters.elementType(endIdx) == Stereocenters.Type.Tricoordinate
            and stereocenters.isStereocenter(begIdx)
            and stereocenters.isStereocenter(endIdx)
        ):
            # only if not in a small ring <7
            if Cycles.smallRingSize(bond, 7) == 0:
                bond.setProperty(StandardGenerator.ANNOTATION_LABEL, "(?)")
    # no defined stereo?
    if not mol.stereoElements().iterator().hasNext():
        return mol

    CdkLabeller.label(mol)

    # update to label appropriately for racmic and relative stereochemistry
    for se in mol.stereoElements():
        if se.getConfigClass() == IStereoElement.TH and se.getGroupInfo() != 0:
            focus = se.getFocus()
            label = focus.getProperty(BaseMol.CIP_LABEL_KEY)
            if (
                isinstance(label, Descriptor)
                and label != Descriptor.ns
                and label != Descriptor.Unknown
            ):
                if (se.getGroupInfo() & IStereoElement.GRP_RAC) != 0:
                    inv = None
                    if label == Descriptor.R:
                        inv = Descriptor.S
                    elif label == Descriptor.S:
                        inv = Descriptor.R
                    if inv is not None:
                        focus.setProperty(
                            BaseMol.CIP_LABEL_KEY, label.toString() + inv.name()
                        )
                elif (se.getGroupInfo() & IStereoElement.GRP_REL) != 0:
                    if label in [Descriptor.R, Descriptor.S]:
                        focus.setProperty(BaseMol.CIP_LABEL_KEY, label.toString() + "*")

    for atom in mol.atoms():
        if atom.getProperty(BaseMol.CONF_INDEX) is not None:
            atom.setProperty(
                StandardGenerator.ANNOTATION_LABEL,
                StandardGenerator.ITALIC_DISPLAY_PREFIX
                + str(atom.getProperty(BaseMol.CONF_INDEX)),
            )
        elif atom.getProperty(BaseMol.CIP_LABEL_KEY) is not None:
            atom.setProperty(
                StandardGenerator.ANNOTATION_LABEL,
                StandardGenerator.ITALIC_DISPLAY_PREFIX
                + str(atom.getProperty(BaseMol.CIP_LABEL_KEY)),
            )
    for bond in mol.bonds():
        if bond.getProperty(BaseMol.CIP_LABEL_KEY) is not None:
            bond.setProperty(
                StandardGenerator.ANNOTATION_LABEL,
                StandardGenerator.ITALIC_DISPLAY_PREFIX
                + bond.getProperty(BaseMol.CIP_LABEL_KEY),
            )

    return mol


def getCXSMILES(smiles: str):
    """This function takes the user input SMILES and creates a
    CXSMILES string with 2D atom coordinates
    Args:
            smiles (string): SMILES string given by the user.
    Returns:
            smiles (string): CXSMILES string.

    """
    moleculeSDG = getCDKSDG(smiles)
    SmiFlavor = JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute | SmiFlavor.CxSmilesWithCoords
    )
    CXSMILES = SmilesGenerator.create(moleculeSDG)
    return str(CXSMILES)
