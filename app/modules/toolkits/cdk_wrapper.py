from __future__ import annotations

import os
from typing import List
from typing import Union

import pystow
from jpype import getDefaultJVMPath
from jpype import isJVMStarted
from jpype import JClass
from jpype import JPackage
from jpype import JVMNotFoundException
from jpype import startJVM


def setup_jvm():
    try:
        jvmPath = getDefaultJVMPath()
    except JVMNotFoundException:
        print("If you see this message, for some reason JPype cannot find jvm.dll.")
        print(
            "This indicates that the environment variable JAVA_HOME is not set properly."
        )
        print("You can set it or set it manually in the code")
        jvmPath = "Define/path/or/set/JAVA_HOME/variable/properly"

    print(jvmPath)

    if not isJVMStarted():
        paths = {
            "cdk-2.10": "https://github.com/cdk/cdk/releases/download/cdk-2.10/cdk-2.10.jar",
            "SugarRemovalUtility-jar-with-dependencies": "https://github.com/JonasSchaub/SugarRemoval/releases/download/v1.3.2/SugarRemovalUtility-jar-with-dependencies.jar",
            "centres": "https://github.com/SiMolecule/centres/releases/download/1.0/centres.jar",
            "opsin-cli-2.8.0-jar-with-dependencies": "https://github.com/dan2097/opsin/releases/download/2.8.0/opsin-cli-2.8.0-jar-with-dependencies.jar",
        }

        jar_paths = {
            key: str(pystow.join("STOUT-V2")) + f"/{key}.jar" for key in paths.keys()
        }
        for key, url in paths.items():
            if not os.path.exists(jar_paths[key]):
                pystow.ensure("STOUT-V2", url=url)

        startJVM("-ea", "-Xmx4096M", classpath=[jar_paths[key] for key in jar_paths])


setup_jvm()
cdk_base = "org.openscience.cdk"
opsin_base = JPackage("uk").ac.cam.ch.wwmm.opsin
_nametostruct = opsin_base.NameToStructure.getInstance()
_restoinchi = opsin_base.NameToInchi.convertResultToInChI


def get_CDK_IAtomContainer(smiles: str):
    """This function takes the input SMILES and creates a CDK IAtomContainer.

    Args:
        smiles (str): SMILES string as input.

    Returns:
        mol (object): IAtomContainer with CDK.
    """
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = JClass(
        cdk_base + ".smiles.SmilesParser",
    )(SCOB.getInstance())
    molecule = SmilesParser.parseSmiles(smiles)
    return molecule


def get_CDK_SDG(molecule: any):
    """This function takes the input IAtomContainer and Creates a.

    Structure Diagram Layout using the CDK.

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        mol object: mol object with CDK SDG.
    """
    StructureDiagramGenerator = JClass(
        cdk_base + ".layout.StructureDiagramGenerator",
    )()
    StructureDiagramGenerator.generateCoordinates(molecule)
    molecule_ = StructureDiagramGenerator.getMolecule()

    return molecule_


def get_CDK_SDG_mol(molecule: any, V3000=False) -> str:
    """Returns a mol block string with Structure Diagram Layout for the given.

    SMILES.

    Args:
        molecule (IAtomContainer): molecule given by the user.
        V3000 (bool, optional): Option to return V3000 mol. Defaults to False.

    Returns:
        str: CDK Structure Diagram Layout mol block.
    """
    StringW = JClass("java.io.StringWriter")()
    moleculeSDG = get_CDK_SDG(molecule)
    SDFW = JClass(cdk_base + ".io.SDFWriter")(StringW)
    SDFW.setAlwaysV3000(V3000)
    SDFW.write(moleculeSDG)
    SDFW.flush()
    mol_str = str(StringW.toString())
    return mol_str


def get_murko_framework(molecule: any) -> str:
    """This function takes the user input SMILES and returns.

    the Murko framework

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        smiles (string): Murko Framework as SMILES.
    """

    MurkoFragmenter = JClass(cdk_base + ".fragment.MurckoFragmenter")(True, 3)
    MurkoFragmenter.generateFragments(molecule)
    if len(MurkoFragmenter.getFrameworks()) == 0:
        return "None"

    return str(MurkoFragmenter.getFrameworks()[0])


def get_aromatic_ring_count(molecule) -> int:
    """Calculate the number of aromatic rings present in a given molecule.

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        int: The number of aromatic rings present in the molecule.
    """

    Cycles = JClass(cdk_base + ".graph.Cycles")
    ElectronDonation = JClass(cdk_base + ".aromaticity.ElectronDonation")

    Aromaticity = JClass(cdk_base + ".aromaticity.Aromaticity")(
        ElectronDonation.daylight(),
        Cycles.cdkAromaticSet(),
    )
    Aromaticity.apply(molecule)
    MCBRings = Cycles.mcb(molecule).toRingSet()
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


def get_vander_waals_volume(molecule: any) -> float:
    """Calculate the Van der Waals volume of a given molecule.

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        float: The Van der Waals volume of the molecule.
    """

    AtomContainerManipulator = JClass(
        cdk_base + ".tools.manipulator.AtomContainerManipulator",
    )
    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule)
    VABCVolume = JClass(
        cdk_base + ".geometry.volume.VABCVolume",
    )().calculate(molecule)
    return VABCVolume


def get_CDK_MolecularFormula(molecule: any) -> str:
    """This function takes the input SMILES and creates a CDK IAtomContainer.

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        str : MolecularFormula generated using CDK.
    """
    MolecularFormulaManipulator = JClass(
        cdk_base + ".tools.manipulator.MolecularFormulaManipulator"
    )

    MolecularFormula = MolecularFormulaManipulator.getMolecularFormula(molecule)
    return MolecularFormulaManipulator.getString(MolecularFormula)


def get_CDK_descriptors(molecule: any) -> Union[tuple, str]:
    """Take an input SMILES and generate a selected set of molecular.

    descriptors generated using CDK as a list.

    Args (str):     molecule (IAtomContainer): molecule given by the
    user.

    Returns (list):     A list of calculated descriptors.
    """
    SDGMol = get_CDK_SDG(molecule)
    if SDGMol:
        AtomCountDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.AtomCountDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        HeavyAtomsC = SDGMol.getAtomCount()
        WeightDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.WeightDescriptor")()
            .calculate(SDGMol)
            .getValue()
            .toString()
        )
        TotalExactMass = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator",
        ).getTotalExactMass(SDGMol)
        ALogP = (
            JClass(cdk_base + ".qsar.descriptors.molecular.ALOGPDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        NumRotatableBonds = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.RotatableBondsCountDescriptor",
            )()
            .calculate(SDGMol)
            .getValue()
        )
        TPSADescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.TPSADescriptor")()
            .calculate(SDGMol)
            .getValue()
            .toString()
        )
        HBondAcceptorCountDescriptor = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
            )()
            .calculate(SDGMol)
            .getValue()
        )
        HBondDonorCountDescriptor = (
            JClass(
                cdk_base + ".qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
            )()
            .calculate(SDGMol)
            .getValue()
        )
        RuleOfFiveDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.RuleOfFiveDescriptor")()
            .calculate(SDGMol)
            .getValue()
        )
        AromaticRings = get_aromatic_ring_count(SDGMol)
        QEDWeighted = None
        FormalCharge = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator",
        ).getTotalFormalCharge(SDGMol)
        FractionalCSP3Descriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.FractionalCSP3Descriptor")()
            .calculate(SDGMol)
            .getValue()
            .toString()
        )
        NumRings = (
            JClass(
                cdk_base + ".graph.Cycles",
            )
            .mcb(SDGMol)
            .numberOfCycles()
        )
        VABCVolume = get_vander_waals_volume(SDGMol)

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
            float("{:.2f}".format(float(str(VABCVolume)))),
        )
    else:
        return "Check input and try again!"


def get_tanimoto_similarity_PubChem_CDK(mol1: any, mol2: any) -> str:
    """Calculate the Tanimoto similarity index between two molecules using.

    PubChem fingerprints.

    Args:
        mol1 (IAtomContainer): First molecule given by the user.
        mol2 (IAtomContainer): Second molecule given by the user.

    Returns:
        str: The Tanimoto similarity as a string with 5 decimal places, or an error message.
    """

    Tanimoto = JClass(cdk_base + ".similarity.Tanimoto")
    SCOB = JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    PubchemFingerprinter = JClass(cdk_base + ".fingerprint.PubchemFingerprinter")(
        SCOB.getInstance(),
    )
    CDKHydrogenAdder = JClass(cdk_base + ".tools.CDKHydrogenAdder").getInstance(
        SCOB.getInstance(),
    )
    AtomContainerManipulator = JClass(
        cdk_base + ".tools.manipulator.AtomContainerManipulator",
    )
    Cycles = JClass(cdk_base + ".graph.Cycles")
    ElectronDonation = JClass(cdk_base + ".aromaticity.ElectronDonation")
    Aromaticity = JClass(cdk_base + ".aromaticity.Aromaticity")(
        ElectronDonation.cdk(),
        Cycles.cdkAromaticSet(),
    )
    if mol1 and mol2:
        # Perceive atom types and configure atoms
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1)
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2)

        # add Implicit Hydrogens
        CDKHydrogenAdder.addImplicitHydrogens(mol1)
        CDKHydrogenAdder.addImplicitHydrogens(mol2)

        # Convert implicit to explicit Hydrogens
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol1)
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol2)

        # Apply Aromaticity
        Aromaticity.apply(mol1)
        Aromaticity.apply(mol2)

        # Generate BitSets using PubChemFingerprinter
        fingerprint1 = PubchemFingerprinter.getBitFingerprint(mol1).asBitSet()
        fingerprint2 = PubchemFingerprinter.getBitFingerprint(mol2).asBitSet()

        # Calculate Tanimoto similarity
        Similarity = Tanimoto.calculate(fingerprint1, fingerprint2)

        return "{:.5f}".format(float(str(Similarity)))
    else:
        return "Check the SMILES string for errors"


def get_tanimoto_similarity_ECFP_CDK(
    mol1: any, mol2: any, ECFP: int = 2, bitset_len: int = 2048
) -> str:
    """Calculate the Tanimoto similarity index between two molecules using.

    CircularFingerprinter fingerprints.

    https://cdk.github.io/cdk/2.8/docs/api/org/openscience/cdk/fingerprint/CircularFingerprinter.html

    Args:
        mol1 (IAtomContainer): First molecule given by the user.
        mol2 (IAtomContainer): Second molecule given by the user.
        ECFP (int): The ECFP version to use (2, 4, or 6) defaults to 2.
        bitset_len (int): The length of the bitset.

    Returns:
        str: The Tanimoto similarity as a string with 5 decimal places, or an error message.
    """
    Tanimoto = JClass(cdk_base + ".similarity.Tanimoto")
    CircularFingerprinter = JClass(
        cdk_base + ".fingerprint.CircularFingerprinter",
    )()
    if ECFP == 2:
        fingerprinter_class = CircularFingerprinter.CLASS_ECFP2
    elif ECFP == 4:
        fingerprinter_class = CircularFingerprinter.CLASS_ECFP4
    elif ECFP == 6:
        fingerprinter_class = CircularFingerprinter.CLASS_ECFP6
    else:
        return "only ECFP 2/4/6 allowed"

    CircularFingerprinter_ECFP = JClass(
        cdk_base + ".fingerprint.CircularFingerprinter",
    )(fingerprinter_class, bitset_len)

    if mol1 and mol2:
        fingerprint1 = CircularFingerprinter_ECFP.getBitFingerprint(mol1)
        fingerprint2 = CircularFingerprinter_ECFP.getBitFingerprint(mol2)
        # Calculate Tanimoto similarity
        Similarity = Tanimoto.calculate(fingerprint1, fingerprint2)
        return "{:.5f}".format(float(str(Similarity)))
    else:
        return "Check the SMILES string for errors"


def get_tanimoto_similarity_CDK(
    mol1: any,
    mol2: any,
    fingerprinter: str = "PubChem",
    ECFP: int = 2,
    bitset_len: int = 2048,
) -> float:
    """Calculate the Tanimoto similarity between two molecules using.

    PubChem/CircularFingerprints in CDK.

    Args:
        mol1 (IAtomContainer): First molecule given by the user.
        mol2 (IAtomContainer): Second molecule given by the user.
        fingerprinter (str, optional): The fingerprinter to use. Currently, only "PubChem/ECFP" is supported. Defaults to "PubChem".
        ECFP (int, optional): The ECFP version to use (2, 4, or 6). Defaults to 2.
        bitset_len (int, optional): The length of the bitset. Defaults to 2048.

    Returns:
        float: The Tanimoto similarity score between the two molecules.

    Raises:
        ValueError: If an unsupported fingerprinter is specified.
    """
    if fingerprinter == "PubChem":
        tanimoto = get_tanimoto_similarity_PubChem_CDK(mol1, mol2)
    elif fingerprinter == "ECFP":
        tanimoto = get_tanimoto_similarity_ECFP_CDK(mol1, mol2, ECFP, bitset_len)
    else:
        raise ValueError(
            "Unsupported fingerprinter. Currently, only 'PubChem' and 'ECFP' is supported.",
        )

    return tanimoto


def get_cip_annotation(molecule: any) -> str:
    """Return the CIP (Cahn–Ingold–Prelog) annotations using the CDK CIP.

    toolkit.

    This function takes a SMILES (Simplified Molecular Input Line Entry System) string
    as input and returns a CIP annotated molecule block using the CDK CIP toolkit.

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        str: A CIP annotated molecule block.
    """
    SDGMol = get_CDK_SDG(molecule)
    centres_base = "com.simolecule.centres"
    Cycles = JClass(cdk_base + ".graph.Cycles")
    IBond = JClass(cdk_base + ".interfaces.IBond")
    IStereoElement = JClass(cdk_base + ".interfaces.IStereoElement")
    Stereocenters = JClass(cdk_base + ".stereo.Stereocenters")
    StandardGenerator = JClass(
        cdk_base + ".renderer.generators.standard.StandardGenerator",
    )

    BaseMol = JClass(centres_base + ".BaseMol")
    CdkLabeller = JClass(centres_base + ".CdkLabeller")
    Descriptor = JClass(centres_base + ".Descriptor")

    stereocenters = Stereocenters.of(SDGMol)
    for atom in SDGMol.atoms():
        if (
            stereocenters.isStereocenter(atom.getIndex())
            and stereocenters.elementType(atom.getIndex())
            == Stereocenters.Type.Tetracoordinate
        ):
            atom.setProperty(StandardGenerator.ANNOTATION_LABEL, "(?)")

    # Iterate over bonds
    for bond in SDGMol.bonds():
        if bond.getOrder() != IBond.Order.DOUBLE:
            continue
        begIdx = bond.getBegin().getIndex()
        endIdx = bond.getEnd().getIndex()
        if (
            stereocenters.elementType(
                begIdx,
            )
            == Stereocenters.Type.Tricoordinate
            and stereocenters.elementType(endIdx) == Stereocenters.Type.Tricoordinate
            and stereocenters.isStereocenter(begIdx)
            and stereocenters.isStereocenter(endIdx)
        ):
            # Check if not in a small ring <7
            if Cycles.smallRingSize(bond, 7) == 0:
                bond.setProperty(StandardGenerator.ANNOTATION_LABEL, "(?)")

    # no defined stereo?
    if not SDGMol.stereoElements().iterator().hasNext():
        return SDGMol

    # Call the Java method
    CdkLabeller.label(SDGMol)

    # Update to label appropriately for racemic and relative stereochemistry
    for se in SDGMol.stereoElements():
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
                            BaseMol.CIP_LABEL_KEY,
                            label.toString() + inv.name(),
                        )
                elif (se.getGroupInfo() & IStereoElement.GRP_REL) != 0:
                    if label == Descriptor.R or label == Descriptor.S:
                        focus.setProperty(
                            BaseMol.CIP_LABEL_KEY,
                            label.toString() + "*",
                        )

    # Iterate over atoms
    for atom in SDGMol.atoms():
        if atom.getProperty(BaseMol.CONF_INDEX) is not None:
            atom.setProperty(
                StandardGenerator.ANNOTATION_LABEL,
                StandardGenerator.ITALIC_DISPLAY_PREFIX
                + atom.getProperty(BaseMol.CONF_INDEX).toString(),
            )
        elif atom.getProperty(BaseMol.CIP_LABEL_KEY) is not None:
            atom.setProperty(
                StandardGenerator.ANNOTATION_LABEL,
                StandardGenerator.ITALIC_DISPLAY_PREFIX
                + atom.getProperty(BaseMol.CIP_LABEL_KEY).toString(),
            )

    # Iterate over bonds
    for bond in SDGMol.bonds():
        if bond.getProperty(BaseMol.CIP_LABEL_KEY) is not None:
            bond.setProperty(
                StandardGenerator.ANNOTATION_LABEL,
                StandardGenerator.ITALIC_DISPLAY_PREFIX
                + bond.getProperty(BaseMol.CIP_LABEL_KEY).toString(),
            )

    return SDGMol


def get_CXSMILES(molecule: any) -> str:
    """Generate CXSMILES representation with 2D atom coordinates from the.

    given.

    SMILES.

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        str: CXSMILES representation with 2D atom coordinates.
    """
    SDGMol = get_CDK_SDG(molecule)
    SmiFlavor = JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute | SmiFlavor.CxSmilesWithCoords,
    )
    CXSMILES = SmilesGenerator.create(SDGMol)
    return str(CXSMILES)


def get_canonical_SMILES(molecule: any) -> str:
    """Generate Canonical SMILES representation with 2D atom coordinates from.

    the given SMILES.

    Args:
        molecule (IAtomContainer): molecule given by the user.

    Returns:
        str: Canonical SMILES representation with 2D atom coordinates.
    """
    SDGMol = get_CDK_SDG(molecule)
    SmiFlavor = JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = JClass(
        cdk_base + ".smiles.SmilesGenerator",
    )(SmiFlavor.Absolute)
    CanonicalSMILES = SmilesGenerator.create(SDGMol)
    return str(CanonicalSMILES)


def get_InChI(molecule: any, InChIKey=False) -> str:
    """Generate InChI or InChIKey from the given SMILES string.

    Args:
        molecule (IAtomContainer): molecule given by the user.
        InChIKey (bool): If True, return InChIKey instead of InChI. The default is False.

    Returns:
        str: InChI or InChIKey string.
    """
    SDGMol = get_CDK_SDG(molecule)
    InChIGeneratorFactory = JClass(cdk_base + ".inchi.InChIGeneratorFactory")
    InChI = InChIGeneratorFactory.getInstance().getInChIGenerator(SDGMol).getInchi()
    if InChIKey:
        InChIKey = (
            InChIGeneratorFactory.getInstance().getInChIGenerator(SDGMol).getInchiKey()
        )
        return InChIKey
    return InChI


def get_smiles_opsin(input_text: str) -> str:
    """Convert IUPAC chemical name to SMILES notation using OPSIN.

    Parameters:
    - input_text (str): The IUPAC chemical name to be converted.

    Returns:
    - str: The SMILES notation corresponding to the given IUPAC name.

    Raises:
    - Exception: If the IUPAC name is not valid or if there are issues in the conversion process. The exception message will guide the user to check the data again.
    """
    try:
        OpsinResult = _nametostruct.parseChemicalName(input_text)
        if str(OpsinResult.getStatus()) == "FAILURE":
            raise Exception(
                (
                    "Failed to convert '%s' to format '%s'\n%s using OPSIN"
                    % (input_text, format, OpsinResult.getMessage())
                ),
            )
        return str(OpsinResult.getSmiles())
    except Exception:
        return str(
            "Failed to convert '%s' to format '%s'\n%s using OPSIN"
            % (input_text, format, OpsinResult.getMessage()),
        )


async def get_CDK_HOSE_codes(
    molecule: any,
    noOfSpheres: int,
    ringsize: bool,
) -> List[str]:
    """Generate CDK-generated HOSECodes for the given SMILES.

    Args:
        molecule (IAtomContainer): molecule given by the user.
        noOfSpheres (int): Number of spheres for HOSECode generation.
        ringsize (bool): Whether to consider ring size for HOSECode generation.

    Returns:
        List[str]: List of CDK-generated HOSECodes.
    """
    HOSECodeGenerator = JClass(cdk_base + ".tools.HOSECodeGenerator")()
    HOSECodes = []
    atoms = molecule.atoms()
    for atom in atoms:
        moleculeHOSECode = HOSECodeGenerator.getHOSECode(
            molecule,
            atom,
            noOfSpheres,
            ringsize,
        )
        HOSECodes.append(str(moleculeHOSECode))
    return HOSECodes
