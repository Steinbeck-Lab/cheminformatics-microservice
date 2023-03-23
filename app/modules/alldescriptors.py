from rdkit.Chem import Descriptors, QED, Lipinski, rdMolDescriptors, rdmolops
from app.modules.rdkitmodules import checkRo5Violations, checkSMILES
from app.modules.cdkmodules import getCDKSDG, JClass, cdk_base


def getAllRDKitDescriptors(smiles: str):
    """Take an input SMILES and generate a selected set of molecular
    descriptors generated using RDKit as a list.
    Args (str): SMILES string
    Returns (list): a list of calculated descriptors
    """
    mol = checkSMILES(smiles)
    if mol:
        AtomC = rdMolDescriptors.CalcNumAtoms(mol)
        BondC = mol.GetNumBonds()
        HeavyAtomsC = rdMolDescriptors.CalcNumHeavyAtoms(mol)
        MolWt = "%.2f" % Descriptors.MolWt(mol)
        ExactMolWt = "%.2f" % Descriptors.ExactMolWt(mol)
        ALogP = "%.2f" % QED.properties(mol).ALOGP
        NumRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        PSA = "%.2f" % rdMolDescriptors.CalcTPSA(mol)
        HBA = Descriptors.NumHAcceptors(mol)
        HBD = Descriptors.NumHDonors(mol)
        Lipinski_HBA = Lipinski.NumHAcceptors(mol)
        Lipinski_HBD = Lipinski.NumHDonors(mol)
        Ro5Violations = checkRo5Violations(mol)
        AromaticRings = rdMolDescriptors.CalcNumAromaticRings(mol)
        QEDWeighted = "%.2f" % QED.qed(mol)
        FormalCharge = "%.2f" % rdmolops.GetFormalCharge(mol)
        fsp3 = "%.3f" % rdMolDescriptors.CalcFractionCSP3(mol)
        NumRings = rdMolDescriptors.CalcNumRings(mol)
        return (
            AtomC,
            BondC,
            HeavyAtomsC,
            MolWt,
            ExactMolWt,
            ALogP,
            NumRotatableBonds,
            PSA,
            HBA,
            HBD,
            Lipinski_HBA,
            Lipinski_HBD,
            Ro5Violations,
            AromaticRings,
            QEDWeighted,
            FormalCharge,
            fsp3,
            NumRings,
        )
    else:
        return "Error reading SMILES string, check again."


def getAllCDKDescriptors(smiles: str):
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
        BondCountDescriptor = JClass(cdk_base+".qsar.descriptors.molecular.BondCountDescriptor")().calculate(Mol).getValue()
        HeavyAtomsC = Mol.getAtomCount()
        WeightDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.WeightDescriptor")()
            .calculate(Mol)
            .getValue().toString()
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
            .getValue().toString()
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
                cdk_base + ".qsar.descriptors.molecular.HBondDonorCountDescriptor"
            )()
            .calculate(Mol)
            .getValue()
        )
        RuleOfFiveDescriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.RuleOfFiveDescriptor")()
            .calculate(Mol)
            .getValue()
        )
        AromaticRings = None
        QEDWeighted = None
        FormalCharge = JClass(
            cdk_base + ".tools.manipulator.AtomContainerManipulator"
        ).getTotalFormalCharge(Mol)
        FractionalCSP3Descriptor = (
            JClass(cdk_base + ".qsar.descriptors.molecular.FractionalCSP3Descriptor")()
            .calculate(Mol)
            .getValue().toString()
        )
        NumRings = JClass(cdk_base + ".graph.Cycles").mcb(Mol).numberOfCycles()

        return (
            str(AtomCountDescriptor),
            str(BondCountDescriptor),
            HeavyAtomsC,
            "{:.2f}".format(float(str(WeightDescriptor))),
            "{:.2f}".format(float(str(TotalExactMass))),
            "{:.2f}".format(float(str(ALogP).split(",")[0])),
            str(NumRotatableBonds),
            "{:.2f}".format(float(str(TPSADescriptor))),
            str(HBondAcceptorCountDescriptor),
            str(HBondDonorCountDescriptor),
            str(HBondAcceptorCountDescriptor),
            str(HBondDonorCountDescriptor),
            str(RuleOfFiveDescriptor),
            str(AromaticRings),
            str(QEDWeighted),
            FormalCharge,
            "{:.2f}".format(float(str(FractionalCSP3Descriptor))),
            NumRings,
        )


def getCDKRDKitcombinedDescriptors(smiles: str):
    """Take an input SMILES and generates a selected set of molecular
    descriptors using CDK and RDKit and returns as a dictionary.
    Args (str): SMILES string.
    Returns (dist): a dictionary of calculated descriptors.
    """
    RDKitDescriptors = getAllRDKitDescriptors(smiles)
    CDKDescriptors = getAllCDKDescriptors(smiles)
    AllDescriptors = {
        "Atom count",
        "Bond count",
        "Heavy atom count",
        "Molecular weight",
        "Exact molecular weight",
        "ALogP",
        "Rotatable bond count",
        "Topological polar surface area",
        "Hydrogen bond acceptors",
        "Hydrogen bond donors",
        "Hydrogen bond acceptors(Lipinski)",
        "Hydrogen bond donors(Lipinski)",
        "Lipinski's rule of five violations",
        "Aromatic rings count",
        "QED drug likeliness",
        "Formal Charge",
        "FractionCSP3",
        "Number of Minimal Rings",
    }

    if len(AllDescriptors) == len(RDKitDescriptors) == len(CDKDescriptors):
        combinedDict = dict(zip(AllDescriptors, zip(RDKitDescriptors, CDKDescriptors)))
        return combinedDict
    else:
        return "Error dictionary lenth invalid"
