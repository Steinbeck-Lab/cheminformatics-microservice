import app.modules.toolkits.cdk_wrapper as cdk


def getSugarInfo(smiles: str):
    """This function uses the sugar removal utility and checks
    whether a molecule has circular or linear sugars

    Args:
        smiles (str): SMILES string given by the user.
    Returns:
        (boolean): True or false values whtehr or not molecule has sugar.

    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    SCOB = cdk.JClass(cdk.cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = cdk.JClass(cdk.cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    molecule = SmilesParser.parseSmiles(smiles)

    sru_base = "de.unijena.cheminf.deglycosylation"

    SugarRemovalUtility = cdk.JClass(sru_base + ".SugarRemovalUtility")(
        SCOB.getInstance()
    )
    hasCircularOrLinearSugars = SugarRemovalUtility.hasCircularOrLinearSugars(molecule)

    if hasCircularOrLinearSugars:
        hasLinearSugar = SugarRemovalUtility.hasLinearSugars(molecule)
        hasCircularSugars = SugarRemovalUtility.hasCircularSugars(molecule)
        return hasLinearSugar, hasCircularSugars
    else:
        return (False, False)


def removeLinearSugar(smiles: str):
    """This fucntion detects and removes linear sugars from a give
    SMILES string. Uses the CDK based sugar removal utility.

    Args:
        smiles (str): SMILES string given by the user.
    Returns:
        smiles (str): SMILES string without linear sugars.

    """
    cdk_base = "org.openscience.cdk"
    SCOB = cdk.JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = cdk.JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    SmiFlavor = cdk.JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = cdk.JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute
    )

    molecule = SmilesParser.parseSmiles(smiles)

    sru_base = "de.unijena.cheminf.deglycosylation"

    SugarRemovalUtility = cdk.JClass(sru_base + ".SugarRemovalUtility")(
        SCOB.getInstance()
    )
    hasLinearSugar = SugarRemovalUtility.hasLinearSugars(molecule)

    if hasLinearSugar:
        MoleculeWithoutSugars = SugarRemovalUtility.removeLinearSugars(molecule, True)
        L_SMILES = SmilesGenerator.create(MoleculeWithoutSugars)
        return str(L_SMILES)
    else:
        return "No Linear sugar found"


def removeCircularSugar(smiles: str):
    """This fucntion detects and removes circular sugars from a give
    SMILES string. Uses the CDK based sugar removal utility.

    Args:
        smiles (str): SMILES string given by the user.
    Returns:
        smiles (str): SMILES string without circular sugars.

    """
    cdk_base = "org.openscience.cdk"
    SCOB = cdk.JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = cdk.JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    SmiFlavor = cdk.JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = cdk.JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute
    )

    molecule = SmilesParser.parseSmiles(smiles)

    sru_base = "de.unijena.cheminf.deglycosylation"

    SugarRemovalUtility = cdk.JClass(sru_base + ".SugarRemovalUtility")(
        SCOB.getInstance()
    )
    hasCircularSugar = SugarRemovalUtility.hasCircularSugars(molecule)

    if hasCircularSugar:
        SugarRemovalUtility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(True)
        MoleculeWithoutSugars = SugarRemovalUtility.removeCircularSugars(molecule, True)
        C_SMILES = SmilesGenerator.create(MoleculeWithoutSugars)
        return str(C_SMILES)
    else:
        return "No Circular sugars found"


def removeLinearandCircularSugar(smiles: str):
    """This fucntion detects and removes linear and circular sugars from a give
    SMILES string. Uses the CDK based sugar removal utility.

    Args:
        smiles (str): SMILES string given by the user.
    Returns:
        smiles (str): SMILES string without linear and circular sugars.

    """
    cdk_base = "org.openscience.cdk"
    SCOB = cdk.JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = cdk.JClass(cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    SmiFlavor = cdk.JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = cdk.JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute
    )

    molecule = SmilesParser.parseSmiles(smiles)

    sru_base = "de.unijena.cheminf.deglycosylation"

    SugarRemovalUtility = cdk.JClass(sru_base + ".SugarRemovalUtility")(
        SCOB.getInstance()
    )
    hasCircularOrLinearSugars = SugarRemovalUtility.hasCircularOrLinearSugars(molecule)

    if hasCircularOrLinearSugars:
        SugarRemovalUtility.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(True)
        MoleculeWithoutSugars = SugarRemovalUtility.removeCircularAndLinearSugars(
            molecule, True
        )
        try:
            S_SMILES = SmilesGenerator.create(MoleculeWithoutSugars)
        except Exception as e:
            print(e)
            return "Error generating SMILES"
        return str(S_SMILES)
    else:
        return "No Linear or Circular sugars found"
