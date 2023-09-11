import app.modules.toolkits.cdk_wrapper as cdk


def get_sugar_info(molecule: any) -> tuple:
    """
    Analyzes a molecule represented by a SMILES string to determine if it contains sugars.
    This function utilizes the Sugar Removal Utility to check for the presence of circular or linear sugars.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        tuple: A tuple containing two boolean values indicating whether the molecule has linear sugars
               and whether the molecule has circular sugars. If no sugars are found, both values will be False.
    """
    SCOB = cdk.JClass(cdk.cdk_base + ".silent.SilentChemObjectBuilder")

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


def remove_linear_sugar(molecule: any) -> str:
    """
    Detects and removes linear sugars from a given SMILES string using the CDK-based
    sugar removal utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        str: The SMILES string with linear sugars removed, or a message indicating no linear sugar found.

    Raises:
        ValueError: If there is an issue with parsing the input SMILES string.
    """

    cdk_base = "org.openscience.cdk"
    SCOB = cdk.JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmiFlavor = cdk.JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = cdk.JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute
    )

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


def remove_circular_sugar(molecule: any) -> str:
    """
    Detects and removes circular sugars from a given SMILES string using the CDK-based sugar removal utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.

    Returns:
        str: SMILES string with circular sugars removed, or a message if no circular sugars are found.
    """
    cdk_base = "org.openscience.cdk"
    SCOB = cdk.JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmiFlavor = cdk.JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = cdk.JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute
    )

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


def remove_linear_and_circular_sugar(molecule: any):
    """This fucntion detects and removes linear and circular sugars from a give
    SMILES string. Uses the CDK based sugar removal utility.

    Args:
        molecule (IAtomContainer): CDK molecule object.
    Returns:
        smiles (str): SMILES string without linear and circular sugars.

    """
    cdk_base = "org.openscience.cdk"
    SCOB = cdk.JClass(cdk_base + ".silent.SilentChemObjectBuilder")
    SmiFlavor = cdk.JClass(cdk_base + ".smiles.SmiFlavor")
    SmilesGenerator = cdk.JClass(cdk_base + ".smiles.SmilesGenerator")(
        SmiFlavor.Absolute
    )

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
        if MoleculeWithoutSugars.isEmpty():
            raise Exception(
                "Molecule is empty after removal of circular and linear sugars"
            )
        else:
            try:
                S_SMILES = SmilesGenerator.create(MoleculeWithoutSugars)
                return str(S_SMILES)
            except Exception as e:
                raise Exception(f"{str(e)}")
    else:
        return "No Linear or Circular sugars found"
