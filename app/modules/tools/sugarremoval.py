import app.modules.toolkits.cdk_wrapper as cdk

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
    SCOB = cdk.JClass(cdk.cdk_base + ".silent.SilentChemObjectBuilder")
    SmilesParser = cdk.JClass(cdk.cdk_base + ".smiles.SmilesParser")(SCOB.getInstance())
    molecule = SmilesParser.parseSmiles(smiles)

    sru_base = "de.unijena.cheminf.deglycosylation"

    SugarRemovalUtility = cdk.JClass(sru_base + ".SugarRemovalUtility")(SCOB.getInstance())
    hasCircularOrLinearSugars = SugarRemovalUtility.hasCircularOrLinearSugars(molecule)

    if hasCircularOrLinearSugars:
        hasLinearSugar = SugarRemovalUtility.hasLinearSugars(molecule)
        hasCircularSugars = SugarRemovalUtility.hasCircularSugars(molecule)
        return hasLinearSugar, hasCircularSugars
    else:
        return (False, False)