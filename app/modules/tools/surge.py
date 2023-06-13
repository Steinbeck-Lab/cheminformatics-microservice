from subprocess import Popen, PIPE


def generateStructures(mf: str):
    """This function uses the sugar removal utility and checks
    whether a molecule has ring or linear sugars
    Args:
        smiles (string): SMILES string given by the user.
    Returns:
        (boolean): True or false values whtehr or not molecule has sugar.
    """
    smiles = []
    process = Popen(
        ["surge", "-P", "-T", "-B1,2,3,4,5,7,9", "-t0", "-f0", "-S", mf],
        stdout=PIPE,
        stderr=PIPE,
    )
    for line in iter(process.stdout.readline, b""):
        smiles.append(line.decode("utf-8").rstrip())
    return smiles
