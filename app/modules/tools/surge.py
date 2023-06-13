from subprocess import Popen, PIPE


def generateStructures(molecular_formula: str):
    """This function uses surge - chemical structure generator that generates
    structures based on the canonical generation path method
    Args:
        molecular_formula (string): molecular_formula string given by the user.
    Returns:
        (array): array of SMILEs generated for the given MF
    """
    smiles = []
    process = Popen(
        ["surge", "-P", "-T", "-B1,2,3,4,5,7,9", "-t0", "-f0", "-S", molecular_formula],
        stdout=PIPE,
        stderr=PIPE,
    )
    for line in iter(process.stdout.readline, b""):
        smiles.append(line.decode("utf-8").rstrip())
    return smiles
