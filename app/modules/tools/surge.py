import re
from subprocess import Popen, PIPE


def getHeavyAtomCount(formula: str) -> int:
    """
    Calculate the heavy atom count from a given molecular formula.

    Args:
        formula (str): The molecular formula of the molecule.

    Returns:
        int: The number of heavy atoms in the molecule.
    """

    elements = re.findall(r"[A-Z][a-z]*\d*", formula)
    heavy_atom_count = 0

    for element in elements:
        element_name = re.findall(r"[A-Z][a-z]*", element)[0]
        element_count = re.findall(r"\d+", element)

        if len(element_count) == 0:
            atom_count = 1
        else:
            atom_count = int(element_count[0])

        if element_name != "H":
            heavy_atom_count += atom_count

    return heavy_atom_count


def generateStructures(molecular_formula: str):
    """This function uses surge - chemical structure generator that generates
    structures based on the canonical generation path method
    Args:
        molecular_formula (string): molecular_formula string given by the user.
    Returns:
        (array): array of SMILEs generated for the given MF
    """
    smiles = []
    if getHeavyAtomCount(molecular_formula) <= 10:
        process = Popen(
            [
                "surge",
                "-P",
                "-T",
                "-B1,2,3,4,5,7,9",
                "-t0",
                "-f0",
                "-S",
                molecular_formula,
            ],
            stdout=PIPE,
            stderr=PIPE,
        )
        for line in iter(process.stdout.readline, b""):
            smiles.append(line.decode("utf-8").rstrip())
        return smiles
    else:
        return "The molecular formula contains more heavy atoms than allowed (10 Heavy Atoms max)."
