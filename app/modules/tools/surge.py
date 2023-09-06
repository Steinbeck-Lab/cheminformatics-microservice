import re
from subprocess import Popen, PIPE
from typing import Union


def get_heavy_atom_count(formula: str) -> int:
    """
    Calculate the heavy atom count from a given molecular formula.

    Args:
        formula (str): The molecular formula of the molecule.

    Returns:
        count (int): The number of heavy atoms in the molecule.

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


def generate_structures_SURGE(molecular_formula: str) -> Union[list, str]:
    """
    Generate chemical structures using the surge tool based on the canonical generation path method.

    Args:
        molecular_formula (str): Molecular formula provided by the user.

    Returns:
        list: List of SMILES strings representing generated chemical structures.
            If the molecular formula contains more than 10 heavy atoms, a message
            indicating the limitation is returned instead.
    """

    smiles = []
    if get_heavy_atom_count(molecular_formula) <= 10:
        try:
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
            stdout, stderr = process.communicate()

            if process.returncode == 0:
                output_lines = stdout.decode("utf-8").splitlines()
                smiles = [line.strip() for line in output_lines]
                return smiles
            else:
                return f"Error running surge: {stderr.decode('utf-8')}"
        except Exception as e:
            return f"An error occurred: {str(e)}"
    else:
        return "The molecular formula contains more heavy atoms than allowed (10 Heavy Atoms max)."
