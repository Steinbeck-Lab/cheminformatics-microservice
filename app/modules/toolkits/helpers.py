from __future__ import annotations

from chembl_structure_pipeline import standardizer
from rdkit import Chem

from app.exception_handlers import InvalidInputException
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from app.modules.toolkits.cdk_wrapper import get_CDK_SDG_mol
from app.modules.toolkits.openbabel_wrapper import get_ob_mol


def parse_input(input: str, framework: str = "rdkit", standardize: bool = False):
    """Parse and check if the input is valid.

    Args:
        input (str): Input string.

    Returns:
        Mol or None: Valid molecule object or None if an error occurs.
            If an error occurs during SMILES parsing, an error message is returned.
    """

    # auto detect the format
    format = "SMILES"

    if format == "SMILES":
        return parse_SMILES(input, framework, standardize)


def parse_SMILES(smiles: str, framework: str = "rdkit", standardize: bool = False):
    """Check whether the input SMILES string is valid.

    If not, attempt to standardize
    the molecule using the ChEMBL standardization pipeline.

    Args:
        smiles (str): Input SMILES string.
        framework (str): Framework to use for parsing. Default is "rdkit".
        standardize (bool): Whether to standardize the molecule. Default is False.

    Returns:
        Chem.Mol or None: Valid molecule object or None if an error occurs.
            If an error occurs during SMILES parsing, an error message is returned.
    """
    try:
        smiles = smiles.replace(" ", "+")
        if framework == "rdkit":
            if "R" in smiles:
                mol = get_CDK_IAtomContainer(smiles)
                mol_str = get_CDK_SDG_mol(mol)
                mol = Chem.MolFromMolBlock(mol_str)
            else:
                mol = Chem.MolFromSmiles(smiles)
            if standardize:
                mol_block = Chem.MolToMolBlock(mol)
                standardized_mol = standardizer.standardize_molblock(mol_block)
                mol = Chem.MolFromMolBlock(standardized_mol)
        elif framework == "cdk":
            mol = get_CDK_IAtomContainer(smiles)
        elif framework == "openbabel":
            mol = get_ob_mol(smiles)
        else:
            raise ValueError(f"Invalid framework specified: {framework}")

        if mol:
            return mol
        else:
            mol = get_CDK_IAtomContainer(smiles)
            mol_str = get_CDK_SDG_mol(mol)
            return Chem.MolFromMolBlock(mol_str)
    except Exception:
        raise InvalidInputException(name="smiles", value=smiles)


def _is_plausible_atom_row(line: str) -> bool:
    """Return True if *line* looks like a valid XYZ atom row.

    A plausible row, after stripping, has at least 4 whitespace-separated
    tokens where the first token starts with an ASCII letter (element symbol)
    and tokens 2-4 are parseable as floats (the x, y, z coordinates).
    """
    parts = line.strip().split()
    if len(parts) < 4:
        return False
    if not parts[0][0].isascii() or not parts[0][0].isalpha():
        return False
    try:
        float(parts[1])
        float(parts[2])
        float(parts[3])
    except ValueError:
        return False
    return True


def split_xyz_frames(xyz_data: str) -> list[str]:
    """Split a multi-frame XYZ string into a list of single-frame XYZ strings.

    XYZ format: each record is ``<atom_count>\\n<comment>\\n<atom_count atom rows>``.
    A multi-record file simply concatenates these blocks.

    The splitter is permissive:

    * Accepts CRLF or LF line endings (output is always LF).
    * Recovers the next valid frame after a truncated or garbage block by
      advancing one line at a time rather than aborting the whole parse.
    * Validates each candidate atom row before accepting a block, so a fake
      ``atom_count`` cannot greedily consume real lines from the next frame.
    * Drops zero-atom frames (degenerate).
    * Normalizes the count line (strips surrounding whitespace).
    * Each emitted frame is itself a valid XYZ string ending in ``\\n``.

    Args:
        xyz_data: Raw text, possibly containing multiple concatenated XYZ records.

    Returns:
        A list of self-contained single-frame XYZ strings. Empty list if the
        input contains no recognizable frames.
    """
    if not xyz_data or not xyz_data.strip():
        return []

    lines = xyz_data.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    frames: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        token = lines[i].strip()
        if not token.isdecimal():
            i += 1
            continue
        atom_count = int(token)
        if atom_count == 0:
            # Skip the count line + the comment line; emit nothing.
            i += 2
            continue
        # Need atom_count + 2 (count + comment + atoms) lines starting at i.
        end = i + 2 + atom_count
        if end > n:
            # Truncated frame: advance by one line and keep scanning.
            i += 1
            continue
        # Validate every candidate atom row before accepting the block.
        atom_rows = lines[i + 2 : end]
        if not all(_is_plausible_atom_row(row) for row in atom_rows):
            i += 1
            continue
        # Normalize the count line; keep all other lines verbatim.
        block_lines = [lines[i].strip()] + lines[i + 1 : end]
        # Each emitted frame ends in exactly one newline.
        frames.append("\n".join(block_lines) + "\n")
        i = end
    return frames
