from __future__ import annotations

import re

import selfies as sf
from chembl_structure_pipeline import standardizer
from rdkit import Chem

from app.exception_handlers import InvalidInputException
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer_from_molblock
from app.modules.toolkits.cdk_wrapper import get_CDK_SDG_mol
from app.modules.toolkits.cdk_wrapper import get_canonical_SMILES
from app.modules.toolkits.cdk_wrapper import get_smiles_opsin
from app.modules.toolkits.openbabel_wrapper import get_ob_mol
from app.modules.toolkits.rdkit_wrapper import convert_xyz_to_mol

STRUCTURE_INPUT_FORMATS = frozenset(
    {"smiles", "inchi", "selfies", "iupac", "molblock", "xyz"}
)
_SUPPORTED_INPUT_FORMATS = STRUCTURE_INPUT_FORMATS | {"auto"}

_SELFIES_RE = re.compile(r"^(\[[\w@=#+\-/\\%;.,^*]+\])+$")
_MOLBLOCK_COUNTS_RE = re.compile(r"\n\d+\s+\d+\s")


def _normalize_format_name(fmt: str) -> str:
    return fmt.strip().lower()


def _clean_molblock(data: str) -> str:
    if "$$$$" in data:
        return data.split("$$$$")[0].strip()
    return data.strip()


def _looks_like_inchi(text: str) -> bool:
    return text.strip().startswith("InChI=")


def _looks_like_selfies(text: str) -> bool:
    return bool(_SELFIES_RE.match(text.strip()))


def _looks_like_molblock(text: str) -> bool:
    return "M  END" in text or bool(_MOLBLOCK_COUNTS_RE.search(text))


def _looks_like_xyz(text: str) -> bool:
    return bool(split_xyz_frames(text))


def _looks_like_iupac(text: str) -> bool:
    trimmed = text.strip()
    if not trimmed:
        return False
    return bool(
        re.match(r"^[a-z0-9][a-z0-9 ,\-().]*$", trimmed)
        and not re.search(r"[a-z]\d", trimmed)
        and not re.search(r"\d[a-z]", trimmed)
    )


def _can_parse_as_smiles(text: str) -> bool:
    normalized = text.replace(" ", "+")
    if Chem.MolFromSmiles(normalized):
        return True
    try:
        return get_CDK_IAtomContainer(normalized) is not None
    except Exception:
        return False


def _detect_with_confidence(text: str) -> tuple[str, str]:
    """Return ``(detected_format, confidence)`` for a structure string."""
    if not text or not text.strip():
        raise InvalidInputException(name="input", value=text or "")

    trimmed = text.strip()

    if _looks_like_inchi(trimmed):
        return "inchi", "high"
    if _looks_like_selfies(trimmed):
        return "selfies", "high"
    if _looks_like_molblock(text):
        return "molblock", "high"
    if _looks_like_xyz(text):
        return "xyz", "high"
    if _can_parse_as_smiles(trimmed):
        return "smiles", "medium"
    if _looks_like_iupac(trimmed):
        return "iupac", "low"
    return "iupac", "low"


def detect_input_format(text: str) -> str:
    """Detect the chemical structure format of *text*.

    Detection order: InChI, SELFIES, molblock, XYZ, SMILES (try-parse), IUPAC.
    """
    detected_format, _ = _detect_with_confidence(text)
    return detected_format


def detect_input_format_with_confidence(text: str) -> tuple[str, str]:
    """Return ``(detected_format, confidence)`` where confidence is high/medium/low."""
    return _detect_with_confidence(text)


def input_to_smiles(value: str, input_format: str) -> str:
    """Convert structure *value* in *input_format* to a SMILES string."""
    fmt = _normalize_format_name(input_format)
    if fmt == "auto":
        fmt = detect_input_format(value)

    if fmt == "smiles":
        return value
    if fmt == "iupac":
        smiles = get_smiles_opsin(value)
        if not smiles:
            raise ValueError(f"Failed to convert IUPAC name '{value}' to SMILES")
        return smiles
    if fmt == "selfies":
        smiles = sf.decoder(value)
        if not smiles:
            raise ValueError(f"Failed to decode SELFIES '{value}' to SMILES")
        return smiles
    if fmt == "inchi":
        mol = Chem.inchi.MolFromInchi(value)
        if not mol:
            raise ValueError(f"Failed to convert InChI '{value}' to molecule")
        return Chem.MolToSmiles(mol)
    if fmt == "molblock":
        cleaned = _clean_molblock(value)
        mol = Chem.MolFromMolBlock(cleaned)
        if mol:
            return Chem.MolToSmiles(mol)
        mol_container = get_CDK_IAtomContainer_from_molblock(cleaned)
        if mol_container is None:
            raise ValueError("Failed to parse molblock")
        smiles = get_canonical_SMILES(mol_container)
        if not smiles:
            raise ValueError("Failed to parse molblock")
        return str(smiles)
    if fmt == "xyz":
        frames = split_xyz_frames(value)
        if not frames:
            raise ValueError("Failed to parse XYZ")
        mol, _ = convert_xyz_to_mol(frames[0])
        return Chem.MolToSmiles(mol)

    raise ValueError(f"Unsupported input format: {input_format}")


def resolve_input_smiles(value: str, auto_detect: bool = False) -> str:
    """Return SMILES for *value*, auto-detecting format when requested."""
    if auto_detect:
        return input_to_smiles(value, "auto")
    return value


def parse_structure_query(
    smiles: str,
    framework: str = "rdkit",
    auto_detect: bool = False,
    standardize: bool = False,
):
    """Parse a structure query parameter, optionally auto-detecting its format."""
    input_format = "auto" if auto_detect else "smiles"
    return parse_input(smiles, framework, standardize, input_format=input_format)


def parse_input(
    input: str,
    framework: str = "rdkit",
    standardize: bool = False,
    input_format: str = "smiles",
):
    """Parse and check if the input is valid.

    Args:
        input (str): Input string.
        framework (str): Toolkit framework (rdkit, cdk, openbabel).
        standardize (bool): Whether to standardize via ChEMBL pipeline.
        input_format (str): Input format or ``auto`` to detect.

    Returns:
        Mol or None: Valid molecule object or None if an error occurs.
    """
    fmt = _normalize_format_name(input_format)
    if fmt == "auto":
        fmt = detect_input_format(input)
    elif fmt not in _SUPPORTED_INPUT_FORMATS - {"auto"}:
        raise ValueError(f"Unsupported input format: {input_format}")

    if fmt == "smiles":
        return parse_SMILES(input, framework, standardize)

    smiles = input_to_smiles(input, fmt)
    return parse_SMILES(smiles, framework, standardize)


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
