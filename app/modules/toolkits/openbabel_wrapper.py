from __future__ import annotations

from openbabel import openbabel as ob
from openbabel import pybel

from app.exception_handlers import InvalidInputException


def get_ob_canonical_SMILES(smiles: str) -> str:
    """Convert a SMILES string to Canonical SMILES.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        str: Canonical SMILES string.
    """
    smiles = smiles.replace(" ", "+")

    # Create an Open Babel molecule object
    mol = ob.OBMol()

    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "can")
    conv.ReadString(mol, smiles)

    if mol.NumAtoms() <= 0:
        raise InvalidInputException(name="smiles", value=smiles)
    else:
        canSMILES = conv.WriteString(mol)
        canSMILES = canSMILES.strip()  # Remove leading/trailing whitespace
        return canSMILES


def get_ob_InChI(smiles: str, InChIKey: bool = False) -> str:
    """Convert a SMILES string to InChI.

    Args:
        smiles (str): Input SMILES string.
        InChIKey (bool, optional): Whether to return InChIKey. Defaults to False.

    Returns:
        str: InChI string or InChIKey string if InChIKey is True.
    """

    smiles = smiles.replace(" ", "+")

    # Create an Open Babel molecule object
    mol = ob.OBMol()

    # Create OBConversion
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "inchi")
    conv.ReadString(mol, smiles)

    if mol.NumAtoms() <= 0:
        raise InvalidInputException(name="smiles", value=smiles)
    else:
        inchi = conv.WriteString(mol)
        inchi = inchi.strip()  # Remove leading/trailing whitespace
        if InChIKey:
            conv.SetOptions("K", conv.OUTOPTIONS)
            inchikey_ = conv.WriteString(mol).rstrip()
            return inchikey_
        return inchi


def get_ob_mol(smiles: str, threeD: bool = False, depict: bool = False) -> str:
    """Convert a SMILES string to a 2D/3D mol block.

    Args:
        smiles (str): Input SMILES string.
        threeD (bool, optional): Generate 3D structure. Defaults to False.
        depict (bool, optional): Generate 3D structure for depiction. Defaults to False.

    Returns:
        str: Mol block (2D/3D).
    """
    smiles = smiles.replace(" ", "+")

    if threeD:
        try:
            mol = pybel.readstring("smi", smiles)
        except OSError as e:
            raise InvalidInputException(name="smiles", value=smiles) from e
        else:
            mol.addh()
            mol.make3D()
            gen3d = ob.OBOp.FindType("gen3D")
            gen3d.Do(mol.OBMol, "--best")
            if depict:
                return mol.write("mol")
            else:
                mol.removeh()
                return mol.write("mol")

    # Create an Open Babel molecule object
    mol = ob.OBMol()

    conv = ob.OBConversion()
    conv.SetInAndOutFormats("smi", "mol")
    conv.ReadString(mol, smiles)

    if mol.NumAtoms() <= 0:
        raise InvalidInputException(name="smiles", value=smiles)
    else:
        # Generate 2D coordinates
        obBuilder = ob.OBBuilder()
        obBuilder.Build(mol)

        mol_block = conv.WriteString(mol)
        mol_block = mol_block.strip()  # Remove leading/trailing whitespace
        return mol_block


def get_ob_xyz_conversions(xyz_data: str) -> dict:
    """Parse an XYZ block once via OpenBabel and emit all standard formats.

    OpenBabel performs distance-based bond perception automatically when
    parsing XYZ. It does not accept an explicit charge parameter, so this
    path is best for neutral species; charged species should prefer the
    RDKit pathway, which uses the xyz2mol algorithm with charge support.

    Args:
        xyz_data (str): Plain-text XYZ block.

    Returns:
        dict: Keys ``canonicalsmiles``, ``inchi``, ``inchikey``, and
            ``molblock`` (V2000 MOL block with the original 3D coordinates).

    Raises:
        InvalidInputException: If the XYZ block cannot be parsed by
            OpenBabel or contains no atoms.
    """
    if not xyz_data or not xyz_data.strip():
        raise InvalidInputException(name="xyz", value=xyz_data)

    try:
        mol = pybel.readstring("xyz", xyz_data)
    except (OSError, IOError) as exc:
        raise InvalidInputException(name="xyz", value=xyz_data) from exc

    if mol.OBMol.NumAtoms() <= 0:
        raise InvalidInputException(name="xyz", value=xyz_data)

    # OpenBabel's "can" writer appends "\t<title>" from the XYZ comment line.
    # Keep only the SMILES token so the response matches the RDKit branch.
    can_raw = mol.write("can").strip()
    canonical = can_raw.split("\t", 1)[0].split()[0] if can_raw else ""

    return {
        "canonicalsmiles": canonical,
        "inchi": mol.write("inchi").strip(),
        "inchikey": mol.write("inchikey").strip(),
        "molblock": mol.write("mol").strip(),
    }
