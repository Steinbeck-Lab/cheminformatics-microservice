import pytest
from rdkit import Chem

from app.exception_handlers import InvalidInputException
from app.modules.toolkits.helpers import _is_plausible_atom_row
from app.modules.toolkits.helpers import _looks_like_iupac
from app.modules.toolkits.helpers import _native_rdkit_mol
from app.modules.toolkits.helpers import detect_input_format
from app.modules.toolkits.helpers import detect_input_format_with_confidence
from app.modules.toolkits.helpers import input_to_smiles
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.helpers import parse_structure_query
from app.modules.toolkits.helpers import resolve_input_smiles
from app.modules.toolkits.helpers import split_xyz_frames

CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
CAFFEINE_INCHI = "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"
CAFFEINE_SELFIES = "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]"
ETHANOL_SELFIES = "[C][C][O]"
WATER_XYZ = (
    "3\n"
    "water\n"
    "O    0.0000    0.0000    0.0000\n"
    "H    0.7572    0.5860    0.0000\n"
    "H   -0.7572    0.5860    0.0000\n"
)

MOLBLOCK = """
  CDK     08302311362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.9743    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    1.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    2.8125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6238    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  END
"""


class TestDetectInputFormat:
    def test_detect_inchi(self):
        assert detect_input_format(CAFFEINE_INCHI) == "inchi"

    def test_detect_selfies(self):
        assert detect_input_format(CAFFEINE_SELFIES) == "selfies"

    def test_detect_molblock(self):
        assert detect_input_format(MOLBLOCK) == "molblock"

    def test_detect_xyz(self):
        assert detect_input_format(WATER_XYZ) == "xyz"

    def test_detect_smiles(self):
        assert detect_input_format(CAFFEINE_SMILES) == "smiles"

    def test_detect_iupac_name(self):
        assert detect_input_format("ethanol") == "iupac"

    def test_detect_empty_raises(self):
        with pytest.raises(InvalidInputException):
            detect_input_format("")

    def test_confidence_high_for_inchi(self):
        fmt, confidence = detect_input_format_with_confidence(CAFFEINE_INCHI)
        assert fmt == "inchi"
        assert confidence == "high"

    def test_confidence_low_for_iupac(self):
        fmt, confidence = detect_input_format_with_confidence("ethanol")
        assert fmt == "iupac"
        assert confidence == "low"


class TestParseInputAutoDetect:
    def test_default_smiles_behavior_unchanged(self):
        mol = parse_input(CAFFEINE_SMILES)
        assert mol is not None
        assert mol.GetNumAtoms() == 14

    def test_explicit_smiles_rejects_inchi(self):
        with pytest.raises(InvalidInputException):
            parse_input(CAFFEINE_INCHI, input_format="smiles")

    def test_auto_detect_inchi(self):
        mol = parse_input(CAFFEINE_INCHI, input_format="auto")
        assert mol is not None
        assert mol.GetNumAtoms() == 14

    def test_auto_detect_selfies(self):
        mol = parse_input(ETHANOL_SELFIES, input_format="auto")
        assert mol is not None
        assert mol.GetNumAtoms() == 3

    def test_auto_detect_molblock(self):
        mol = parse_input(MOLBLOCK, input_format="auto")
        assert mol is not None
        assert mol.GetNumAtoms() == 4

    def test_auto_detect_xyz(self):
        mol = parse_input(WATER_XYZ, input_format="auto")
        assert mol is not None
        assert mol.GetNumAtoms() == 3

    def test_parse_structure_query_auto_detect_flag(self):
        mol = parse_structure_query(CAFFEINE_INCHI, auto_detect=True)
        assert mol is not None

    def test_input_to_smiles_round_trip_formats(self):
        smiles = input_to_smiles(CAFFEINE_INCHI, "inchi")
        assert Chem.MolFromSmiles(smiles) is not None

        smiles = input_to_smiles(CAFFEINE_SELFIES, "selfies")
        assert Chem.MolFromSmiles(smiles) is not None

        smiles = input_to_smiles(MOLBLOCK, "molblock")
        assert Chem.MolFromSmiles(smiles) is not None

        smiles = input_to_smiles(WATER_XYZ, "xyz")
        assert Chem.MolFromSmiles(smiles) is not None

    def test_detect_fallback_to_iupac_for_unparseable_gibberish(self):
        # Uppercase start fails the IUPAC heuristic too, but detection must
        # still land somewhere: the low-confidence IUPAC default.
        assert detect_input_format("Q9Z!!!") == "iupac"

    def test_parse_input_rejects_unsupported_format(self):
        with pytest.raises(ValueError):
            parse_input(CAFFEINE_SMILES, input_format="not-a-real-format")

    def test_resolve_input_smiles_auto_detect(self):
        smiles = resolve_input_smiles(CAFFEINE_INCHI, auto_detect=True)
        assert Chem.MolFromSmiles(smiles) is not None

    def test_resolve_input_smiles_no_auto_detect_passthrough(self):
        assert resolve_input_smiles(CAFFEINE_SMILES) == CAFFEINE_SMILES

    def test_molblock_sdf_multi_record_delimiter_is_trimmed(self):
        # A single-record SDF file: the real molblock followed by the "$$$$"
        # terminator and trailing junk that must never reach the parser.
        sdf_style = MOLBLOCK + "> <PROPERTY>\nvalue\n\n$$$$\ngarbage after delimiter\n"
        mol = parse_input(sdf_style, input_format="molblock")
        assert mol is not None
        assert mol.GetNumAtoms() == 4

    def test_native_rdkit_mol_standardize_molblock(self):
        mol = parse_input(MOLBLOCK, input_format="molblock", standardize=True)
        assert mol is not None
        assert mol.GetNumAtoms() == 4

    def test_native_rdkit_mol_standardize_inchi(self):
        mol = parse_input(CAFFEINE_INCHI, input_format="inchi", standardize=True)
        assert mol is not None
        assert mol.GetNumAtoms() == 14

    def test_native_rdkit_mol_standardize_xyz(self):
        # ChEMBL standardization normalizes explicit hydrogens away, same as
        # it already does for standardize=True on plain SMILES input.
        mol = parse_input(WATER_XYZ, input_format="xyz", standardize=True)
        assert mol is not None
        assert mol.GetNumAtoms() == 1


class TestNativeRdkitMolFallback:
    def test_xyz_unparseable_text_returns_none(self):
        assert _native_rdkit_mol("not an xyz block", "xyz") is None

    def test_xyz_bond_perception_failure_returns_none(self):
        # An unrecognized element symbol makes both DetermineBonds tiers
        # raise; the caller must see None, not a propagated exception.
        assert _native_rdkit_mol("1\nbad\nXx 0 0 0\n", "xyz") is None

    def test_unsupported_format_returns_none(self):
        assert _native_rdkit_mol(CAFFEINE_SMILES, "smiles") is None


class TestInputToSmilesDirectErrors:
    def test_xyz_unparseable_text_raises(self):
        with pytest.raises(ValueError):
            input_to_smiles("not an xyz block", "xyz")

    def test_unsupported_format_raises(self):
        with pytest.raises(ValueError):
            input_to_smiles(CAFFEINE_SMILES, "not-a-real-format")

    def test_molblock_total_failure_raises(self):
        # Malformed beyond what either RDKit or the CDK fallback can parse.
        garbage_molblock = "\nfake\n\n  x  y  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
        with pytest.raises(ValueError):
            input_to_smiles(garbage_molblock, "molblock")


class TestFormatDetectionHelpers:
    def test_looks_like_iupac_rejects_empty_string(self):
        assert _looks_like_iupac("") is False

    def test_plausible_atom_row_rejects_too_few_tokens(self):
        assert _is_plausible_atom_row("O 0.0") is False

    def test_plausible_atom_row_rejects_non_alpha_element(self):
        assert _is_plausible_atom_row("1 0.0 0.0 0.0") is False

    def test_plausible_atom_row_rejects_non_numeric_coordinate(self):
        assert _is_plausible_atom_row("O x y z") is False

    def test_split_xyz_frames_empty_input(self):
        assert split_xyz_frames("") == []
        assert split_xyz_frames("   ") == []

    def test_split_xyz_frames_skips_zero_atom_frame(self):
        data = "0\nempty frame\n" + WATER_XYZ
        frames = split_xyz_frames(data)
        assert len(frames) == 1

    def test_split_xyz_frames_recovers_after_truncated_frame(self):
        # Declares far more atoms than exist anywhere in the remaining text,
        # so `end > n` triggers the truncated-frame recovery path.
        data = "999\ntruncated\nO 0 0 0\n" + WATER_XYZ
        frames = split_xyz_frames(data)
        assert len(frames) == 1

    def test_split_xyz_frames_recovers_after_invalid_atom_row(self):
        data = "3\nbad row\nnotanatom !!! @@@ ###\nH 0 0 0\nH 0 0 0\n" + WATER_XYZ
        frames = split_xyz_frames(data)
        assert len(frames) == 1
