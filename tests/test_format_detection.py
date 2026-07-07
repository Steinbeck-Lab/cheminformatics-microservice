import pytest
from rdkit import Chem

from app.exception_handlers import InvalidInputException
from app.modules.toolkits.helpers import detect_input_format
from app.modules.toolkits.helpers import detect_input_format_with_confidence
from app.modules.toolkits.helpers import input_to_smiles
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.helpers import parse_structure_query

CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
CAFFEINE_INCHI = "InChI=1S/C8H10N4O2/c1-3(9-6)10-5-4(8(13)12)11(2)7(5)14/h3H,1-2H3,(H,9,10,12,13,14)"
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
