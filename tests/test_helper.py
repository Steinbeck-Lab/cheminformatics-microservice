import pytest

from app.modules.toolkits.helpers import parse_input, InvalidInputException


@pytest.fixture
def test_smiles():
    return "C1=CC=CC=C1"


def test_valid_smiles_default_framework(test_smiles):
    mol = parse_input(test_smiles)
    assert mol is not None
    assert mol.GetNumAtoms() == 6


def test_valid_smiles_rdkit_standardize(test_smiles):
    mol = parse_input(test_smiles, framework="rdkit", standardize=True)
    assert mol is not None
    assert mol.GetNumAtoms() == 6


def test_valid_smiles_cdk(test_smiles):
    mol = parse_input(test_smiles, framework="cdk")
    assert mol is not None


def test_valid_smiles_openbabel(test_smiles):
    mol = parse_input(test_smiles, framework="openbabel")
    assert mol is not None


def test_invalid_smiles():
    smiles = "invalid_smiles"
    with pytest.raises(InvalidInputException):
        parse_input(smiles)


def test_parse_SMILES_invalid_framework(test_smiles):
    with pytest.raises(InvalidInputException):
        parse_input(test_smiles, framework="framework")


def test_smiles_with_r_groups():
    smiles = "C1=CC=C[R]C=C1"
    mol = parse_input(smiles, standardize=False)
    assert mol is not None
    assert mol.GetNumAtoms() > 0


def test_smiles_with_spaces():
    smiles = "[NH4 ]"
    mol = parse_input(smiles)
    assert mol is not None
    assert mol.GetNumAtoms() == 1
