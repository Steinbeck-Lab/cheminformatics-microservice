import pytest
import selfies as sf
from app.modules.depiction import getRDKitDepiction, getCDKDepiction
from app.modules.npscorer import getNPScore
from app.modules.toolkits.helpers import parseInput

@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

@pytest.fixture
def test_RDKit_Mol():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    return parseInput(smiles, "rdkit", False)

@pytest.fixture
def test_CDK_Mol():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    return parseInput(smiles, "cdk", False)


def test_npscore(test_RDKit_Mol):
    expected_result = "-1.09"
    actual_result = getNPScore(test_RDKit_Mol)
    assert expected_result == actual_result


# RDKit Depiction tests
def test_getRDKitDepiction(test_RDKit_Mol):
    svg = getRDKitDepiction(test_RDKit_Mol)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getRDKitDepiction_kekulize(test_RDKit_Mol):
    svg = getRDKitDepiction(test_RDKit_Mol, kekulize=False)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getRDKitDepiction_rotate(test_RDKit_Mol):
    svg = getRDKitDepiction(test_RDKit_Mol, rotate=90)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getRDKitDepiction_size(test_RDKit_Mol):
    svg = getRDKitDepiction(test_RDKit_Mol, molSize=(512, 512))
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


# CDK depiction tests
def test_getCDKDepiction(test_CDK_Mol):
    svg = getCDKDepiction(test_CDK_Mol)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getCDKDepiction_unicolor(test_CDK_Mol):
    svg = getCDKDepiction(test_CDK_Mol, unicolor=True)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getCDKDepiction_rotate(test_CDK_Mol):
    svg = getCDKDepiction(test_CDK_Mol, rotate=90)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getCDKDepiction_size(test_CDK_Mol):
    svg = getCDKDepiction(test_CDK_Mol, molSize=(512, 512))
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_smilestoselfies():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    expected_result = "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]"
    actual_result = sf.encoder(smiles)
    assert expected_result == actual_result


def test_selfiestosmiles():
    selfies = "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]"
    expected_result = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    actual_result = sf.decoder(selfies)
    assert expected_result == actual_result
