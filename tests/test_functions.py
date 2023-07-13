import pytest
import selfies as sf
from app.modules.depiction import getRDKitDepiction, getCDKDepiction
from app.modules.npscorer import getNPScore


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


def test_npscore(test_smiles):
    expected_result = "-1.09"
    actual_result = getNPScore(test_smiles)
    assert expected_result == actual_result


# RDKit Depiction tests
def test_getRDKitDepiction(test_smiles):
    svg = getRDKitDepiction(test_smiles)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getRDKitDepiction_kekulize(test_smiles):
    svg = getRDKitDepiction(test_smiles, kekulize=False)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getRDKitDepiction_rotate(test_smiles):
    svg = getRDKitDepiction(test_smiles, rotate=90)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getRDKitDepiction_size(test_smiles):
    svg = getRDKitDepiction(test_smiles, molSize=(512, 512))
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


# CDK depiction tests
def test_getCDKDepiction(test_smiles):
    svg = getCDKDepiction(test_smiles)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getCDKDepiction_unicolor(test_smiles):
    svg = getCDKDepiction(test_smiles, unicolor=True)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getCDKDepiction_rotate(test_smiles):
    svg = getCDKDepiction(test_smiles, rotate=90)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_getCDKDepiction_size(test_smiles):
    svg = getCDKDepiction(test_smiles, molSize=(512, 512))
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
