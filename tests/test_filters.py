import pytest
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import (
    get_sas_score,
    get_PAINS,
    get_GhoseFilter,
    get_VeberFilter,
    get_REOSFilter,
    get_RuleofThree,
)


# Define fixtures for example molecules
@pytest.fixture
def molecule1():
    return parse_input("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "rdkit", False)


@pytest.fixture
def molecule2():
    return parse_input("CCCCCCCCCCCCCC", "rdkit", False)


@pytest.fixture
def molecule3():
    return parse_input("C", "rdkit", False)


# Test get_sas_score
def test_get_sas_score(molecule1):
    assert isinstance(get_sas_score(molecule1), float)


# Test get_PAINS
def test_get_PAINS(molecule2):
    assert get_PAINS(molecule2) is False


# Test get_GhoseFilter
def test_get_GhoseFilter(molecule1, molecule2):
    assert get_GhoseFilter(molecule1) is True
    assert get_GhoseFilter(molecule2) is False


# Test get_VeberFilter
def test_get_VeberFilter(molecule1, molecule2):
    assert get_VeberFilter(molecule1) is True
    assert get_VeberFilter(molecule2) is False


# Test get_REOSFilter
def test_get_REOSFilter(molecule1, molecule2):
    assert get_REOSFilter(molecule1) is True
    assert get_REOSFilter(molecule2) is False


# Test get_RuleofThree
def test_get_RuleofThree(molecule3, molecule2):
    assert get_RuleofThree(molecule3) is True
    assert get_RuleofThree(molecule2) is False


if __name__ == "__main__":
    pytest.main()
