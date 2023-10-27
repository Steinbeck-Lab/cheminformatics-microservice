import pytest
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import (
    check_RO5_violations,
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
    return parse_input(
        "O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2", "rdkit", False
    )


@pytest.fixture
def molecule2():
    return parse_input("CCCCCCCCCCCCCC", "rdkit", False)


@pytest.fixture
def molecule3():
    return parse_input("C", "rdkit", False)


# Test get_sas_score
def test_get_sas_score(molecule1):
    assert isinstance(get_sas_score(molecule1), float)


def test_lipinski(molecule1, molecule2):
    # Test for a molecule with no violations
    assert check_RO5_violations(molecule1) == 0

    # Test for a molecule with one violation (MolLogP > 5)
    assert check_RO5_violations(molecule2) == 1


# Test get_PAINS
def test_get_PAINS(molecule1, molecule2):
    assert get_PAINS(molecule1) == ("PAINS filters (family A)", "Hzone_phenol_a(479)")
    assert "PAINS filters" in str(get_PAINS(molecule1))
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
