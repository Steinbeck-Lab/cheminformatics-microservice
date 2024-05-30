import pytest
from app.modules.classyfire import classify, result
import asyncio


@pytest.fixture
def valid_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def invalid_smiles():
    return "invalid_smiles"

"""
def test_valid_classyfire(valid_smiles):
    loop = asyncio.get_event_loop()
    result_ = loop.run_until_complete(classify(valid_smiles))
    assert result_["query_type"] == "STRUCTURE"
    id_ = result_["id"]
    classified = loop.run_until_complete(result(id_))
    assert classified["classification_status"] == "Done"
    assert classified["entities"][0]["class"]["name"] == "Imidazopyrimidines"


def test_invalid_classyfire(invalid_smiles):
    loop = asyncio.get_event_loop()
    result_ = loop.run_until_complete(classify(invalid_smiles))
    assert result_["query_input"] == "invalid_smiles"
    id_ = result_["id"]
    classified = loop.run_until_complete(result(id_))
    assert classified["classification_status"] == "Done"
    assert (
        classified["invalid_entities"][0]["report"][0]
        == "Cannot process the input SMILES string, please check again"
    )
"""
