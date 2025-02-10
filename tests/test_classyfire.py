import pytest
from app.modules.classyfire import classify, result
import asyncio


@pytest.fixture
def valid_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def invalid_smiles():
    return "invalid_smiles"


async def test_valid_classyfire(valid_smiles):
    result_ = await classify(valid_smiles)
    assert result_["query_type"] == "STRUCTURE"
    id_ = result_["id"]

    while True:
        classified = await result(id_)
        if classified["classification_status"] == "Done":
            break
        await asyncio.sleep(2)

    assert classified["classification_status"] == "Done"
    assert classified["entities"][0]["class"]["name"] == "Imidazopyrimidines"


async def test_invalid_classyfire(invalid_smiles):
    result_ = await classify(invalid_smiles)
    assert result_["query_input"] == "invalid_smiles"
    id_ = result_["id"]

    while True:
        classified = await result(id_)
        if classified["classification_status"] == "Done":
            break
        await asyncio.sleep(2)

    assert classified["classification_status"] == "Done"
    assert (
        classified["invalid_entities"][0]["report"][0]
        == "Cannot process the input SMILES string, please check again"
    )
