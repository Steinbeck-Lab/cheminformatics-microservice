from __future__ import annotations

import pytest
from DECIMER import predict_SMILES
from rdkit import Chem


@pytest.fixture
def test_smiles():
    return "CN1C(=O)C2=C(N=CN2C)N(C)C1=O"


def test_imagetosmiles(test_smiles):
    img_path = "tests/caffeine.png"
    expected_result = test_smiles
    actual_result = predict_SMILES(img_path)
    mol = Chem.MolFromSmiles(actual_result)
    if mol:
        actual_result_canonical = Chem.MolToSmiles(
            mol,
            isomericSmiles=True,
            kekuleSmiles=True,
        )
    assert expected_result == actual_result_canonical
