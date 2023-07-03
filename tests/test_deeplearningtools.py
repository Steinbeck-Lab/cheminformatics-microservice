import pytest
from STOUT import translate_forward, translate_reverse
from DECIMER import predict_SMILES
from rdkit import Chem


@pytest.fixture
def test_smiles():
    return "CN1C(=O)C2=C(N=CN2C)N(C)C1=O"


def test_smilestoiupac(test_smiles):
    smiles = test_smiles
    expected_result = "1,3,7-trimethylpurine-2,6-dione"
    actual_result = translate_forward(smiles)
    assert expected_result == actual_result


def test_iupactosmiles(test_smiles):
    iupac_name = "1,3,7-trimethylpurine-2,6-dione"
    expected_result = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
    actual_result = translate_reverse(iupac_name)
    assert expected_result == actual_result


def test_imagetosmiles(test_smiles):
    img_path = "tests/caffeine.png"
    expected_result = test_smiles
    actual_result = predict_SMILES(img_path)
    mol = Chem.MolFromSmiles(actual_result)
    if mol:
        actual_result_canonical = Chem.MolToSmiles(
            mol, isomericSmiles=True, kekuleSmiles=True
        )
    assert expected_result == actual_result_canonical
