import pytest
import selfies as sf
from STOUT import translate_forward, translate_reverse
from DECIMER import predict_SMILES


def test_smilestoiupac():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    expected_result = "1,3,7-trimethylpurine-2,6-dione"
    actual_result = translate_forward(smiles)
    assert expected_result == actual_result


def test_iupactosmiles():
    iupac_name = "1,3,7-trimethylpurine-2,6-dione"
    expected_result = "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
    actual_result = translate_reverse(iupac_name)
    assert expected_result == actual_result


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


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_imagetosmiles():
    img_path = "tests/caffeine.png"
    expected_result = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    actual_result = predict_SMILES(img_path)
    assert expected_result == actual_result
