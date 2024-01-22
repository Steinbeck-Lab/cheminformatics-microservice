from __future__ import annotations

import pytest
import selfies as sf
from rdkit import Chem

from app.modules.all_descriptors import get_all_cdk_descriptors
from app.modules.all_descriptors import get_all_rdkit_descriptors
from app.modules.all_descriptors import get_cdk_rdkit_combined_descriptors
from app.modules.all_descriptors import get_tanimoto_similarity
from app.modules.depiction import get_cdk_depiction
from app.modules.depiction import get_rdkit_depiction
from app.modules.npscorer import get_np_score
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import check_RO5_violations
from app.modules.toolkits.rdkit_wrapper import get_3d_conformers
from app.modules.toolkits.rdkit_wrapper import get_ertl_functional_groups
from app.modules.toolkits.rdkit_wrapper import get_tanimoto_similarity_rdkit


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def tanimoto_smiles():
    return "CC,CCO,C"


@pytest.fixture
def test_RDKit_Mol(test_smiles):
    return parse_input(test_smiles, "rdkit", False)


@pytest.fixture
def test_CDK_Mol(test_smiles):
    return parse_input(test_smiles, "cdk", False)


# Example valid molecules
mol1 = Chem.MolFromSmiles("CCO")  # Ethanol
mol2 = Chem.MolFromSmiles("CC")  # Ethane

# Example invalid molecules
invalid_mol1 = None
invalid_mol2 = Chem.MolFromSmiles("Invalid_SMILES")

mol_with_violations = Chem.MolFromSmiles(
    "O=C1OC=2C(=C(O)C(=C(O)C2C(=C1)C=3C=CC=CC3)CC=C(C)C)C(=O)C(C)CC",
)
mol_without_violations = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")


def test_npscore(test_RDKit_Mol):
    expected_result = "-1.09"
    actual_result = get_np_score(test_RDKit_Mol)
    assert expected_result == actual_result


# RDKit Depiction tests
def test_get_rdkit_depiction(test_RDKit_Mol):
    svg = get_rdkit_depiction(test_RDKit_Mol)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_get_rdkit_depiction_kekulize(test_RDKit_Mol):
    svg = get_rdkit_depiction(test_RDKit_Mol, kekulize=False)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_get_rdkit_depiction_rotate(test_RDKit_Mol):
    svg = get_rdkit_depiction(test_RDKit_Mol, rotate=90)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_get_rdkit_depiction_size(test_RDKit_Mol):
    svg = get_rdkit_depiction(test_RDKit_Mol, molSize=(512, 512))
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


# CDK depiction tests
def test_get_cdk_depiction(test_CDK_Mol):
    svg = get_cdk_depiction(test_CDK_Mol)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_get_cdk_depiction_unicolor(test_CDK_Mol):
    svg = get_cdk_depiction(test_CDK_Mol, unicolor=True)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_get_cdk_depiction_rotate(test_CDK_Mol):
    svg = get_cdk_depiction(test_CDK_Mol, rotate=90)
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_get_cdk_depiction_size(test_CDK_Mol):
    svg = get_cdk_depiction(test_CDK_Mol, molSize=(512, 512))
    assert isinstance(svg, str)
    assert "svg" in svg
    assert "Error" not in svg


def test_smilestoselfies(test_smiles):
    expected_result = "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]"
    actual_result = sf.encoder(test_smiles)
    assert expected_result == actual_result


def test_selfiestosmiles(test_smiles):
    selfies = "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]"
    expected_result = test_smiles
    actual_result = sf.decoder(selfies)
    assert expected_result == actual_result


def test_all_rdkit_descriptors(test_smiles):
    mol = parse_input(test_smiles, "rdkit", False)
    descriptors = get_all_rdkit_descriptors(mol)
    expected_result = (
        24,
        15,
        14,
        194.19,
        194.08038,
        -1.03,
        0,
        61.82,
        6,
        0,
        6,
        0,
        0,
        2,
        0.54,
        0,
        0.38,
        2,
        None,
    )
    assert expected_result == descriptors


def test_all_cdk_descriptors(test_CDK_Mol):
    descriptors = get_all_cdk_descriptors(test_CDK_Mol)
    expected_result = (
        24,
        15,
        14,
        194.19,
        194.08038,
        -0.96,
        0,
        56.22,
        6,
        0,
        6,
        0,
        0,
        2,
        "None",
        0,
        0.38,
        2,
        162.33307773672266,
    )
    assert expected_result == descriptors


def test_all_combined_descriptors(test_smiles):
    descriptors = get_cdk_rdkit_combined_descriptors(test_smiles)
    expected_result = {
        "Atom count": (24, 24),
        "Bond count": (15, 15),
        "Heavy atom count": (14, 14),
        "Molecular weight": (194.19, 194.19),
        "Exact molecular weight": (194.08038, 194.08038),
        "Calculated LogP": (-1.03, -0.96),
        "Rotatable bond count": (0, 0),
        "Topological polar surface area": (61.82, 56.22),
        "Hydrogen bond acceptors": (6, 6),
        "Hydrogen bond donors": (0, 0),
        "Hydrogen bond acceptors (Lipinski)": (6, 6),
        "Hydrogen bond donors (Lipinski)": (0, 0),
        "Lipinski's rule of five violations": (0, 0),
        "Aromatic rings count": (2, 2),
        "QED drug likeliness": (0.54, "None"),
        "Formal Charge": (0, 0),
        "FractionCSP3": (0.38, 0.38),
        "Number of Minimal Rings": (2, 2),
        "Van der Waals Volume": (None, 162.33307773672266),
    }
    assert expected_result == descriptors


def test_tanimoto_similarity_rdkit(tanimoto_smiles):
    matrix = get_tanimoto_similarity(tanimoto_smiles, toolkit="rdkit")
    assert len(matrix) == 260
    expected_result = "<table><tr><th></th><th>0</th><th>1</th><th>2</th></tr><tr><td>0</td><td>1.0</td><td>0.14285714285714285</td><td>0.0</td></tr><tr><td>1</td><td>0.14285714285714285</td><td>1.0</td><td>0.0</td></tr><tr><td>2</td><td>0.0</td><td>0.0</td><td>1.0</td></tr></table>"
    assert expected_result == matrix


def test_tanimoto_similarity_cdk(tanimoto_smiles):
    matrix = get_tanimoto_similarity(tanimoto_smiles, toolkit="cdk")
    assert len(matrix) == 264
    expected_result = "<table><tr><th></th><th>0</th><th>1</th><th>2</th></tr><tr><td>0</td><td>1.00000</td><td>0.42857</td><td>0.33333</td></tr><tr><td>1</td><td>0.42857</td><td>1.00000</td><td>0.14286</td></tr><tr><td>2</td><td>0.33333</td><td>0.14286</td><td>1.00000</td></tr></table>"
    assert expected_result == matrix


def test_invalid_toolkit(tanimoto_smiles):
    with pytest.raises(ValueError):
        get_tanimoto_similarity(tanimoto_smiles, toolkit="invalid_toolkit")


def test_valid_ecfp_similarity():
    similarity = get_tanimoto_similarity_rdkit(
        mol1,
        mol2,
        fingerprinter="ECFP",
    )
    assert isinstance(similarity, float)
    assert 0.0 <= similarity <= 1.0


def test_valid_rdkit_similarity():
    similarity = get_tanimoto_similarity_rdkit(
        mol1,
        mol2,
        fingerprinter="RDKit",
    )
    assert isinstance(similarity, float)
    assert 0.0 <= similarity <= 1.0


def test_valid_atompairs_similarity():
    similarity = get_tanimoto_similarity_rdkit(
        mol1,
        mol2,
        fingerprinter="Atompairs",
    )
    assert isinstance(similarity, float)
    assert 0.0 <= similarity <= 1.0


def test_valid_maccs_similarity():
    similarity = get_tanimoto_similarity_rdkit(
        mol1,
        mol2,
        fingerprinter="MACCS",
    )
    assert isinstance(similarity, float)
    assert 0.0 <= similarity <= 1.0


def test_invalid_molecule():
    result = get_tanimoto_similarity_rdkit(
        invalid_mol1,
        mol2,
        fingerprinter="ECFP",
    )
    assert isinstance(result, str)
    assert "Check SMILES strings for Errors" in result


def test_unsupported_fingerprinter():
    result = get_tanimoto_similarity_rdkit(
        mol1,
        mol2,
        fingerprinter="InvalidFingerprinter",
    )
    assert isinstance(result, str)
    assert "Unsupported fingerprinter!" in result


def test_check_RO5_violations():
    violations = check_RO5_violations(mol_with_violations)
    assert violations == 1

    violations = check_RO5_violations(mol_without_violations)
    assert violations == 0


def test_get_3d_conformers():
    mol_with_hydrogens = get_3d_conformers(mol_with_violations, depict=False)
    assert mol_with_hydrogens is not None

    mol_without_hydrogens = get_3d_conformers(
        mol_without_violations,
        depict=False,
    )
    assert mol_without_hydrogens is not None

    mol_molblock = get_3d_conformers(mol_with_violations, depict=True)
    assert isinstance(mol_molblock, str)


def test_valid_rdkit_smiles(test_smiles):
    mol = parse_input(test_smiles, framework="rdkit")
    assert isinstance(mol, Chem.Mol)


def test_invalid_rdkit_smiles():
    with pytest.raises(Exception):
        parse_input(invalid_mol1, framework="rdkit")


def test_valid_cdk_smiles(test_smiles):
    mol = parse_input(test_smiles, framework="cdk")
    assert mol is not None


def test_valid_openbabel_smiles(test_smiles):
    mol = parse_input(test_smiles, framework="openbabel")
    assert mol is not None


def test_get_ertl_functional_groups_valid_molecule(test_smiles):
    mol = parse_input(test_smiles, framework="rdkit")

    result = get_ertl_functional_groups(mol)

    assert isinstance(result, list)
    assert len(result) > 0
    assert str(result[0]) == "IFG(atomIds=(1,), atoms='n', type='cn(c)C')"


def test_get_ertl_functional_groups_no_fragments():
    mol = parse_input("CC", framework="rdkit")
    result = get_ertl_functional_groups(mol)

    assert isinstance(result, list)
    assert len(result) == 1
    assert result[0] == {"None": "No fragments found"}
