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
from app.modules.toolkits.cdk_wrapper import JVMNotFoundException
from app.modules.toolkits.cdk_wrapper import setup_jvm
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import check_RO5_violations
from app.modules.toolkits.rdkit_wrapper import get_3d_conformers
from app.modules.toolkits.rdkit_wrapper import get_ertl_functional_groups_ifg
from app.modules.toolkits.rdkit_wrapper import get_tanimoto_similarity_rdkit
from app.modules.toolkits.rdkit_wrapper import has_cis_trans_stereochemistry


@pytest.fixture
def test_smiles():
    return "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"


@pytest.fixture
def test_smiles_descriptors():
    return "CCC"


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

# Test molecules for has_cis_trans_stereochemistry
mol_with_cis_trans = Chem.MolFromSmiles("C/C=C/C")  # E-isomer
mol_without_cis_trans = Chem.MolFromSmiles("CCC")  # No double bond


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
    svg = get_rdkit_depiction(test_RDKit_Mol, mol_size=(512, 512))
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


def test_all_rdkit_descriptors(test_smiles_descriptors):
    mol = parse_input(test_smiles_descriptors, "rdkit", False)
    descriptors = get_all_rdkit_descriptors(mol)

    # Check specific values that should be consistent across platforms
    assert descriptors[0] == 11, f"AtomC expected 11, got {descriptors[0]}"
    assert descriptors[1] == 2, f"HeavyAtomsC expected 2, got {descriptors[1]}"
    assert descriptors[2] == 3, f"First part expected 3, got {descriptors[2]}"
    assert descriptors[3] == 44.1, f"MolWt expected 44.1, got {descriptors[3]}"
    assert (
        descriptors[4] == 44.0626
    ), f"ExactMolWt expected 44.0626, got {descriptors[4]}"
    assert descriptors[5] == 1.42, f"ALogP expected 1.42, got {descriptors[5]}"
    assert descriptors[6] == 0, f"NumRotatableBonds expected 0, got {descriptors[6]}"

    # PSA can be 0.0 or 0 depending on platform
    assert descriptors[7] in [0, 0.0], f"PSA expected 0 or 0.0, got {descriptors[7]}"

    assert descriptors[8] == 0, f"HBA expected 0, got {descriptors[8]}"
    assert descriptors[9] == 0, f"HBD expected 0, got {descriptors[9]}"
    assert descriptors[10] == 0, f"Lipinski_HBA expected 0, got {descriptors[10]}"
    assert descriptors[11] == 0, f"Lipinski_HBD expected 0, got {descriptors[11]}"
    assert descriptors[12] == 0, f"Ro5Violations expected 0, got {descriptors[12]}"
    assert descriptors[13] == 0, f"AromaticRings expected 0, got {descriptors[13]}"
    assert descriptors[14] == 0.39, f"QEDWeighted expected 0.39, got {descriptors[14]}"
    assert descriptors[15] == 0, f"FormalCharge expected 0, got {descriptors[15]}"
    assert descriptors[16] == 1.0, f"fsp3 expected 1.0, got {descriptors[16]}"
    assert descriptors[17] == 0, f"NumRings expected 0, got {descriptors[17]}"

    # Check VABCVolume with tolerance for platform differences
    volume = descriptors[-1]
    assert isinstance(volume, float), f"Volume should be float, got {type(volume)}"
    assert (
        60.0 <= volume <= 65.0
    ), f"Volume {volume} outside expected range [60.0, 65.0]"


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


def test_all_combined_descriptors(test_smiles_descriptors):
    descriptors = get_cdk_rdkit_combined_descriptors(test_smiles_descriptors)
    expected_result = {
        "Atom count": (11, 11),
        "Bond count": (2, 2),
        "Heavy atom count": (3, 3),
        "Molecular weight": (44.1, 44.1),
        "Exact molecular weight": (44.0626, 44.0626),
        "Calculated LogP": (1.42, 1.74),
        "Rotatable bond count": (0, 0),
        "Topological polar surface area": (0.0, 0.0),
        "Hydrogen bond acceptors": (0, 0),
        "Hydrogen bond donors": (0, 0),
        "Hydrogen bond acceptors (Lipinski)": (0, 0),
        "Hydrogen bond donors (Lipinski)": (0, 0),
        "Lipinski's rule of five violations": (0, 0),
        "Aromatic rings count": (0, 0),
        "QED drug likeliness": (0.39, "None"),
        "Formal Charge": (0, 0),
        "FractionCSP3": (1.0, 1.0),
        "Number of Minimal Rings": (0, 0),
    }

    # Check Van der Waals Volume separately with tolerance for platform differences
    vdw_volume = descriptors.pop("Van der Waals Volume")
    cdk_volume, rdkit_volume = vdw_volume

    # CDK volume should be in range [60.0, 65.0]
    assert isinstance(
        cdk_volume, (int, float)
    ), f"CDK volume should be numeric, got {type(cdk_volume)}"
    assert (
        60.0 <= cdk_volume <= 65.0
    ), f"CDK volume {cdk_volume} outside expected range [60.0, 65.0]"

    # RDKit volume should be approximately 60.45
    assert isinstance(
        rdkit_volume, float
    ), f"RDKit volume should be float, got {type(rdkit_volume)}"
    assert (
        59.0 <= rdkit_volume <= 65.0
    ), f"RDKit volume {rdkit_volume} outside expected range [60.0, 65.0]"

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

    result = get_ertl_functional_groups_ifg(mol)

    assert isinstance(result, list)
    assert len(result) > 0
    # Check for new structured format
    first_group = result[0]
    assert isinstance(first_group, dict)
    assert "atomIds" in first_group
    assert "atoms" in first_group
    assert "type" in first_group
    assert "description" in first_group


def test_get_ertl_functional_groups_no_fragments():
    mol = parse_input("CC", framework="rdkit")
    result = get_ertl_functional_groups_ifg(mol)

    assert isinstance(result, list)
    assert len(result) == 1
    assert result[0] == {"None": "No fragments found"}


def test_setup_jvm_exception(monkeypatch, capsys):
    def mock_get_default_jvm_path():
        raise JVMNotFoundException

    monkeypatch.setattr(
        "app.modules.toolkits.cdk_wrapper.getDefaultJVMPath", mock_get_default_jvm_path
    )

    setup_jvm()
    captured = capsys.readouterr()
    assert "If you see this message" in captured.out
    assert (
        "This indicates that the environment variable JAVA_HOME is not set properly"
        in captured.out
    )
    assert "You can set it or set it manually in the code" in captured.out


# =============================================
# has_cis_trans_stereochemistry Function Tests
# =============================================


def test_has_cis_trans_stereochemistry_with_stereo():
    """Test has_cis_trans_stereochemistry with E/Z stereochemistry."""
    result = has_cis_trans_stereochemistry(mol_with_cis_trans)
    assert isinstance(result, bool)
    assert result is True


def test_has_cis_trans_stereochemistry_without_stereo():
    """Test has_cis_trans_stereochemistry without double bonds."""
    result = has_cis_trans_stereochemistry(mol_without_cis_trans)
    assert isinstance(result, bool)
    assert result is False


def test_has_cis_trans_stereochemistry_none_molecule():
    """Test has_cis_trans_stereochemistry with None molecule."""
    result = has_cis_trans_stereochemistry(None)
    assert isinstance(result, bool)
    assert result is False
