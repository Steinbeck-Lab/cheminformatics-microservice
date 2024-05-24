from __future__ import annotations

import pytest
from rdkit import Chem

from app.exception_handlers import InvalidInputException
from app.modules.toolkits.openbabel_wrapper import get_ob_canonical_SMILES
from app.modules.toolkits.openbabel_wrapper import get_ob_InChI
from app.modules.toolkits.openbabel_wrapper import get_ob_mol
from app.modules.toolkits.rdkit_wrapper import check_RO5_violations
from app.modules.toolkits.rdkit_wrapper import get_3d_conformers
from app.modules.toolkits.rdkit_wrapper import has_stereo_defined
from app.modules.toolkits.rdkit_wrapper import is_valid_molecule


@pytest.fixture
def invalid_smiles():
    return "invalid"


def test_invalid_canonical_smiles(invalid_smiles):
    with pytest.raises(InvalidInputException):
        get_ob_canonical_SMILES(invalid_smiles)


def test_invalid_inchi(invalid_smiles):
    with pytest.raises(InvalidInputException):
        get_ob_InChI(invalid_smiles)


def test_invalid_smiles_3d():
    smiles = "CCC[R]"
    with pytest.raises(InvalidInputException):
        get_ob_mol(smiles, threeD=True)


def test_invalid_2d(invalid_smiles):
    with pytest.raises(InvalidInputException):
        get_ob_mol(invalid_smiles)


def test_mol_weight_violation():
    smiles = "FC=1C(F)=C(F)C(=C(F)C1F)CON=CC=C(C)C(Br)CC(Br)C(=C)C"
    mol = Chem.MolFromSmiles(smiles)
    violations = check_RO5_violations(mol)
    assert violations == 2


def test_num_h_acceptors_violation():
    smiles = "C(=O)(O)C(=O)OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)O"
    mol = Chem.MolFromSmiles(smiles)
    violations = check_RO5_violations(mol)
    assert violations == 1


def test_num_h_donors_violation():
    smiles = "NCCCNCCCCNCCCNCNCCCNCCCCNCCCNC"
    mol = Chem.MolFromSmiles(smiles)
    violations = check_RO5_violations(mol)
    assert violations == 1


def test_exception_handling():
    smiles = "C1#CC#C1"
    mol = Chem.MolFromSmiles(smiles)
    result = get_3d_conformers(mol)
    assert result is not None


def test_valid_smiles():
    smiles = "C1=CC=CC=C1"
    result = is_valid_molecule(smiles)
    assert result == "smiles"


def test_valid_molblock():
    molblock = """
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1261    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8522    0.7071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8522    1.4142    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1261    2.1213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4000    2.1213    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
"""
    result = is_valid_molecule(molblock)
    assert result == "mol"


def test_invalid_molblock():
    molblock = """
     RDKit          2D
M  END
"""
    result = is_valid_molecule(molblock)
    assert result is False


def test_invalid_input(invalid_smiles):
    result = is_valid_molecule(invalid_smiles)
    assert result is False


def test_no_stereochemistry():
    smiles = "C1=CC=CC=C1"
    mol = Chem.MolFromSmiles(smiles)
    assert has_stereo_defined(mol) is False


def test_chiral_center():
    smiles = "C[C@H](Cl)Br"
    mol = Chem.MolFromSmiles(smiles)
    assert has_stereo_defined(mol) is True


def test_tetrahedral_stereochemistry():
    smiles = "C[C@@H]1CCCC[C@H]1Br"
    mol = Chem.MolFromSmiles(smiles)
    assert has_stereo_defined(mol) is True


def test_invalid_molecule():
    mol = None
    assert has_stereo_defined(mol) is False
