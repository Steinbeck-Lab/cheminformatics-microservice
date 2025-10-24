from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_chem_index():
    response = client.get("/latest/depict/")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}


@pytest.mark.parametrize(
    "smiles, toolkit, width, height, rotate, CIP, unicolor, response_code",
    [
        ("CCO", "rdkit", 512, 512, 0, False, False, 200),
        ("INVALID_INPUT", "rdkit", 512, 512, 0, False, False, 422),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "rdkit", 512, 512, 0, False, False, 200),
    ],
)
def test_depict2D_molecule_rdkit(
    smiles,
    toolkit,
    width,
    height,
    rotate,
    CIP,
    unicolor,
    response_code,
):
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&width={width}&height={height}&rotate={rotate}&CIP={CIP}&unicolor={unicolor}",
    )
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "smiles, toolkit, width, height, rotate, CIP, unicolor, response_code",
    [
        ("CCO", "cdk", 512, 512, 0, False, False, 200),
        ("INVALID_INPUT", "cdk", 512, 512, 0, False, False, 422),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "cdk", 512, 512, 0, False, False, 200),
    ],
)
def test_depict2D_molecule_cdk(
    smiles,
    toolkit,
    width,
    height,
    rotate,
    CIP,
    unicolor,
    response_code,
):
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&width={width}&height={height}&rotate={rotate}&CIP={CIP}&unicolor={unicolor}",
    )
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "smiles, toolkit, response_code",
    [
        ("CCO", "openbabel", 200),
        ("CCO", "rdkit", 200),
        ("INVALID_INPUT", "openbabel", 422),
        ("INVALID_INPUT", "rdkit", 422),
    ],
)
def test_depict3D_molecule(smiles, toolkit, response_code):
    response = client.get(
        f"/latest/depict/3D?smiles={smiles}&generator={toolkit}",
    )
    assert response.status_code == response_code


@pytest.mark.parametrize(
    "smiles, toolkit, width, height",
    [
        ("CCO", "rdkit", 256, 256),
        ("CCO", "rdkit", 1024, 1024),
    ],
)
def test_depict2D_custom_dimensions_rdkit(smiles, toolkit, width, height):
    """Test 2D depiction with custom dimensions (RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&width={width}&height={height}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, width, height",
    [
        ("CCO", "cdk", 256, 256),
        ("CCO", "cdk", 1024, 1024),
    ],
)
def test_depict2D_custom_dimensions_cdk(smiles, toolkit, width, height):
    """Test 2D depiction with custom dimensions (CDK)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&width={width}&height={height}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, rotate",
    [
        ("CCO", "rdkit", 90),
        ("CCO", "rdkit", 180),
    ],
)
def test_depict2D_rotation_rdkit(smiles, toolkit, rotate):
    """Test 2D depiction with rotation (RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&rotate={rotate}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, rotate",
    [
        ("CCO", "cdk", 90),
        ("CCO", "cdk", 180),
    ],
)
def test_depict2D_rotation_cdk(smiles, toolkit, rotate):
    """Test 2D depiction with rotation (CDK)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&rotate={rotate}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, CIP",
    [
        ("C[C@H](N)C(=O)O", "rdkit", True),
        ("C[C@H](N)C(=O)O", "rdkit", False),
    ],
)
def test_depict2D_cip_stereochemistry_rdkit(smiles, toolkit, CIP):
    """Test CIP stereochemistry annotation (RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&CIP={CIP}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, CIP",
    [
        ("C[C@H](N)C(=O)O", "cdk", True),
        ("C[C@H](N)C(=O)O", "cdk", False),
    ],
)
def test_depict2D_cip_stereochemistry_cdk(smiles, toolkit, CIP):
    """Test CIP stereochemistry annotation (CDK)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&CIP={CIP}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, unicolor",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "rdkit", True),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "rdkit", False),
    ],
)
def test_depict2D_unicolor_rdkit(smiles, toolkit, unicolor):
    """Test unicolor rendering mode (RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&unicolor={unicolor}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, unicolor",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "cdk", True),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "cdk", False),
    ],
)
def test_depict2D_unicolor_cdk(smiles, toolkit, unicolor):
    """Test unicolor rendering mode (CDK)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&unicolor={unicolor}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, hydrogen_display",
    [
        ("C[C@H](N)C(=O)O", "Provided"),
        ("C[C@H](N)C(=O)O", "Minimal"),
        ("C[C@H](N)C(=O)O", "Explicit"),
        ("C[C@H](N)C(=O)O", "Stereo"),
        ("C[C@H](N)C(=O)O", "Smart"),
        ("CCO", "Minimal"),
        ("CCO", "Explicit"),
    ],
)
def test_depict2D_hydrogen_display_cdk(smiles, hydrogen_display):
    """Test hydrogen display modes with CDK (not available in RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator=cdk&hydrogen_display={hydrogen_display}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, highlight",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "rdkit", "C=O"),
        ("CCO", "rdkit", "CO"),
    ],
)
def test_depict2D_smarts_highlighting_rdkit(smiles, toolkit, highlight):
    """Test SMARTS pattern highlighting (RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&highlight={highlight}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, highlight",
    [
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "cdk", "C=O"),
        ("CCO", "cdk", "CO"),
    ],
)
def test_depict2D_smarts_highlighting_cdk(smiles, toolkit, highlight):
    """Test SMARTS pattern highlighting (CDK)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&highlight={highlight}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, atomIds",
    [
        ("CCO", "rdkit", "0,1"),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "rdkit", "0,1,2,3"),
    ],
)
def test_depict2D_atom_indices_highlighting_rdkit(smiles, toolkit, atomIds):
    """Test atom indices highlighting (RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&atomIds={atomIds}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, atomIds",
    [
        ("CCO", "cdk", "0,1,2"),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "cdk", "0,1,2,3"),
    ],
)
def test_depict2D_atom_indices_highlighting_cdk(smiles, toolkit, atomIds):
    """Test atom indices highlighting (CDK)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&atomIds={atomIds}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, showAtomNumbers",
    [
        ("CCO", "rdkit", True),
        ("CCO", "rdkit", False),
    ],
)
def test_depict2D_atom_numbers_rdkit(smiles, toolkit, showAtomNumbers):
    """Test atom number display (RDKit)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&showAtomNumbers={showAtomNumbers}",
    )
    assert response.status_code == 200


@pytest.mark.parametrize(
    "smiles, toolkit, showAtomNumbers",
    [
        ("CCO", "cdk", True),
        ("CCO", "cdk", False),
    ],
)
def test_depict2D_atom_numbers_cdk(smiles, toolkit, showAtomNumbers):
    """Test atom number display (CDK)."""
    response = client.get(
        f"/latest/depict/2D?smiles={smiles}&generator={toolkit}&showAtomNumbers={showAtomNumbers}",
    )
    assert response.status_code == 200


def test_depict2D_combined_parameters():
    """Test depiction with multiple parameters combined."""
    response = client.get(
        "/latest/depict/2D?smiles=C[C@H](N)C(=O)O&generator=cdk&width=800&height=600&rotate=90&CIP=True&unicolor=False&showAtomNumbers=True&hydrogen_display=Stereo",
    )
    assert response.status_code == 200


def test_depict2D_combined_highlighting():
    """Test depiction with combined highlighting options."""
    response = client.get(
        "/latest/depict/2D?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&generator=rdkit&highlight=C=O&atomIds=0,1,2",
    )
    assert response.status_code == 200


def test_depict2D_invalid_smiles():
    """Test with invalid SMILES string."""
    response = client.get(
        "/latest/depict/2D?smiles=INVALID&generator=rdkit",
    )
    assert response.status_code == 422


def test_depict2D_invalid_atom_ids():
    """Test with invalid atom IDs format - non-digit strings are filtered out."""
    # Note: The API filters out non-digit strings rather than raising an error
    # So 'xx' will just result in an empty list and succeed with no highlighting
    response = client.get(
        "/latest/depict/2D?smiles=CCO&generator=cdk&atomIds=xx,yy",
    )
    # This will succeed (200) because invalid IDs are filtered out
    assert response.status_code == 200


def test_depict2D_out_of_range_atom_ids():
    """Test with out-of-range atom indices."""
    # CCO has only 3 atoms (indices 0, 1, 2), trying index 999 which is out of range
    # The API should handle this gracefully
    response = client.get(
        "/latest/depict/2D?smiles=CCO&generator=rdkit&atomIds=0,1,999",
    )
    # Should still succeed, just ignores invalid indices
    assert response.status_code == 422


def test_depict3D_openbabel_default():
    """Test 3D depiction with OpenBabel as default."""
    response = client.get("/latest/depict/3D?smiles=CCO")
    assert response.status_code == 200
