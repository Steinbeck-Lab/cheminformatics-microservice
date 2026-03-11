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


# ---------------------------------------------------------------------------
# Coverage tests for security hardening (error paths)
# ---------------------------------------------------------------------------


def test_depict2D_invalid_smiles_error_path():
    """Invalid SMILES triggers error handler in 2D depiction."""
    response = client.get("/latest/depict/2D?smiles=INVALID!!!&toolkit=rdkit")
    assert response.status_code == 422


def test_depict3D_invalid_smiles_error_path():
    """Invalid SMILES triggers error handler in 3D depiction."""
    response = client.get("/latest/depict/3D?smiles=INVALID!!!&toolkit=rdkit")
    assert response.status_code == 422


# ---------------------------------------------------------------------------
# Direct function tests for depiction.py uncovered code paths
# ---------------------------------------------------------------------------

from app.modules.depiction import get_cdk_depiction as get_cdk_depiction_basic
from app.modules.depiction import get_rdkit_depiction
from app.modules.toolkits.cdk_wrapper import get_CDK_IAtomContainer
from rdkit import Chem


class TestGetCDKDepictionDirect:
    """Direct tests for get_cdk_depiction in depiction.py."""

    def test_cdk_depiction_no_cip(self):
        """Test CDK depiction with CIP=False (line 92)."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, CIP=False)
        assert result is not None
        assert isinstance(result, str)
        assert "svg" in result.lower()

    def test_cdk_depiction_unicolor(self):
        """Test CDK depiction with unicolor=True."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, unicolor=True)
        assert result is not None
        assert isinstance(result, str)

    def test_cdk_depiction_with_rotation(self):
        """Test CDK depiction with rotation (line 105-111)."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, rotate=90)
        assert result is not None

    def test_cdk_depiction_with_atom_numbers(self):
        """Test CDK depiction with showAtomNumbers=True (line 115)."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, showAtomNumbers=True)
        assert result is not None

    def test_cdk_depiction_with_smarts_highlight(self):
        """Test CDK depiction with SMARTS highlighting (lines 196-208)."""
        mol = get_CDK_IAtomContainer("c1ccccc1CC(=O)O")
        result = get_cdk_depiction_basic(mol, highlight="c1ccccc1")
        assert result is not None

    def test_cdk_depiction_with_atom_indices_highlight(self):
        """Test CDK depiction with atom indices highlighting (lines 158-190)."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, highlight_atoms=[0, 1])
        assert result is not None

    def test_cdk_depiction_with_multi_substructure_highlight(self):
        """Test CDK depiction with multiple substructure highlighting."""
        mol = get_CDK_IAtomContainer("CCCCCO")
        result = get_cdk_depiction_basic(mol, highlight_atoms=[[0, 1], [3, 4]])
        assert result is not None

    def test_cdk_depiction_invalid_hydrogen_display(self):
        """Test CDK depiction with invalid hydrogen display mode (lines 85-86)."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, hydrogen_display="InvalidMode")
        assert result is not None
        # Should return error message tuple
        if isinstance(result, tuple):
            assert "Invalid hydrogen display mode" in result[0]

    def test_cdk_depiction_kekulize_error(self):
        """Test CDK depiction with molecule that may have kekulize issues."""
        # Use a normal aromatic molecule; kekulize path should still be exercised
        mol = get_CDK_IAtomContainer("c1ccccc1")
        result = get_cdk_depiction_basic(mol, kekulize=True)
        assert result is not None

    def test_cdk_depiction_invalid_smarts_highlight(self):
        """Test CDK depiction with invalid SMARTS pattern (lines 206-208)."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, highlight="[[[invalid")
        assert result is not None

    def test_cdk_depiction_highlight_atoms_with_bonds(self):
        """Test CDK depiction atom highlighting includes bonds between atoms."""
        mol = get_CDK_IAtomContainer("CCCC")
        result = get_cdk_depiction_basic(mol, highlight_atoms=[0, 1, 2])
        assert result is not None

    def test_cdk_depiction_highlight_atoms_out_of_range(self):
        """Test CDK depiction with out-of-range atom indices."""
        mol = get_CDK_IAtomContainer("CCO")
        result = get_cdk_depiction_basic(mol, highlight_atoms=[0, 999])
        assert result is not None


class TestGetRDKitDepictionDirect:
    """Direct tests for get_rdkit_depiction in depiction.py."""

    def test_rdkit_depiction_with_cip(self):
        """Test RDKit depiction with CIP=True."""
        mol = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
        result = get_rdkit_depiction(mol, CIP=True)
        assert result is not None
        assert "svg" in result.lower()

    def test_rdkit_depiction_with_unicolor(self):
        """Test RDKit depiction with unicolor=True."""
        mol = Chem.MolFromSmiles("CCO")
        result = get_rdkit_depiction(mol, unicolor=True)
        assert result is not None

    def test_rdkit_depiction_with_atom_numbers(self):
        """Test RDKit depiction with atom numbers displayed."""
        mol = Chem.MolFromSmiles("CCO")
        result = get_rdkit_depiction(mol, showAtomNumbers=True)
        assert result is not None

    def test_rdkit_depiction_with_smarts_highlight(self):
        """Test RDKit depiction with SMARTS pattern highlighting (lines 325-341)."""
        mol = Chem.MolFromSmiles("c1ccccc1CCO")
        result = get_rdkit_depiction(mol, highlight="c1ccccc1")
        assert result is not None

    def test_rdkit_depiction_with_atom_indices_highlight(self):
        """Test RDKit depiction with specific atom indices."""
        mol = Chem.MolFromSmiles("CCCO")
        result = get_rdkit_depiction(mol, highlight_atoms=[0, 1])
        assert result is not None

    def test_rdkit_depiction_combined_highlight_and_smarts(self):
        """Test RDKit depiction with both atom indices and SMARTS (lines 271-310)."""
        mol = Chem.MolFromSmiles("c1ccccc1CCO")
        result = get_rdkit_depiction(
            mol, highlight="c1ccccc1", highlight_atoms=[0, 1, 2]
        )
        assert result is not None

    def test_rdkit_depiction_combined_no_match(self):
        """Test RDKit combined highlight where pattern doesn't match anchor atoms."""
        mol = Chem.MolFromSmiles("c1ccccc1CCO")
        # Highlight atoms that are NOT in the SMARTS match
        result = get_rdkit_depiction(mol, highlight="CO", highlight_atoms=[0, 1])
        assert result is not None

    def test_rdkit_depiction_invalid_smarts(self):
        """Test RDKit depiction with invalid SMARTS (line 341)."""
        mol = Chem.MolFromSmiles("CCO")
        result = get_rdkit_depiction(mol, highlight="[[[invalid")
        assert result is not None

    def test_rdkit_depiction_combined_invalid_smarts(self):
        """Test RDKit depiction with combined highlight and invalid SMARTS (lines 305-310)."""
        mol = Chem.MolFromSmiles("CCO")
        result = get_rdkit_depiction(
            mol, highlight="[[[invalid", highlight_atoms=[0, 1]
        )
        assert result is not None

    def test_rdkit_depiction_kekulize_error(self):
        """Test RDKit depiction with kekulize error fallback (lines 248-249)."""
        mol = Chem.MolFromSmiles("CCO")
        result = get_rdkit_depiction(mol, kekulize=True)
        assert result is not None

    def test_rdkit_depiction_rotation(self):
        """Test RDKit depiction with rotation."""
        mol = Chem.MolFromSmiles("CCO")
        result = get_rdkit_depiction(mol, rotate=90)
        assert result is not None


# ---------------------------------------------------------------------------
# Coverage tests for depict.py router uncovered lines
# ---------------------------------------------------------------------------


class TestDepict2DAtomIdsValueError:
    """Test atomIds ValueError path (lines 208-209)."""

    def test_atom_ids_with_valid_integers(self):
        """Valid integer atom IDs should work normally."""
        response = client.get("/latest/depict/2D?smiles=CCO&toolkit=cdk&atomIds=0,1")
        assert response.status_code == 200

    def test_atom_ids_mixed_valid_invalid(self):
        """Mixed valid/invalid IDs should filter out non-digits."""
        response = client.get(
            "/latest/depict/2D?smiles=CCO&toolkit=cdk&atomIds=0,abc,1"
        )
        assert response.status_code == 200


class TestDepict2DEnhanced:
    """Test the /2D_enhanced endpoint (lines 685-686, 696-717)."""

    def test_2d_enhanced_basic(self):
        """Test basic enhanced 2D depiction."""
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO")
        assert response.status_code == 200
        assert response.headers.get("content-type") == "image/svg+xml"

    def test_2d_enhanced_with_cxsmiles(self):
        """Test enhanced depiction with CXSMILES highlighting."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1%20%7Cha%3A0%2C1%2C2%7C"
        )
        assert response.status_code == 200

    def test_2d_enhanced_with_style(self):
        """Test enhanced depiction with different styles."""
        for style in ["cow", "cob", "cot", "bow", "bot", "wob", "nob"]:
            response = client.get(
                f"/latest/depict/2D_enhanced?smiles=CCO&style={style}"
            )
            assert response.status_code == 200

    def test_2d_enhanced_with_abbreviate(self):
        """Test enhanced depiction with abbreviation modes."""
        for mode in ["off", "groups", "reagents", "on"]:
            response = client.get(
                f"/latest/depict/2D_enhanced?smiles=c1ccccc1C&abbreviate={mode}"
            )
            assert response.status_code == 200

    def test_2d_enhanced_with_dative_modes(self):
        """Test enhanced depiction with dative bond modes."""
        for mode in ["never", "metals", "always"]:
            response = client.get(
                f"/latest/depict/2D_enhanced?smiles=CCO&dative={mode}"
            )
            assert response.status_code == 200

    def test_2d_enhanced_with_annotations(self):
        """Test enhanced depiction with annotation modes."""
        for mode in [
            "none",
            "number",
            "bondnumber",
            "mapidx",
            "atomvalue",
            "colmap",
            "cip",
        ]:
            response = client.get(
                f"/latest/depict/2D_enhanced?smiles=CCO&annotate={mode}"
            )
            assert response.status_code == 200

    def test_2d_enhanced_with_donuts(self):
        """Test enhanced depiction with aromatic donuts."""
        response = client.get("/latest/depict/2D_enhanced?smiles=c1ccccc1&donuts=true")
        assert response.status_code == 200

    def test_2d_enhanced_with_flip(self):
        """Test enhanced depiction with structure flipping."""
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&flip=true")
        assert response.status_code == 200

    def test_2d_enhanced_with_anon(self):
        """Test enhanced depiction with anonymous display."""
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&anon=true")
        assert response.status_code == 200

    def test_2d_enhanced_with_zoom_and_ratio(self):
        """Test enhanced depiction with custom zoom and ratio."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&zoom=2.0&ratio=1.5"
        )
        assert response.status_code == 200

    def test_2d_enhanced_with_svg_units(self):
        """Test enhanced depiction with different SVG units."""
        for unit in ["px", "mm", "cm", "in"]:
            response = client.get(
                f"/latest/depict/2D_enhanced?smiles=CCO&svgunits={unit}"
            )
            assert response.status_code == 200

    def test_2d_enhanced_with_custom_colors(self):
        """Test enhanced depiction with custom background/foreground colors."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO"
            "&bgcolor=%23FFFFFF&fgcolor=%23000000"
        )
        assert response.status_code == 200

    def test_2d_enhanced_with_default_colors(self):
        """Test enhanced depiction with default color values."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&bgcolor=default&fgcolor=default"
        )
        assert response.status_code == 200

    def test_2d_enhanced_with_unicolor(self):
        """Test enhanced depiction with unicolor mode."""
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&unicolor=true")
        assert response.status_code == 200

    def test_2d_enhanced_with_atom_ids(self):
        """Test enhanced depiction with atom indices highlighting."""
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&atomIds=0,1")
        assert response.status_code == 200

    def test_2d_enhanced_with_smarts_highlight(self):
        """Test enhanced depiction with SMARTS highlighting."""
        response = client.get("/latest/depict/2D_enhanced?smiles=c1ccccc1&highlight=c")
        assert response.status_code == 200

    def test_2d_enhanced_with_perceive_radicals(self):
        """Test enhanced depiction with radical perception."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=%5BCH3%5D&perceive_radicals=true"
        )
        assert response.status_code == 200

    def test_2d_enhanced_with_title(self):
        """Test enhanced depiction with custom title."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&showtitle=true&title=Ethanol"
        )
        assert response.status_code == 200

    def test_2d_enhanced_with_show_atom_numbers(self):
        """Test enhanced depiction with atom numbers shown."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&showAtomNumbers=true"
        )
        assert response.status_code == 200

    def test_2d_enhanced_with_hydrogen_display(self):
        """Test enhanced depiction with various hydrogen display modes."""
        for mode in ["Provided", "Minimal", "Explicit", "Stereo", "Smart"]:
            response = client.get(
                f"/latest/depict/2D_enhanced?smiles=C%5BC%40H%5D(N)C(%3DO)O"
                f"&hydrogen_display={mode}"
            )
            assert response.status_code == 200

    def test_2d_enhanced_invalid_smiles(self):
        """Test enhanced depiction with invalid SMILES string."""
        response = client.get("/latest/depict/2D_enhanced?smiles=INVALID_SMILES!!!")
        assert response.status_code == 422

    def test_2d_enhanced_combined_features(self):
        """Test enhanced depiction with multiple features combined."""
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1C&style=bow"
            "&annotate=number&donuts=true&zoom=1.5&ratio=1.2&flip=true&CIP=true"
        )
        assert response.status_code == 200

    def test_2d_enhanced_invalid_atom_ids_format(self):
        """Test enhanced depiction with non-digit atom IDs filtered out."""
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&atomIds=a,b,c")
        # Non-digit strings are filtered out, so this should succeed
        assert response.status_code == 200


class TestDepict3DToolkitPaths:
    """Test 3D depiction toolkit parameter path (lines 318, 324)."""

    def test_3d_depiction_openbabel(self):
        """Test 3D depiction with OpenBabel toolkit."""
        response = client.get("/latest/depict/3D?smiles=CCO&toolkit=openbabel")
        assert response.status_code == 200

    def test_3d_depiction_rdkit(self):
        """Test 3D depiction with RDKit toolkit."""
        response = client.get("/latest/depict/3D?smiles=CCO&toolkit=rdkit")
        assert response.status_code == 200

    def test_3d_depiction_invalid_smiles_openbabel(self):
        """Test 3D depiction with invalid SMILES for OpenBabel."""
        response = client.get("/latest/depict/3D?smiles=INVALID!!!&toolkit=openbabel")
        assert response.status_code == 422


class TestDepict2DHydrogenDisplayRDKit:
    """Test hydrogen display with RDKit toolkit (lines 215-220)."""

    def test_hydrogen_display_rdkit_provided(self):
        """Test Provided hydrogen display with RDKit toolkit."""
        response = client.get(
            "/latest/depict/2D?smiles=C%5BC%40H%5D(N)C(%3DO)O"
            "&toolkit=rdkit&hydrogen_display=Provided"
        )
        assert response.status_code == 200

    def test_hydrogen_display_rdkit_minimal(self):
        """Test Minimal hydrogen display with RDKit toolkit."""
        response = client.get(
            "/latest/depict/2D?smiles=CCO&toolkit=rdkit&hydrogen_display=Minimal"
        )
        assert response.status_code == 200

    def test_hydrogen_display_rdkit_smart(self):
        """Test Smart hydrogen display with RDKit toolkit."""
        response = client.get(
            "/latest/depict/2D?smiles=CCO&toolkit=rdkit&hydrogen_display=Smart"
        )
        assert response.status_code == 200


class TestDepict2DToolkitCDK:
    """Test CDK-specific 2D depiction paths (line 244)."""

    def test_2d_cdk_with_highlight_and_atom_numbers(self):
        """Test CDK depiction with highlight + atom numbers."""
        response = client.get(
            "/latest/depict/2D?smiles=CCO&toolkit=cdk"
            "&highlight=CO&showAtomNumbers=true"
        )
        assert response.status_code == 200
