"""Tests for the CDX/CDXML → MOL conversion feature.

Covers:
  - Unit tests for ``convert_cdx_to_mol`` in rdkit_wrapper
  - API endpoint tests for POST /convert/cdx-to-mol
"""

from __future__ import annotations

import os
import tempfile

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from app.main import app
from app.modules.toolkits.rdkit_wrapper import convert_cdx_to_mol

client = TestClient(app)

# Absolute path so tests run regardless of cwd
CDX_PATH = os.path.join(os.path.dirname(__file__), "ketcher.cdx")

# ── Fixtures ────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def ketcher_cdx_bytes() -> bytes:
    """Raw bytes of the real binary ketcher.cdx test file (naphthalene)."""
    with open(CDX_PATH, "rb") as fh:
        return fh.read()


@pytest.fixture(scope="module")
def ketcher_cdxml_bytes() -> bytes:
    """CDXML equivalent of the ketcher.cdx file, generated via pybel."""
    from openbabel import pybel

    mol = pybel.readstring("smi", "c1ccc2ccccc2c1")  # naphthalene
    with tempfile.NamedTemporaryFile(suffix=".cdxml", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        mol.write("cdxml", tmp_path, overwrite=True)
        with open(tmp_path, "rb") as fh:
            return fh.read()
    finally:
        os.unlink(tmp_path)


# ── Unit: convert_cdx_to_mol ────────────────────────────────────────────────


class TestConvertCdxToMolUnit:
    """Direct unit tests for the rdkit_wrapper.convert_cdx_to_mol function."""

    def test_binary_cdx_returns_string(self, ketcher_cdx_bytes):
        result = convert_cdx_to_mol(ketcher_cdx_bytes, fmt="cdx")
        assert isinstance(result, str)

    def test_binary_cdx_contains_mol_end(self, ketcher_cdx_bytes):
        result = convert_cdx_to_mol(ketcher_cdx_bytes, fmt="cdx")
        assert "M  END" in result

    def test_binary_cdx_rdkit_parseable(self, ketcher_cdx_bytes):
        molblock = convert_cdx_to_mol(ketcher_cdx_bytes, fmt="cdx")
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        assert mol is not None

    def test_binary_cdx_atom_count(self, ketcher_cdx_bytes):
        """ketcher.cdx contains naphthalene (10 heavy atoms)."""
        molblock = convert_cdx_to_mol(ketcher_cdx_bytes, fmt="cdx")
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        assert mol.GetNumAtoms() == 10

    def test_binary_cdx_smiles_roundtrip(self, ketcher_cdx_bytes):
        """SMILES round-trip must match canonical naphthalene."""
        molblock = convert_cdx_to_mol(ketcher_cdx_bytes, fmt="cdx")
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        result_smi = Chem.MolToSmiles(mol)
        expected_smi = Chem.MolToSmiles(Chem.MolFromSmiles("c1ccc2ccccc2c1"))
        assert result_smi == expected_smi

    def test_cdxml_returns_mol_block(self, ketcher_cdxml_bytes):
        result = convert_cdx_to_mol(ketcher_cdxml_bytes, fmt="cdxml")
        assert "M  END" in result

    def test_cdxml_rdkit_parseable(self, ketcher_cdxml_bytes):
        molblock = convert_cdx_to_mol(ketcher_cdxml_bytes, fmt="cdxml")
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_invalid_bytes_raises_value_error(self):
        with pytest.raises(ValueError, match="No valid molecules"):
            convert_cdx_to_mol(b"this is not a cdx file", fmt="cdx")

    def test_empty_bytes_raises_value_error(self):
        with pytest.raises(ValueError):
            convert_cdx_to_mol(b"", fmt="cdx")

    def test_molblock_does_not_have_leading_strip(self, ketcher_cdx_bytes):
        """rstrip() must preserve the leading blank mol-name line required by
        V2000 format — the 4th line must start with a space-padded atom count,
        not a coordinate value."""
        molblock = convert_cdx_to_mol(ketcher_cdx_bytes, fmt="cdx")
        lines = molblock.splitlines()
        # V2000 counts line is line index 3 (0-based)
        assert len(lines) >= 4
        counts_line = lines[3]
        # Must be parseable as a counts line, not a coordinate line
        assert "V2000" in counts_line or counts_line.strip()[:6].strip().isdigit()


# ── API: POST /convert/cdx-to-mol ───────────────────────────────────────────


class TestCdxToMolEndpoint:
    """Integration tests for POST /latest/convert/cdx-to-mol."""

    # ── Happy path ──────────────────────────────────────────────────────────

    def test_status_200_with_valid_cdx(self, ketcher_cdx_bytes):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("ketcher.cdx", ketcher_cdx_bytes, "chemical/x-cdx")},
        )
        assert response.status_code == 200

    def test_response_json_has_molblock_key(self, ketcher_cdx_bytes):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("ketcher.cdx", ketcher_cdx_bytes, "chemical/x-cdx")},
        )
        assert "molblock" in response.json()

    def test_returned_molblock_has_mol_end(self, ketcher_cdx_bytes):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("ketcher.cdx", ketcher_cdx_bytes, "chemical/x-cdx")},
        )
        molblock = response.json()["molblock"]
        assert "M  END" in molblock

    def test_returned_molblock_rdkit_parseable(self, ketcher_cdx_bytes):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("ketcher.cdx", ketcher_cdx_bytes, "chemical/x-cdx")},
        )
        molblock = response.json()["molblock"]
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        assert mol is not None

    def test_atom_count_matches_naphthalene(self, ketcher_cdx_bytes):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("ketcher.cdx", ketcher_cdx_bytes, "chemical/x-cdx")},
        )
        molblock = response.json()["molblock"]
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        assert mol.GetNumAtoms() == 10

    def test_smiles_roundtrip_via_endpoint(self, ketcher_cdx_bytes):
        """CDX → MOL (endpoint) → SMILES must equal canonical naphthalene."""
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("ketcher.cdx", ketcher_cdx_bytes, "chemical/x-cdx")},
        )
        molblock = response.json()["molblock"]
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        result_smi = Chem.MolToSmiles(mol)
        expected_smi = Chem.MolToSmiles(Chem.MolFromSmiles("c1ccc2ccccc2c1"))
        assert result_smi == expected_smi

    def test_cdxml_file_accepted(self, ketcher_cdxml_bytes):
        """CDXML files (text/XML ChemDraw format) must also be accepted."""
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("molecule.cdxml", ketcher_cdxml_bytes, "chemical/x-cdxml")},
        )
        assert response.status_code == 200
        molblock = response.json()["molblock"]
        mol = Chem.MolFromMolBlock(molblock, sanitize=True)
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_mol_to_smiles_pipeline(self, ketcher_cdx_bytes):
        """MOL returned by CDX endpoint can be fed into the molblock endpoint."""
        # Step 1: CDX → MOL
        cdx_resp = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("ketcher.cdx", ketcher_cdx_bytes, "chemical/x-cdx")},
        )
        assert cdx_resp.status_code == 200
        molblock = cdx_resp.json()["molblock"]

        # Step 2: MOL → SMILES
        mol_resp = client.post(
            "/latest/convert/molblock?toolkit=rdkit",
            data=molblock,
            headers={"Content-Type": "text/plain"},
        )
        assert mol_resp.status_code == 200
        smiles = mol_resp.text.strip('"')
        assert len(smiles) > 0

    # ── Error handling ──────────────────────────────────────────────────────

    def test_wrong_extension_returns_400(self):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("molecule.png", b"fake", "image/png")},
        )
        assert response.status_code == 400

    def test_wrong_extension_txt_returns_400(self):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("molecule.txt", b"fake", "text/plain")},
        )
        assert response.status_code == 400

    def test_wrong_extension_error_message(self):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("molecule.sdf", b"fake", "chemical/x-mdl-sdfile")},
        )
        assert response.status_code == 400
        assert "cdx" in response.json()["detail"].lower()

    def test_garbage_bytes_returns_422(self):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("bad.cdx", b"this is not a cdx file", "chemical/x-cdx")},
        )
        assert response.status_code == 422

    def test_empty_file_returns_422(self):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("empty.cdx", b"", "chemical/x-cdx")},
        )
        assert response.status_code == 422

    def test_random_bytes_returns_422(self):
        import os as _os

        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": ("random.cdx", _os.urandom(256), "chemical/x-cdx")},
        )
        assert response.status_code == 422

    # ── Parametrized matrix ─────────────────────────────────────────────────

    @pytest.mark.parametrize(
        "filename, content, expected_status",
        [
            ("mol.cdx", open(CDX_PATH, "rb").read(), 200),
            ("mol.cdx", b"not a cdx file", 422),
            ("mol.cdx", b"", 422),
            ("mol.txt", open(CDX_PATH, "rb").read(), 400),
            ("mol.sdf", open(CDX_PATH, "rb").read(), 400),
            ("mol.mol", b"fake", 400),
        ],
    )
    def test_parametrized_inputs(self, filename, content, expected_status):
        response = client.post(
            "/latest/convert/cdx-to-mol",
            files={"file": (filename, content, "application/octet-stream")},
        )
        assert response.status_code == expected_status
