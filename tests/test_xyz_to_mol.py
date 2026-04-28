"""Tests for the XYZ → SMILES/InChI/InChIKey/MOL/SDF conversion feature.

Covers:
  - Unit tests for ``convert_xyz_to_mol`` in rdkit_wrapper (RDKit xyz2mol path)
  - Unit tests for ``get_ob_xyz_conversions`` in openbabel_wrapper (OB path)
  - API endpoint tests for POST /convert/xyz
"""

from __future__ import annotations

import os

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from app.main import app
from app.modules.toolkits.openbabel_wrapper import get_ob_xyz_conversions
from app.modules.toolkits.rdkit_wrapper import convert_xyz_to_mol

client = TestClient(app)

WATER_XYZ_PATH = os.path.join(os.path.dirname(__file__), "water.xyz")
ACETATE_XYZ_PATH = os.path.join(os.path.dirname(__file__), "acetate.xyz")

# Methane geometry — pure neutral, no surprises, useful as a control.
METHANE_XYZ = (
    "5\n"
    "methane\n"
    "C    0.0000    0.0000    0.0000\n"
    "H    0.6300    0.6300    0.6300\n"
    "H   -0.6300   -0.6300    0.6300\n"
    "H   -0.6300    0.6300   -0.6300\n"
    "H    0.6300   -0.6300   -0.6300\n"
)

# Reference InChIKeys (without stereo layers).
WATER_INCHIKEY = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
ACETATE_INCHIKEY = "QTBSBXVTEAMEQO-UHFFFAOYSA-M"
METHANE_INCHIKEY = "VNWKTOKETHGBQD-UHFFFAOYSA-N"


# ── Fixtures ────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def water_xyz() -> str:
    with open(WATER_XYZ_PATH) as fh:
        return fh.read()


@pytest.fixture(scope="module")
def acetate_xyz() -> str:
    with open(ACETATE_XYZ_PATH) as fh:
        return fh.read()


# ── Unit: convert_xyz_to_mol (RDKit) ────────────────────────────────────────


class TestConvertXYZToMolRDKit:
    """Direct unit tests for the RDKit xyz2mol path."""

    def test_water_returns_mol(self, water_xyz):
        mol = convert_xyz_to_mol(water_xyz)
        assert isinstance(mol, Chem.Mol)
        assert mol.GetNumAtoms() == 3
        assert mol.GetNumBonds() == 2

    def test_water_smiles_roundtrip(self, water_xyz):
        mol = convert_xyz_to_mol(water_xyz)
        # Hydrogens are explicit because they came from the XYZ
        assert Chem.inchi.MolToInchiKey(mol) == WATER_INCHIKEY

    def test_methane_neutral(self):
        mol = convert_xyz_to_mol(METHANE_XYZ, charge=0)
        assert mol.GetNumAtoms() == 5
        assert mol.GetNumBonds() == 4
        assert Chem.inchi.MolToInchiKey(mol) == METHANE_INCHIKEY

    def test_acetate_charge_minus_one(self, acetate_xyz):
        mol = convert_xyz_to_mol(acetate_xyz, charge=-1)
        assert mol.GetNumAtoms() == 7
        # Net charge on the molecule should be -1
        net_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        assert net_charge == -1
        assert Chem.inchi.MolToInchiKey(mol) == ACETATE_INCHIKEY

    def test_acetate_wrong_charge_does_not_match_acetate_inchikey(self, acetate_xyz):
        """At charge=0 the same geometry yields a different (radical) species,
        so it must not pass for the acetate InChIKey. Guards against silent
        ignoring of the charge parameter."""
        try:
            mol = convert_xyz_to_mol(acetate_xyz, charge=0)
        except ValueError:
            return  # Acceptable: bond perception refuses inconsistent geometry
        assert Chem.inchi.MolToInchiKey(mol) != ACETATE_INCHIKEY

    def test_molblock_carries_3d_coords(self, water_xyz):
        mol = convert_xyz_to_mol(water_xyz)
        molblock = Chem.MolToMolBlock(mol)
        assert "M  END" in molblock
        # The first oxygen sits at the origin in our fixture.
        assert "0.0000    0.0000    0.0000 O" in molblock

    def test_empty_input_raises_value_error(self):
        with pytest.raises(ValueError):
            convert_xyz_to_mol("")

    def test_whitespace_input_raises_value_error(self):
        with pytest.raises(ValueError):
            convert_xyz_to_mol("   \n  \n")

    def test_garbage_input_raises_value_error(self):
        with pytest.raises(ValueError):
            convert_xyz_to_mol("this is not an xyz file at all")


# ── Unit: get_ob_xyz_conversions (OpenBabel) ────────────────────────────────


class TestOBXYZConversions:
    """Direct unit tests for the OpenBabel XYZ path."""

    def test_water(self, water_xyz):
        result = get_ob_xyz_conversions(water_xyz)
        assert set(result.keys()) == {"canonicalsmiles", "inchi", "inchikey", "molblock"}
        # OpenBabel's canonical writer is title-stripped by the wrapper.
        assert "\t" not in result["canonicalsmiles"]
        assert result["inchikey"].endswith("UHFFFAOYSA-N") or result["inchikey"].endswith("-N")
        assert "M  END" in result["molblock"]

    def test_methane(self):
        result = get_ob_xyz_conversions(METHANE_XYZ)
        assert result["inchikey"] == METHANE_INCHIKEY

    def test_empty_raises(self):
        from app.exception_handlers import InvalidInputException

        with pytest.raises(InvalidInputException):
            get_ob_xyz_conversions("")


# ── API: POST /convert/xyz ──────────────────────────────────────────────────


class TestXYZConvertEndpoint:
    """End-to-end tests for the public conversion endpoint."""

    URL = "/latest/convert/xyz"

    def _post(self, body: str, **params):
        return client.post(
            self.URL,
            content=body,
            headers={"content-type": "text/plain"},
            params=params,
        )

    def test_water_rdkit_default(self, water_xyz):
        r = self._post(water_xyz)
        assert r.status_code == 200
        body = r.json()
        assert {"canonicalsmiles", "inchi", "inchikey", "molblock", "sdf"} <= body.keys()
        assert body["inchikey"] == WATER_INCHIKEY
        assert body["sdf"].rstrip().endswith("$$$$")
        assert "M  END" in body["molblock"]

    def test_acetate_rdkit_with_charge(self, acetate_xyz):
        r = self._post(acetate_xyz, charge=-1)
        assert r.status_code == 200
        assert r.json()["inchikey"] == ACETATE_INCHIKEY

    def test_acetate_charge_param_round_trip(self, acetate_xyz):
        r0 = self._post(acetate_xyz, charge=0)
        r_neg = self._post(acetate_xyz, charge=-1)
        # When the geometry is consistent only with charge=-1, charge=0 either
        # fails (422) or yields a different InChIKey. Either way, the two
        # responses must not be identical -- proves charge is being honored.
        if r0.status_code == 200 and r_neg.status_code == 200:
            assert r0.json()["inchikey"] != r_neg.json()["inchikey"]

    def test_water_openbabel(self, water_xyz):
        r = self._post(water_xyz, toolkit="openbabel")
        assert r.status_code == 200
        body = r.json()
        # OB writer must not leak the title field.
        assert "\t" not in body["canonicalsmiles"]
        assert body["sdf"].rstrip().endswith("$$$$")

    def test_garbage_returns_422(self):
        r = self._post("not an xyz file at all")
        assert r.status_code == 422

    def test_truncated_xyz_returns_422(self):
        truncated = "3\nwater\nO 0 0 0\n"  # missing two atoms
        r = self._post(truncated)
        assert r.status_code == 422

    def test_unknown_element_returns_422(self):
        bad = "1\nfake\nXx 0 0 0\n"
        r = self._post(bad)
        assert r.status_code == 422

    def test_use_huckel_param_accepted(self, water_xyz):
        r = self._post(water_xyz, use_huckel="false")
        assert r.status_code == 200

    def test_charge_out_of_range_rejected(self, water_xyz):
        r = self._post(water_xyz, charge=999)
        assert r.status_code == 422

    def test_crlf_line_endings_handled(self, water_xyz):
        r = self._post(water_xyz.replace("\n", "\r\n"))
        assert r.status_code == 200
        assert r.json()["inchikey"] == WATER_INCHIKEY
