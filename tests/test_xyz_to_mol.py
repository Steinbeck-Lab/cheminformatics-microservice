"""Tests for the XYZ batch conversion feature.

Covers:
  - convert_xyz_to_mol two-tier perception (rdkit_wrapper)
  - get_ob_xyz_conversions single-frame helper (openbabel_wrapper)
  - POST /convert/xyz batch endpoint with per-frame error isolation
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
METALLOCENES_PATH = os.path.join(
    os.path.dirname(__file__), "metallocenes_truncated.xyz"
)

METHANE_XYZ = (
    "5\n"
    "methane\n"
    "C    0.0000    0.0000    0.0000\n"
    "H    0.6300    0.6300    0.6300\n"
    "H   -0.6300   -0.6300    0.6300\n"
    "H   -0.6300    0.6300   -0.6300\n"
    "H    0.6300   -0.6300   -0.6300\n"
)

WATER_INCHIKEY = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
ACETATE_INCHIKEY = "QTBSBXVTEAMEQO-UHFFFAOYSA-M"
METHANE_INCHIKEY = "VNWKTOKETHGBQD-UHFFFAOYSA-N"


@pytest.fixture(scope="module")
def water_xyz() -> str:
    with open(WATER_XYZ_PATH) as fh:
        return fh.read()


@pytest.fixture(scope="module")
def acetate_xyz() -> str:
    with open(ACETATE_XYZ_PATH) as fh:
        return fh.read()


@pytest.fixture(scope="module")
def metallocenes_xyz() -> str:
    with open(METALLOCENES_PATH) as fh:
        return fh.read()


# ── Unit: convert_xyz_to_mol returns (mol, method) ──────────────────────────


class TestConvertXYZToMolRDKit:
    def test_water_returns_bond_orders(self, water_xyz):
        mol, method = convert_xyz_to_mol(water_xyz)
        assert method == "bond_orders"
        assert mol.GetNumBonds() == 2
        assert Chem.inchi.MolToInchiKey(mol) == WATER_INCHIKEY

    def test_methane_returns_bond_orders(self):
        mol, method = convert_xyz_to_mol(METHANE_XYZ)
        assert method == "bond_orders"
        assert Chem.inchi.MolToInchiKey(mol) == METHANE_INCHIKEY

    def test_acetate_with_charge_returns_bond_orders(self, acetate_xyz):
        mol, method = convert_xyz_to_mol(acetate_xyz, charge=-1)
        assert method == "bond_orders"
        assert Chem.inchi.MolToInchiKey(mol) == ACETATE_INCHIKEY

    def test_metallocene_falls_back_to_connectivity_only(self, metallocenes_xyz):
        from app.modules.toolkits.helpers import split_xyz_frames

        frames = split_xyz_frames(metallocenes_xyz)
        assert len(frames) == 3
        for frame in frames:
            mol, method = convert_xyz_to_mol(frame)
            assert method == "connectivity_only"
            assert mol.GetNumBonds() > 0

    def test_empty_input_raises_value_error(self):
        with pytest.raises(ValueError):
            convert_xyz_to_mol("")

    def test_garbage_input_raises_value_error(self):
        with pytest.raises(ValueError):
            convert_xyz_to_mol("not an xyz file at all")


# ── Unit: get_ob_xyz_conversions ────────────────────────────────────────────


class TestOBXYZConversions:
    def test_water(self, water_xyz):
        result = get_ob_xyz_conversions(water_xyz)
        assert "\t" not in result["canonicalsmiles"]
        assert "M  END" in result["molblock"]

    def test_methane(self):
        result = get_ob_xyz_conversions(METHANE_XYZ)
        assert result["inchikey"] == METHANE_INCHIKEY

    def test_empty_raises(self):
        from app.exception_handlers import InvalidInputException

        with pytest.raises(InvalidInputException):
            get_ob_xyz_conversions("")


# ── API: POST /convert/xyz batch ────────────────────────────────────────────


class TestXYZConvertEndpoint:
    URL = "/latest/convert/xyz"

    def _post(self, body: str, **params):
        return client.post(
            self.URL,
            content=body,
            headers={"content-type": "text/plain"},
            params=params,
        )

    # --- single frame ---

    def test_single_water_returns_length_one_batch(self, water_xyz):
        r = self._post(water_xyz)
        assert r.status_code == 200, r.text
        body = r.json()
        assert len(body["structures"]) == 1
        s = body["structures"][0]
        assert s["success"] is True
        assert s["index"] == 0
        assert s["method"] == "bond_orders"
        assert s["bond_orders_perceived"] is True
        assert s["inchikey"] == WATER_INCHIKEY
        assert body["sdf"].rstrip().endswith("$$$$")
        assert body["summary"] == {
            "total": 1,
            "successful": 1,
            "failed": 0,
            "bond_orders_count": 1,
            "connectivity_only_count": 0,
        }

    def test_acetate_with_charge(self, acetate_xyz):
        r = self._post(acetate_xyz, charge=-1)
        assert r.status_code == 200
        s = r.json()["structures"][0]
        assert s["inchikey"] == ACETATE_INCHIKEY
        assert s["method"] == "bond_orders"

    # --- multi-frame ---

    def test_two_frame_request_returns_two_results(self, water_xyz, acetate_xyz):
        r = self._post(water_xyz + acetate_xyz, charge=-1)
        # The second frame (acetate, charge=-1) must always succeed; the first
        # (water at charge=-1) may either fail or yield a radical/protonated
        # variant. Either way, the batch should report 2 frames.
        assert r.status_code == 200
        body = r.json()
        assert body["summary"]["total"] == 2
        # Find the acetate result (will be the one with the acetate InChIKey)
        ikeys = [s["inchikey"] for s in body["structures"] if s.get("success")]
        assert ACETATE_INCHIKEY in ikeys

    def test_metallocenes_all_connectivity_only(self, metallocenes_xyz):
        r = self._post(metallocenes_xyz)
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["summary"]["total"] == 3
        assert body["summary"]["successful"] == 3
        assert body["summary"]["connectivity_only_count"] == 3
        assert body["summary"]["bond_orders_count"] == 0
        for s in body["structures"]:
            assert s["method"] == "connectivity_only"
            assert s["bond_orders_perceived"] is False
            assert any("connectivity-only" in w.lower() for w in s["warnings"])

    # --- per-frame failure isolation ---

    def test_bad_frame_does_not_fail_whole_request(self, water_xyz):
        # Bogus 1-atom frame with an invalid element symbol; the splitter
        # accepts it, but RDKit/OpenBabel both fail to perceive bonds.
        bad_frame = "1\nbogus\nXx 0 0 0\n"
        r = self._post(water_xyz + bad_frame + water_xyz)
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["summary"]["total"] >= 2
        assert body["summary"]["successful"] >= 1
        # SDF skips failed frames; verify it parses
        suppl = Chem.SDMolSupplier()
        suppl.SetData(body["sdf"])
        n_records = sum(1 for m in suppl if m is not None)
        assert n_records == body["summary"]["successful"]

    # --- whole-request errors ---

    def test_no_frames_detected_returns_422(self):
        r = self._post("just some text\nwithout an atom count\n")
        assert r.status_code == 422

    def test_all_frames_fail_returns_422(self):
        # Two 1-atom frames with unknown element.
        bad = "1\na\nXx 0 0 0\n" * 2
        r = self._post(bad)
        # Either every frame fails (-> 422) or RDKit's parser rejects them
        # before the loop. Either way, the request must not be 200.
        assert r.status_code in (422, 400)

    def test_too_many_frames_returns_413(self, water_xyz):
        # 51 frames: cap is 50.
        r = self._post(water_xyz * 51)
        assert r.status_code == 413
        assert "Frame count exceeds limit" in r.json()["detail"]

    def test_too_large_body_returns_413(self):
        # The body cap fires before splitting, so any >5 MB content works.
        big = "x" * (5 * 1024 * 1024 + 10)
        r = self._post(big)
        assert r.status_code == 413
        assert "exceeds" in r.json()["detail"]

    def test_charge_out_of_range_rejected(self, water_xyz):
        r = self._post(water_xyz, charge=999)
        assert r.status_code == 422

    def test_crlf_handled(self, water_xyz):
        r = self._post(water_xyz.replace("\n", "\r\n"))
        assert r.status_code == 200
        assert r.json()["structures"][0]["inchikey"] == WATER_INCHIKEY

    # --- toolkit=openbabel ---

    def test_water_openbabel(self, water_xyz):
        r = self._post(water_xyz, toolkit="openbabel")
        assert r.status_code == 200
        s = r.json()["structures"][0]
        assert s["method"] == "bond_orders"  # OB always reports bond_orders
        assert s["bond_orders_perceived"] is True
        assert "\t" not in s["canonicalsmiles"]

    # --- SDF roundtrip ---

    def test_sdf_roundtrip(self, water_xyz):
        r = self._post(water_xyz + water_xyz)
        body = r.json()
        suppl = Chem.SDMolSupplier()
        suppl.SetData(body["sdf"])
        mols = [m for m in suppl if m is not None]
        assert len(mols) == body["summary"]["successful"]

    # --- title preserved in MOL block ---

    def test_title_preserved_in_molblock(self):
        xyz = "3\nMyMolecule\nO 0 0 0\nH 0.76 0.59 0\nH -0.76 0.59 0\n"
        r = self._post(xyz)
        body = r.json()
        assert "MyMolecule" in body["structures"][0]["molblock"]
