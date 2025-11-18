from __future__ import annotations

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


class TestEnhancedDepictEndpoint:
    """Test the enhanced depiction endpoint /latest/depict/2D_enhanced."""

    def test_basic_molecule(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO")
        assert response.status_code == 200
        assert "svg" in response.text.lower()

    def test_caffeine(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        )
        assert response.status_code == 200
        assert "svg" in response.text.lower()

    def test_invalid_smiles(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=INVALID_SMILES{]")
        assert response.status_code == 422


class TestCXSMILESHighlighting:
    """Test CXSMILES highlighting functionality."""

    def test_atom_highlighting(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=c1ccccc1 |ha:0,1,2|")
        assert response.status_code == 200
        assert "svg" in response.text.lower()

    def test_bond_highlighting(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO |hb:0,1|")
        assert response.status_code == 200
        assert "svg" in response.text.lower()

    def test_combined_highlighting(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1O |ha:0,1,2,hb:0,1|"
        )
        assert response.status_code == 200
        assert "svg" in response.text.lower()


class TestAbbreviations:
    """Test abbreviation modes."""

    def test_abbreviate_off(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=C1CCOC1&abbreviate=off"
        )
        assert response.status_code == 200

    def test_abbreviate_groups(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1C&abbreviate=groups"
        )
        assert response.status_code == 200

    def test_abbreviate_reagents(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=C1CCOC1&abbreviate=reagents"
        )
        assert response.status_code == 200

    def test_abbreviate_all(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1CC(C)(C)OC(=O)N&abbreviate=on"
        )
        assert response.status_code == 200


class TestDativeBonds:
    """Test dative bond perception modes."""

    def test_dative_never(self):
        # Use a simple molecule instead of invalid [NH3]BF3
        response = client.get("/latest/depict/2D_enhanced?smiles=CS(=O)C&dative=never")
        assert response.status_code == 200

    def test_dative_metals(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=[Co][N+]([O-])(=O)&dative=metals"
        )
        assert response.status_code == 200

    def test_dative_always(self):
        # Use dimethyl sulfoxide instead of invalid [NH3]BF3
        response = client.get("/latest/depict/2D_enhanced?smiles=CS(=O)C&dative=always")
        assert response.status_code == 200


class TestMulticenterBonds:
    """Test multicenter bond display styles."""

    def test_multicenter_provided(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=[Fe]c1ccccc1&multicenter=provided"
        )
        assert response.status_code == 200

    def test_multicenter_dative(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=[Fe]c1ccccc1&multicenter=dative"
        )
        assert response.status_code == 200

    def test_multicenter_dashed(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=[Fe]c1ccccc1&multicenter=dashed"
        )
        assert response.status_code == 200

    def test_multicenter_hidden(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=[Fe]c1ccccc1&multicenter=hidden"
        )
        assert response.status_code == 200


class TestAnnotationModes:
    """Test annotation modes."""

    def test_annotate_none(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&annotate=none")
        assert response.status_code == 200

    def test_annotate_number(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&annotate=number")
        assert response.status_code == 200

    def test_annotate_bondnumber(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&annotate=bondnumber"
        )
        assert response.status_code == 200

    def test_annotate_mapidx(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&annotate=mapidx")
        assert response.status_code == 200

    def test_annotate_atomvalue(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&annotate=atomvalue"
        )
        assert response.status_code == 200

    def test_annotate_colmap(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&annotate=colmap")
        assert response.status_code == 200

    def test_annotate_cip(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=C[C@H](O)CC&annotate=cip"
        )
        assert response.status_code == 200


class TestStylePresets:
    """Test style preset modes."""

    def test_style_cow(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&style=cow")
        assert response.status_code == 200

    def test_style_cob(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&style=cob")
        assert response.status_code == 200

    def test_style_bow(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&style=bow")
        assert response.status_code == 200

    def test_style_nob(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&style=nob")
        assert response.status_code == 200


class TestAromaticDisplay:
    """Test aromatic donut display."""

    def test_donuts_enabled(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=c1ccccc1&donuts=true")
        assert response.status_code == 200

    def test_donuts_disabled(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=c1ccccc1&donuts=false")
        assert response.status_code == 200


class TestAdvancedControls:
    """Test advanced rendering controls."""

    def test_zoom_parameter(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&zoom=1.5")
        assert response.status_code == 200

    def test_ratio_parameter(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&ratio=1.2")
        assert response.status_code == 200

    def test_flip_parameter(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&flip=true")
        assert response.status_code == 200

    def test_anon_parameter(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&anon=true")
        assert response.status_code == 200

    def test_svgunits_px(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&svgunits=px")
        assert response.status_code == 200

    def test_svgunits_mm(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&svgunits=mm")
        assert response.status_code == 200


class TestHydrogenDisplay:
    """Test hydrogen display modes."""

    def test_hydrogen_provided(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&hydrogen_display=Provided"
        )
        assert response.status_code == 200

    def test_hydrogen_minimal(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&hydrogen_display=Minimal"
        )
        assert response.status_code == 200

    def test_hydrogen_explicit(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&hydrogen_display=Explicit"
        )
        assert response.status_code == 200

    def test_hydrogen_stereo(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=C[C@H](O)CC&hydrogen_display=Stereo"
        )
        assert response.status_code == 200

    def test_hydrogen_smart(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=C[C@H](O)CC&hydrogen_display=Smart"
        )
        assert response.status_code == 200


class TestCustomDimensions:
    """Test custom image dimensions."""

    def test_custom_width_height(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&width=1024&height=1024"
        )
        assert response.status_code == 200

    def test_small_dimensions(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&width=256&height=256"
        )
        assert response.status_code == 200


class TestRotation:
    """Test molecule rotation."""

    def test_rotate_90(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&rotate=90")
        assert response.status_code == 200

    def test_rotate_180(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&rotate=180")
        assert response.status_code == 200

    def test_rotate_270(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&rotate=270")
        assert response.status_code == 200


class TestSMARTSHighlighting:
    """Test SMARTS pattern highlighting."""

    def test_smarts_highlight_carbonyl(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&highlight=C=O"
        )
        assert response.status_code == 200

    def test_smarts_highlight_alcohol(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&highlight=CO")
        assert response.status_code == 200


class TestAtomHighlighting:
    """Test atom index highlighting."""

    def test_highlight_single_atom(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&atomIds=0")
        assert response.status_code == 200

    def test_highlight_multiple_atoms(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&atomIds=0,1,2")
        assert response.status_code == 200

    def test_highlight_invalid_format(self):
        # Invalid atomIds are silently filtered out, not rejected
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&atomIds=invalid")
        assert response.status_code == 200


class TestCustomColors:
    """Test custom background and foreground colors."""

    def test_custom_bgcolor(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&bgcolor=%23FFFFFF")
        assert response.status_code == 200

    def test_custom_fgcolor(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO&fgcolor=%23000000")
        assert response.status_code == 200

    def test_both_custom_colors(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=CCO&bgcolor=%23FFFFFF&fgcolor=%23000000"
        )
        assert response.status_code == 200


class TestChiralMolecules:
    """Test chiral molecule depiction."""

    def test_chiral_with_cip(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=C[C@H](O)CC&CIP=true")
        assert response.status_code == 200

    def test_chiral_without_cip(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=C[C@H](O)CC&CIP=false")
        assert response.status_code == 200

    def test_multiple_chiral_centers(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=C[C@H](O)[C@@H](C)Cl&CIP=true"
        )
        assert response.status_code == 200


class TestComplexMolecules:
    """Test complex molecule structures."""

    def test_steroid_structure(self):
        smiles = "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"
        response = client.get(f"/latest/depict/2D_enhanced?smiles={smiles}")
        assert response.status_code == 200

    def test_peptide_with_protecting_groups(self):
        smiles = "CC(C)(C)OC(=O)NCC(=O)NCC(=O)O"
        response = client.get(f"/latest/depict/2D_enhanced?smiles={smiles}")
        assert response.status_code == 200

    def test_organometallic(self):
        smiles = "[Fe]c1ccccc1"
        response = client.get(f"/latest/depict/2D_enhanced?smiles={smiles}")
        assert response.status_code == 200


class TestCombinedFeatures:
    """Test combinations of multiple features."""

    def test_abbreviations_with_annotations(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1C&abbreviate=groups&annotate=number"
        )
        assert response.status_code == 200

    def test_dative_with_style(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=[Co][N+]([O-])(=O)&dative=metals&style=bow"
        )
        assert response.status_code == 200

    def test_donuts_with_highlighting(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1&donuts=true&atomIds=0,1,2"
        )
        assert response.status_code == 200

    def test_all_features_combined(self):
        response = client.get(
            "/latest/depict/2D_enhanced?smiles=c1ccccc1C&abbreviate=groups&"
            "annotate=number&style=bow&donuts=true&zoom=1.5&ratio=1.2"
        )
        assert response.status_code == 200


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_single_atom(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=C")
        assert response.status_code == 200

    def test_disconnected_fragments(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=CCO.C1CCOC1")
        assert response.status_code == 200

    def test_charged_species(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=[NH4+]")
        assert response.status_code == 200

    def test_radical_species(self):
        response = client.get("/latest/depict/2D_enhanced?smiles=[CH3]")
        assert response.status_code == 200
