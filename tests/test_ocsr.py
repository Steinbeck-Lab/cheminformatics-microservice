from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_ocsr_index():
    response = client.get("/latest/ocsr/")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}


@pytest.mark.parametrize(
    "input, response_code",
    [
        (
            '{"path": "https://static-content.springer.com/image/art%3A10.1186%2Fs13321-023-00744-6/MediaObjects/13321_2023_744_Figa_HTML.png","reference": "xyzzy"}',
            200,
        ),
        ("INVALID_INPUT", 422),
    ],
)
def test_process(input, response_code):
    response = client.post(
        "/latest/ocsr/process",
        data=input,
        headers={"Content-Type": "application/json"},
    )
    assert response.status_code == response_code


def test_process_upload():
    file_path = "tests/caffeine.png"
    file_content = open(file_path, "rb").read()

    response = client.post(
        "/latest/ocsr/process-upload",
        files={"file": ("caffeine.png", file_content, "image/png")},
    )
    assert response.status_code == 200


# ---------------------------------------------------------------------------
# _validate_url SSRF protection tests
# ---------------------------------------------------------------------------


class TestValidateUrl:
    """Tests for the _validate_url SSRF protection helper."""

    def test_valid_https_url(self):
        from app.routers.ocsr import _validate_url

        _validate_url("https://example.com/image.png")  # should not raise

    def test_valid_http_url(self):
        from app.routers.ocsr import _validate_url

        _validate_url("http://example.com/image.png")  # should not raise

    def test_rejects_ftp_scheme(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="scheme.*not allowed"):
            _validate_url("ftp://example.com/file")

    def test_rejects_file_scheme(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="scheme.*not allowed"):
            _validate_url("file:///etc/passwd")

    def test_rejects_no_hostname(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="no hostname"):
            _validate_url("https:///path/only")

    def test_rejects_private_ip(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="internal/private"):
            _validate_url("https://192.168.1.1/image.png")

    def test_rejects_loopback_ip(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="internal/private"):
            _validate_url("https://127.0.0.1/image.png")

    def test_rejects_link_local_ip(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="internal/private"):
            _validate_url("https://169.254.169.254/metadata")

    def test_rejects_localhost(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="internal/private"):
            _validate_url("https://localhost/secret")

    def test_rejects_dot_local(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="internal/private"):
            _validate_url("https://myhost.local/secret")

    def test_rejects_10_x_private(self):
        from app.routers.ocsr import _validate_url

        with pytest.raises(ValueError, match="internal/private"):
            _validate_url("https://10.0.0.1/admin")


# ---------------------------------------------------------------------------
# OCSR endpoint SSRF rejection via API
# ---------------------------------------------------------------------------


def test_process_rejects_private_url():
    """OCSR /process should reject private/internal URLs."""
    payload = '{"path": "https://127.0.0.1/evil.png", "reference": "test"}'
    response = client.post(
        "/latest/ocsr/process",
        data=payload,
        headers={"Content-Type": "application/json"},
    )
    assert response.status_code == 422


def test_process_rejects_localhost_url():
    """OCSR /process should reject localhost URLs."""
    payload = '{"path": "https://localhost/evil.png", "reference": "test"}'
    response = client.post(
        "/latest/ocsr/process",
        data=payload,
        headers={"Content-Type": "application/json"},
    )
    assert response.status_code == 422


# ---------------------------------------------------------------------------
# Additional coverage tests for ocsr.py uncovered lines
# ---------------------------------------------------------------------------


class TestOCSRProcessImgParameter:
    """Test /process endpoint with img parameter (lines 154-170)."""

    def test_process_with_img_invalid_url(self):
        """Test img parameter with invalid/private URL (line 154-170)."""
        payload = {
            "img": "https://127.0.0.1/image.png",
            "reference": "test",
        }
        response = client.post(
            "/latest/ocsr/process",
            json=payload,
        )
        assert response.status_code == 422

    def test_process_with_img_nonexistent_url(self):
        """Test img parameter with URL that doesn't exist (lines 168-173)."""
        payload = {
            "img": "https://example.com/nonexistent_image_xyz.png",
            "reference": "test",
        }
        response = client.post(
            "/latest/ocsr/process",
            json=payload,
        )
        assert response.status_code == 422

    def test_process_with_img_ftp_scheme(self):
        """Test img parameter with disallowed scheme."""
        payload = {
            "img": "ftp://example.com/image.png",
            "reference": "test",
        }
        response = client.post(
            "/latest/ocsr/process",
            json=payload,
        )
        assert response.status_code == 422


class TestOCSRProcessPathErrors:
    """Test /process endpoint path parameter error handling (lines 189, 196-198)."""

    def test_process_with_path_nonexistent_url(self):
        """Test path parameter with URL that returns non-200 (line 189)."""
        payload = {
            "path": "https://httpstat.us/404",
            "reference": "test",
        }
        response = client.post(
            "/latest/ocsr/process",
            json=payload,
        )
        # Should fail with 422 since the URL returns 404
        assert response.status_code == 422

    def test_process_with_path_invalid_url_scheme(self):
        """Test path parameter with disallowed URL scheme (line 196-198)."""
        payload = {
            "path": "ftp://example.com/file.png",
            "reference": "test",
        }
        response = client.post(
            "/latest/ocsr/process",
            json=payload,
        )
        assert response.status_code == 422

    def test_process_with_path_private_ip(self):
        """Test path parameter with private IP address."""
        payload = {
            "path": "https://10.0.0.1/image.png",
            "reference": "test",
        }
        response = client.post(
            "/latest/ocsr/process",
            json=payload,
        )
        assert response.status_code == 422

    def test_process_with_no_path_no_img(self):
        """Test /process with neither path nor img."""
        payload = {
            "reference": "test",
        }
        response = client.post(
            "/latest/ocsr/process",
            json=payload,
        )
        assert response.status_code == 422


class TestOCSRProcessUploadErrors:
    """Test /process-upload endpoint error handling (lines 248-258)."""

    def test_process_upload_invalid_image_content(self):
        """Test upload with invalid image content (line 248-253)."""
        response = client.post(
            "/latest/ocsr/process-upload",
            files={"file": ("test.png", b"not-an-image", "image/png")},
        )
        assert response.status_code == 422

    def test_process_upload_empty_file(self):
        """Test upload with empty file content."""
        response = client.post(
            "/latest/ocsr/process-upload",
            files={"file": ("test.png", b"", "image/png")},
        )
        assert response.status_code == 422

    def test_process_upload_text_as_image(self):
        """Test upload with text content instead of image."""
        response = client.post(
            "/latest/ocsr/process-upload",
            files={"file": ("test.txt", b"This is not an image", "text/plain")},
        )
        assert response.status_code == 422
