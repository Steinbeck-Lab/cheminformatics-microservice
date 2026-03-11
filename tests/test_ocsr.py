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
