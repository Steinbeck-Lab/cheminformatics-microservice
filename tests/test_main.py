import os
from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)


def test_root_redirect():
    response = client.get("/latest/docs")
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/html; charset=utf-8"


def test_root_endpoint():
    response = client.get("/")
    assert response.status_code == 200
    expected_url = "http://testserver" + os.getenv("HOMEPAGE_URL", "/latest/docs")
    assert response.url == expected_url


def test_get_health():
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}
