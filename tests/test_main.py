import pytest
from fastapi.testclient import TestClient
from app.main import app
from app.main import HealthCheck

client = TestClient(app)


def test_root_redirect():
    response = client.get("/latest/docs")
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/html; charset=utf-8"


def test_get_health():
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json() == {"status": "OK"}
