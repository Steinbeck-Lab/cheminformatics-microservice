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
        "/latest/ocsr/process", data=input, headers={"Content-Type": "application/json"}
    )
    assert response.status_code == response_code


# Run the tests
if __name__ == "__main__":
    pytest.main()
