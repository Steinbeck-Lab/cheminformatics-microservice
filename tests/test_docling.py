from __future__ import annotations

from importlib.util import find_spec

import pytest
from fastapi.testclient import TestClient

from app.modules.docling_wrapper import combine_to_paragraph
from app.modules.docling_wrapper import extract_paper_content
from app.modules.docling_wrapper import extract_text
from app.modules.docling_wrapper import first_page_text
from app.routers.pdf import MAX_PDF_SIZE


def _text(label, text, page=1):
    return {"label": label, "text": text, "prov": [{"page_no": page}]}


# Same shape as DoclingDocument.export_to_dict() from docling 2.x
DOC = {
    "schema_name": "DoclingDocument",
    "texts": [
        _text("page_header", "J. Chem. Inf. 2026"),
        _text("section_header", "Deep Learning for Chemical Structure Recognition"),
        _text("text", "ABSTRACT: We convert images of molecules into SMILES."),
        _text("section_header", "INTRODUCTION"),
        _text("text", "Structures in papers are often only images."),
        _text("section_header", "RESULTS AND DISCUSSION"),
        _text("text", "The model reaches 95 percent accuracy."),
        _text("text", "Text on the second page.", page=2),
    ],
}


def test_extract_paper_content():
    content = extract_paper_content(DOC)
    assert content == {
        "title": "Deep Learning for Chemical Structure Recognition",
        "abstract": "We convert images of molecules into SMILES.",
        "main_text": "Structures in papers are often only images.",
    }
    assert combine_to_paragraph(content) == (
        "Deep Learning for Chemical Structure Recognition"
        " We convert images of molecules into SMILES."
        " Structures in papers are often only images."
    )


def test_title_label_wins_over_section_header():
    doc = {
        "texts": [
            _text("section_header", "Long section heading"),
            _text("title", "Real Title"),
        ]
    }
    assert extract_paper_content(doc)["title"] == "Real Title"


def test_abstract_after_bare_heading():
    doc = {
        "texts": [
            _text("section_header", "Abstract"),
            _text("text", "The abstract body."),
        ]
    }
    assert extract_paper_content(doc)["abstract"] == "The abstract body."


def test_first_page_text_skips_headers_and_other_pages():
    text = first_page_text(DOC)
    assert "J. Chem. Inf." not in text
    assert "second page" not in text
    assert text.startswith("Deep Learning")


def test_extract_text_falls_back_to_first_page():
    doc = {"texts": [_text("text", "Only a short note.")]}
    assert extract_text(doc) == "Only a short note."


docling_routes = pytest.mark.skipif(
    not find_spec("docling"),
    reason="docling routes are only registered when docling is installed",
)


@pytest.fixture
def client(monkeypatch):
    from app.main import app
    from app.routers import pdf

    # Skip the real model run; the wrapper functions are tested above.
    monkeypatch.setattr(
        pdf, "get_converted_document", lambda content, number_of_pages: DOC
    )
    return TestClient(app)


@docling_routes
def test_extract_text(client):
    response = client.post(
        "/latest/tools/extract_text",
        files={"pdf_file": ("my paper.pdf", b"%PDF-1.7 fake", "application/pdf")},
    )
    assert response.status_code == 200
    assert response.json()["pdf_filename"] == "my_paper.pdf"
    assert response.json()["text"].startswith("Deep Learning")


@docling_routes
def test_extract_json(client):
    response = client.post(
        "/latest/tools/extract_json",
        files={"pdf_file": ("paper.pdf", b"%PDF-1.7 fake", "application/pdf")},
        data={"pages": "2"},
    )
    assert response.status_code == 200
    assert response.json()["schema_name"] == "DoclingDocument"


@docling_routes
@pytest.mark.parametrize(
    "content, code",
    [(b"not a pdf", 400), (b"%PDF-" + b"0" * MAX_PDF_SIZE, 413)],
)
def test_extract_text_rejects_bad_uploads(client, content, code):
    response = client.post(
        "/latest/tools/extract_text",
        files={"pdf_file": ("paper.pdf", content, "application/pdf")},
    )
    assert response.status_code == code
