from __future__ import annotations

import io
import re
import threading
from functools import lru_cache

SKIP_LABELS = {"page_header", "page_footer"}
STOP_HEADERS = ("RESULTS", "DISCUSSION")
CONVERT_TIMEOUT = 120  # seconds Docling may spend on one document
LOCK_TIMEOUT = 30  # seconds a request waits for the converter before giving up

# ponytail: one converter behind a global lock, a converter pool if PDF throughput matters
_convert_lock = threading.Lock()


class ConverterBusyError(RuntimeError):
    """Raised when the converter stays busy for longer than LOCK_TIMEOUT."""


@lru_cache(maxsize=1)
def _get_converter():
    # Imported lazily: docling pulls in torch, and the lite image ships without it.
    from docling.datamodel.base_models import InputFormat
    from docling.datamodel.pipeline_options import PdfPipelineOptions
    from docling.document_converter import DocumentConverter
    from docling.document_converter import PdfFormatOption

    options = PdfPipelineOptions(document_timeout=CONVERT_TIMEOUT)
    return DocumentConverter(
        format_options={InputFormat.PDF: PdfFormatOption(pipeline_options=options)}
    )


def get_converted_document(pdf_bytes: bytes, number_of_pages: int = 3) -> dict:
    """Convert the first pages of a PDF into a Docling document dictionary.

    Args:
        pdf_bytes (bytes): Raw content of the PDF file.
        number_of_pages (int): Number of pages to convert, counted from the first page.

    Returns:
        dict: The DoclingDocument exported as a dictionary.
    """
    from docling.datamodel.base_models import DocumentStream

    source = DocumentStream(name="upload.pdf", stream=io.BytesIO(pdf_bytes))
    if not _convert_lock.acquire(timeout=LOCK_TIMEOUT):
        raise ConverterBusyError("PDF converter is busy")
    try:
        result = _get_converter().convert(source, page_range=(1, number_of_pages))
    finally:
        _convert_lock.release()
    return result.document.export_to_dict()


def _clean(text: str) -> str:
    text = re.sub(r"\s+", " ", text)
    return re.sub(r"\s+([.,;:?!])", r"\1", text).strip()


def _body_texts(doc: dict) -> list:
    return [
        t
        for t in doc.get("texts", [])
        if t.get("text", "").strip() and t.get("label") not in SKIP_LABELS
    ]


def extract_paper_content(doc: dict) -> dict:
    """Extract the title, abstract and the text before the results section.

    Args:
        doc (dict): A Docling document dictionary.

    Returns:
        dict: The keys "title", "abstract" and "main_text".
    """
    texts = _body_texts(doc)
    headers = [t["text"] for t in texts if t.get("label") == "section_header"]

    title = next((t["text"] for t in texts if t.get("label") == "title"), "") or next(
        (
            h
            for h in headers
            if len(h) > 10 and not any(s in h.upper() for s in STOP_HEADERS)
        ),
        "",
    )

    # Texts after the first one that mentions "abstract", up to the results section
    start = next(
        (i for i, t in enumerate(texts) if "abstract" in t["text"].lower()), None
    )
    abstract, main_text = "", []
    if start is not None:
        abstract = re.sub(r"^\s*abstract\s*[:.]?", "", texts[start]["text"], flags=re.I)
        for t in texts[start + 1 :]:
            if t.get("label") != "section_header":
                main_text.append(t["text"])
            elif any(s in t["text"].upper() for s in STOP_HEADERS):
                break
        if not abstract.strip() and main_text:
            # A bare "Abstract" heading: the abstract is the next paragraph.
            abstract = main_text.pop(0)

    return {
        "title": title,
        "abstract": abstract.strip(),
        "main_text": " ".join(main_text),
    }


def combine_to_paragraph(content: dict) -> str:
    """Join the title, abstract and main text into one clean paragraph."""
    return _clean(
        " ".join(content.get(k, "") for k in ("title", "abstract", "main_text"))
    )


def first_page_text(doc: dict) -> str:
    """Return all body text found on the first page of a Docling document."""
    return _clean(
        " ".join(
            t["text"]
            for t in _body_texts(doc)
            if any(p.get("page_no") == 1 for p in t.get("prov", []))
        )
    )


def extract_text(doc: dict) -> str:
    """Return the title, abstract and introduction as one paragraph.

    Falls back to all text on the first page when 10 words or fewer are found.
    """
    text = combine_to_paragraph(extract_paper_content(doc))
    return text if len(text.split()) > 10 else first_page_text(doc)
