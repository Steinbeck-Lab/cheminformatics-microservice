import os
import re

from fastapi import HTTPException, status, File, Form, UploadFile
from docling.datamodel.base_models import InputFormat
from docling.datamodel.pipeline_options import PdfPipelineOptions
from docling.document_converter import (
    DocumentConverter,
    PdfFormatOption,
    ConversionResult,
)
from PyPDF2 import PdfReader, PdfWriter
from app.config import PDF_DIR

artifacts_path = ".cache/docling/models"

pipeline_options = PdfPipelineOptions(artifacts_path=artifacts_path)
doc_converter = DocumentConverter(
    format_options={InputFormat.PDF: PdfFormatOption(pipeline_options=pipeline_options)}
)


def extract_first_three_pages(input_pdf_path, number_of_pages=3):
    """
    Extract a specified number of pages from the beginning of a PDF file.

    Creates a new PDF file containing only the first N pages from the input PDF,
    useful for processing large documents where only the initial pages are needed.

    Args:
        input_pdf_path (str): Path to the source PDF file to extract pages from.
        number_of_pages (int, optional): Number of pages to extract from the beginning. Defaults to 3.

    Returns:
        str: Path to the newly created PDF file containing the extracted pages.
             The output file has "_out" appended before the file extension.

    Raises:
        FileNotFoundError: If the input PDF file does not exist.
        Exception: If PDF reading/writing operations fail.

    Example:
        >>> output_path = extract_first_three_pages("document.pdf", 2)
        >>> print(output_path)
        "document_out.pdf"
    """
    # Open the original PDF
    with open(input_pdf_path, "rb") as input_pdf_file:
        pdf_reader = PdfReader(input_pdf_file)
        pdf_writer = PdfWriter()

        # Determine the number of pages to extract
        num_pages = min(number_of_pages, len(pdf_reader.pages))

        # Add the pages to the writer object
        for page_num in range(num_pages):
            page = pdf_reader.pages[page_num]
            pdf_writer.add_page(page)

        # Construct the output file path
        base_name = os.path.splitext(input_pdf_path)[0]
        output_pdf_file_path = f"{base_name}_out.pdf"

        # Write the pages to a new PDF
        with open(output_pdf_file_path, "wb") as output_pdf_file:
            pdf_writer.write(output_pdf_file)
        print(output_pdf_file_path)
        return output_pdf_file_path


def get_converted_document(path, number_of_pages=3):
    """
    Convert a PDF document to structured JSON format using Docling.

    Extracts the first N pages from a PDF and converts them to a structured
    document format that can be processed for content extraction.

    Args:
        path (str): Path to the PDF file to be converted.
        number_of_pages (int, optional): Number of pages to process. Defaults to 3.

    Returns:
        dict: Document structure in JSON format containing text elements,
              layout information, and metadata from the converted PDF.

    Raises:
        Exception: If document conversion fails or file cannot be processed.

    Example:
        >>> doc_dict = get_converted_document("paper.pdf", 2)
        >>> print(doc_dict.keys())
        dict_keys(['texts', 'schema_name', 'name', ...])
    """
    output_pdf_file_path = extract_first_three_pages(path, number_of_pages)
    converter = DocumentConverter()
    conv_result: ConversionResult = converter.convert(output_pdf_file_path)
    conv_result_dict = conv_result.document.export_to_dict()
    return conv_result_dict


def extract_paper_content(doc_json):
    """
    Extract the title, abstract, and main text up to the results section from a document JSON.

    Args:
        doc_json (dict): The JSON representation of the document

    Returns:
        dict: A dictionary containing the title, abstract, and main text
    """
    title = ""
    abstract = ""
    main_text = []

    # Get all text elements
    texts = doc_json.get("texts", [])

    # Find the title (usually the first section_header with level 1)
    for text in texts:
        if text.get("label") == "section_header" and text.get("level") == 1:
            if (
                "RESULTS" not in text.get("text", "").upper()
                and len(text.get("text", "")) > 10
            ):
                title = text.get("text", "")
                break

    # Find the abstract (usually a caption or text containing "ABSTRACT:")
    for text in texts:
        if "ABSTRACT:" in text.get("text", ""):
            abstract = text.get("text", "").replace("ABSTRACT:", "").strip()
            break

    # Extract main text up to results section
    found_abstract = False
    for text in texts:
        content = text.get("text", "")

        # Skip empty content
        if not content.strip():
            continue

        # Skip titles, headers, page numbers, etc.
        if text.get("label") in ["page_header", "page_footer"]:
            continue

        # Start collecting text after abstract
        if not found_abstract and "ABSTRACT:" in content:
            found_abstract = True
            continue

        # Stop at results section
        if found_abstract and (
            text.get("label") == "section_header"
            and ("RESULTS" in content.upper() or "DISCUSSION" in content.upper())
        ):
            break

        # Collect main text
        if found_abstract:
            main_text.append(content)

    # Join the main text paragraphs
    main_text_str = " ".join(main_text)

    return {"title": title, "abstract": abstract, "main_text": main_text_str}


def extract_from_docling_document(data):
    """
    Extract paper content from a Docling document format.

    Processes a structured Docling document to extract key academic paper components
    including title, abstract, and main text content.

    Args:
        data (dict): The Docling document JSON containing structured document data.

    Returns:
        dict: Dictionary containing extracted content with keys:
              - title: Paper title text
              - abstract: Abstract content
              - main_text: Main body text
              Or error message if format is invalid.

    Example:
        >>> content = extract_from_docling_document(docling_json)
        >>> print(content['title'])
        "A Novel Approach to Chemical Structure Recognition"
    """
    # For Docling format, we need to extract the main document structure
    if "schema_name" in data and data["schema_name"] == "DoclingDocument":
        result = extract_paper_content(data)

        # Check if any of the key components are empty and extract additional info if needed
        if (
            not result.get("title")
            or not result.get("abstract")
            or not result.get("main_text")
        ):
            # Try to extract more text based on document structure
            texts = data.get("texts", [])

            # If title is empty, try to find a likely title
            if not result.get("title"):
                for text in texts:
                    # Look for large font text at the beginning
                    if (
                        text.get("label") == "paragraph"
                        and text.get("page_number") == 1
                        and text.get("font_size", 0) > 12
                    ):  # Assuming larger font indicates title
                        result["title"] = text.get("text", "")
                        break

            # If abstract is empty, look for any text with "abstract" in it
            if not result.get("abstract"):
                for text in texts:
                    content = text.get("text", "").lower()
                    if "abstract" in content:
                        result["abstract"] = text.get("text", "")
                        break

        return result
    else:
        return {"error": "Not a valid Docling document format"}


def combine_to_paragraph(result_dict):
    """
    Combine extracted paper components into a single formatted paragraph.

    Merges title, abstract, and main text from paper extraction results into
    a clean, properly formatted single paragraph with normalized spacing.

    Args:
        result_dict (dict): Dictionary containing 'title', 'abstract', and 'main_text' keys.

    Returns:
        str: Single paragraph containing all paper content with cleaned formatting.
             Returns error message if input is invalid.

    Example:
        >>> content = {"title": "Paper Title", "abstract": "Abstract text", "main_text": "Body text"}
        >>> paragraph = combine_to_paragraph(content)
        >>> print(len(paragraph.split()))
        150
    """
    # Check for valid input
    if not isinstance(result_dict, dict):
        return "Error: Input must be a dictionary."

    # Extract components
    title = result_dict.get("title", "").strip()
    abstract = result_dict.get("abstract", "").strip()
    main_text = result_dict.get("main_text", "").strip()

    # Combine components with appropriate spacing
    combined_text = []

    if title:
        combined_text.append(title)

    if abstract:
        combined_text.append(abstract)

    if main_text:
        combined_text.append(main_text)

    # Join with spaces and clean up any double spaces
    combined_paragraph = " ".join(combined_text)

    combined_paragraph = re.sub(r"\s+", " ", combined_paragraph)
    combined_paragraph = re.sub(r"\s+([.,;:?!])", r"\1", combined_paragraph)

    return combined_paragraph


async def extract_pdf_text(
    pdf_file: UploadFile = File(...),
    pages: int = Form(1, description="Number of pages to process"),
):
    """
    Extract and process text content from an uploaded PDF file.

    Handles PDF upload, converts to structured format, extracts paper content,
    and returns combined text. Falls back to full page extraction if content is insufficient.

    Args:
        pdf_file (UploadFile): The PDF file to process (required).
        pages (int): Number of pages to process from the beginning. Defaults to 1.

    Returns:
        dict: JSON object containing:
              - text: Combined extracted text content
              - pdf_filename: Sanitized filename of the processed PDF

    Raises:
        HTTPException:
            - 400: If uploaded file is not a PDF
            - 500: If PDF processing fails

    Example:
        >>> result = await extract_pdf_text(pdf_file, pages=2)
        >>> print(result['text'][:100])
        "Title: Novel Chemical Analysis Abstract: This paper presents..."
    """
    if not pdf_file.filename.lower().endswith(".pdf"):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Uploaded file must be a PDF",
        )

    try:
        # Use original filename without adding unique IDs
        original_filename = pdf_file.filename
        safe_filename = original_filename.replace(" ", "_")

        # Create the full file path
        file_path = os.path.join(PDF_DIR, safe_filename)

        # Check if file already exists
        if os.path.exists(file_path):
            # Use existing file
            pass
        else:
            # Save the uploaded file with original name
            with open(file_path, "wb") as buffer:
                content = await pdf_file.read()
                buffer.write(content)

        # Process the PDF file
        json_data = get_converted_document(file_path, number_of_pages=pages)
        result = extract_from_docling_document(json_data)
        combined_text = combine_to_paragraph(result)

        # Check if the extracted text has more than 10 words
        if len(combined_text.split()) <= 10:
            # If extracted text has 10 or fewer words, extract the whole first page
            all_text = extract_full_page_text(json_data)
            if all_text:
                combined_text = all_text

        return {"text": combined_text, "pdf_filename": safe_filename}

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error processing PDF: {str(e)}",
        )


def extract_full_page_text(doc_json):
    """
    Extract all text content from the first page of a document.

    Retrieves and concatenates all text elements from the first page,
    excluding headers and footers, with cleaned formatting.

    Args:
        doc_json (dict): The JSON representation of the document structure.

    Returns:
        str: Concatenated and cleaned text from the first page.
             Falls back to all pages if page numbers are unavailable.

    Example:
        >>> text = extract_full_page_text(document_json)
        >>> print(len(text.split()))
        245
    """
    # Get all text elements
    texts = doc_json.get("texts", [])

    # Filter texts from the first page
    page_1_texts = []
    for text in texts:
        # Check if this text element is on page 1
        if text.get("page_number") == 1:  # Assuming page numbers start at 1
            content = text.get("text", "").strip()
            if content and text.get("label") not in ["page_header", "page_footer"]:
                page_1_texts.append(content)

    # If there's no specific page number field, try to get all texts
    # This is a fallback in case the document structure doesn't have page numbers
    if not page_1_texts and texts:
        # Just take all texts as a fallback
        for text in texts:
            content = text.get("text", "").strip()
            if content and text.get("label") not in ["page_header", "page_footer"]:
                page_1_texts.append(content)

    # Join all text elements with spaces
    full_text = " ".join(page_1_texts)

    # Clean up the text
    full_text = re.sub(r"\s+", " ", full_text)
    full_text = re.sub(r"\s+([.,;:?!])", r"\1", full_text)

    return full_text
