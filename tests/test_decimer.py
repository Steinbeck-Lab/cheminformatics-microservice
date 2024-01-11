from __future__ import annotations

import os

import pytest

from app.modules.decimer import convert_image
from app.modules.decimer import get_predicted_segments
from app.modules.decimer import get_predicted_segments_from_file
from app.modules.decimer import get_segments


# Define a directory for temporary test files
TEST_FILES_DIR = "tests"


@pytest.fixture(scope="module")
def sample_gif_path():
    return os.path.join(TEST_FILES_DIR, "segment_sample.gif")


@pytest.fixture(scope="module")
def sample_png_path():
    return os.path.join(TEST_FILES_DIR, "segment_sample.png")


@pytest.fixture(scope="module")
def sample_image_path():
    return os.path.join(TEST_FILES_DIR, "segment_sample.png")


# Test the convert_image function
def test_convert_image(sample_gif_path, sample_png_path):
    converted_path = convert_image(sample_gif_path)
    assert os.path.isfile(converted_path)
    assert converted_path == sample_png_path


# Test the get_segments function
def test_get_segments(sample_gif_path):
    image_name, segments = get_segments(sample_gif_path)
    assert image_name == "segment_sample.gif"
    assert len(segments) > 0


# Test the get_predicted_segments function
def test_get_predicted_segments(sample_gif_path):
    predicted_smiles = get_predicted_segments(sample_gif_path)
    assert isinstance(predicted_smiles, str)
    assert len(predicted_smiles) > 0


# Test the get_predicted_segments_from_file function
def test_get_predicted_segments_from_file(sample_image_path):
    with open(sample_image_path, "rb") as f:
        content = f.read()
    predicted_smiles = get_predicted_segments_from_file(
        content,
        "caffeine.png",
    )
    assert isinstance(predicted_smiles, str)
    assert len(predicted_smiles) > 0
