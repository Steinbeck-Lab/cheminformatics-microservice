from __future__ import annotations

import os
import tempfile
from unittest.mock import patch

import pytest
from PIL import Image

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


@pytest.fixture(scope="module")
def small_image_path():
    """Small image (400x300) - should trigger direct prediction"""
    return os.path.join(TEST_FILES_DIR, "small_molecule.png")


@pytest.fixture(scope="module")
def tiny_image_path():
    """Tiny image (200x150) - should trigger direct prediction"""
    return os.path.join(TEST_FILES_DIR, "tiny_molecule.png")


@pytest.fixture(scope="module")
def caffeine_image_path():
    """Caffeine image for testing"""
    return os.path.join(TEST_FILES_DIR, "caffeine.png")


# Test the convert_image function
def test_convert_image(sample_gif_path, sample_png_path):
    converted_path = convert_image(sample_gif_path)
    assert os.path.isfile(converted_path)
    assert converted_path == sample_png_path
    # Clean up the converted file
    if os.path.exists(converted_path):
        os.remove(converted_path)


# Test the get_segments function with GIF
def test_get_segments_gif(sample_gif_path):
    image_name, segments = get_segments(sample_gif_path)
    assert image_name == "segment_sample.gif"
    assert isinstance(segments, list)


# Test the get_segments function with PNG
def test_get_segments_png(sample_png_path):
    image_name, segments = get_segments(sample_png_path)
    assert image_name == "segment_sample.png"
    assert isinstance(segments, list)


# Test the get_predicted_segments function
@patch("app.modules.decimer.predict_SMILES")
def test_get_predicted_segments(mock_predict_smiles, sample_png_path):
    mock_predict_smiles.return_value = "CCO"
    predicted_smiles = get_predicted_segments(sample_png_path)
    assert isinstance(predicted_smiles, str)
    assert len(predicted_smiles) > 0


# Test get_predicted_segments_from_file with large image (should use segmentation)
@patch("app.modules.decimer.get_predicted_segments")
def test_get_predicted_segments_from_file_large_image(
    mock_get_predicted_segments, caffeine_image_path
):
    """Test that large images (>=500 pixels) use segmentation approach"""
    mock_get_predicted_segments.return_value = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    with open(caffeine_image_path, "rb") as f:
        content = f.read()

    predicted_smiles = get_predicted_segments_from_file(content, "test_large.png")

    assert isinstance(predicted_smiles, str)
    assert len(predicted_smiles) > 0
    mock_get_predicted_segments.assert_called_once()


# Test get_predicted_segments_from_file with small image (should use direct prediction)
@patch("app.modules.decimer.predict_SMILES")
def test_get_predicted_segments_from_file_small_image(
    mock_predict_smiles, small_image_path
):
    """Test that small images (<500 pixels) use direct prediction"""
    mock_predict_smiles.return_value = "C"

    with open(small_image_path, "rb") as f:
        content = f.read()

    predicted_smiles = get_predicted_segments_from_file(content, "test_small.png")

    assert isinstance(predicted_smiles, str)
    assert len(predicted_smiles) > 0
    mock_predict_smiles.assert_called_once()


# Test get_predicted_segments_from_file with tiny image (should use direct prediction)
@patch("app.modules.decimer.predict_SMILES")
def test_get_predicted_segments_from_file_tiny_image(
    mock_predict_smiles, tiny_image_path
):
    """Test that tiny images (<500 pixels) use direct prediction"""
    mock_predict_smiles.return_value = "CO"

    with open(tiny_image_path, "rb") as f:
        content = f.read()

    predicted_smiles = get_predicted_segments_from_file(content, "test_tiny.png")

    assert isinstance(predicted_smiles, str)
    assert len(predicted_smiles) > 0
    mock_predict_smiles.assert_called_once()


# Test error handling in get_predicted_segments_from_file
def test_get_predicted_segments_from_file_cleanup():
    """Test that temporary files are always cleaned up, even on errors"""
    test_content = b"invalid image content"
    test_filename = "test_cleanup.png"

    # This should fail but still clean up the file
    try:
        get_predicted_segments_from_file(test_content, test_filename)
    except Exception:
        pass  # Expected to fail with invalid image content

    # File should not exist after function completes
    assert not os.path.exists(test_filename)


# Test image size detection logic
def test_image_size_detection():
    """Test that the image size detection works correctly"""
    # Create temporary images with known sizes
    with tempfile.NamedTemporaryFile(
        suffix=".png", delete=False
    ) as tmp_large, tempfile.NamedTemporaryFile(
        suffix=".png", delete=False
    ) as tmp_small:

        try:
            # Create large image (600x600)
            large_img = Image.new("RGB", (600, 600), "white")
            large_img.save(tmp_large.name)

            # Create small image (300x300)
            small_img = Image.new("RGB", (300, 300), "white")
            small_img.save(tmp_small.name)

            # Test with large image content
            with open(tmp_large.name, "rb") as f:
                large_content = f.read()

            # Test with small image content
            with open(tmp_small.name, "rb") as f:
                small_content = f.read()

            # Mock the prediction functions to verify which path is taken
            with patch("app.modules.decimer.predict_SMILES") as mock_direct, patch(
                "app.modules.decimer.get_predicted_segments"
            ) as mock_segment:

                mock_direct.return_value = "direct_prediction"
                mock_segment.return_value = "segmented_prediction"

                # Test large image uses segmentation
                result_large = get_predicted_segments_from_file(
                    large_content, "test_large_600x600.png"
                )
                assert result_large == "segmented_prediction"
                mock_segment.assert_called()
                mock_direct.assert_not_called()

                # Reset mocks
                mock_direct.reset_mock()
                mock_segment.reset_mock()

                # Test small image uses direct prediction
                result_small = get_predicted_segments_from_file(
                    small_content, "test_small_300x300.png"
                )
                assert result_small == "direct_prediction"
                mock_direct.assert_called()
                mock_segment.assert_not_called()

        finally:
            # Clean up temporary files
            for tmp_file in [tmp_large.name, tmp_small.name]:
                if os.path.exists(tmp_file):
                    os.remove(tmp_file)
