from __future__ import annotations

import os

import cv2
from DECIMER import predict_SMILES
from decimer_segmentation import segment_chemical_structures_from_file
from PIL import Image


def convert_image(path: str) -> str:
    """Convert a GIF image to PNG format, resize, and place on a white.

    background.

    Args:
        path (str): The path to the GIF image file.

    Returns:
        str: The path of the converted and processed PNG image file.
    """
    # Open the image and convert to RGBA
    img = Image.open(path).convert("RGBA")

    # Calculate new dimensions for resizing
    new_size = (int(float(img.width) * 2), int(float(img.height) * 2))

    # Resize the image using Lanczos resampling
    resized_image = img.resize(new_size, resample=Image.LANCZOS)

    # Calculate background size
    background_size = (
        int(float(resized_image.width) * 2),
        int(float(resized_image.height) * 2),
    )

    # Create a new image with white background
    new_im = Image.new(resized_image.mode, background_size, "white")

    # Calculate position to paste the resized image in the center
    paste_pos = (
        int((new_im.size[0] - resized_image.size[0]) / 2),
        int((new_im.size[1] - resized_image.size[1]) / 2),
    )

    # Paste the resized image onto the new background
    new_im.paste(resized_image, paste_pos)

    # Save the processed image in PNG format
    new_path = path.replace("gif", "png")
    new_im.save(new_path, optimize=True, quality=100)

    return new_path


def get_segments(path: str) -> tuple:
    """Takes an image file path and returns a set of paths and image names of.

    segmented images.

    Args:
        input_path (str): the path of an image.

    Returns:
        image_name (str): image file name.
        segments (list): a set of segmented images.
    """

    image_name = os.path.split(path)[1]
    if image_name[-3:].lower() == "gif":
        new_path = convert_image(path)
        segments = segment_chemical_structures_from_file(new_path)
        return image_name, segments
    else:
        segments = segment_chemical_structures_from_file(path)
        return image_name, segments


def get_predicted_segments(path: str) -> str:
    """Get predicted SMILES representations for segments within an image.

    This function takes an image path, extracts segments, predicts SMILES representations
    for each segment, and returns a concatenated string of predicted SMILES.

    Args:
        path (str): Path to the input image file.

    Returns:
        str: Predicted SMILES representations joined by '.' if segments are detected,
             otherwise returns a single predicted SMILES for the whole image.
    """
    smiles_predicted = []
    image_name, segments = get_segments(path)

    if len(segments) == 0:
        smiles = predict_SMILES(path)
        return smiles
    else:
        for segment_index in range(len(segments)):
            segmentname = f"{image_name[:-5]}_{segment_index}.png"
            segment_path = os.path.join(segmentname)
            cv2.imwrite(segment_path, segments[segment_index])
            smiles = predict_SMILES(segment_path)
            smiles_predicted.append(smiles)
            os.remove(segment_path)
        return ".".join(smiles_predicted)


def get_predicted_segments_from_file(content: any, filename: str) -> tuple:
    """Takes an image file path and returns a set of paths and image names of.

    segmented images.

    Args:
        input_path (str): the path of an image.

    Returns:
        image_name (str): image file name.
        segments (list): a set of segmented images.
    """

    with open(filename, "wb") as f:
        f.write(content)
        smiles = get_predicted_segments(filename)
        os.remove(filename)
        return smiles
