import os
import cv2
from PIL import Image
from decimer_segmentation import segment_chemical_structures_from_file
from DECIMER import predict_SMILES


def convert_image(path: str):
    """Takes an image filepath of GIF image and returns Hi Res PNG image.
    Args:
        input_path (str): path of an image.

    Returns:
        segment_paths (list): a list of paths of segmented images.
    """
    img = Image.open(path).convert("RGBA")
    new_size = int((float(img.width) * 2)), int((float(img.height) * 2))
    resized_image = img.resize(new_size, resample=Image.LANCZOS)
    background_size = int((float(resized_image.width) * 2)), int(
        (float(resized_image.height) * 2)
    )
    new_im = Image.new(resized_image.mode, background_size, "white")
    paste_pos = (
        int((new_im.size[0] - resized_image.size[0]) / 2),
        int((new_im.size[1] - resized_image.size[1]) / 2),
    )
    new_im.paste(resized_image, paste_pos)
    new_im.save(path.replace("gif", "png"), optimize=True, quality=100)
    return path.replace("gif", "png")


def get_segments(path: str):
    """Takes an image filepath and returns a set of paths and image name of segmented images.
    Args:
        input_path (str): path of an image.

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


def getPredictedSegments(path: str):
    """Takes an image filepath and returns predicted SMILES for segmented images.
    Args:
        input_path (str): path of an image.

    Returns:
        predictions (list): a list of SMILES of the segmented images.
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
