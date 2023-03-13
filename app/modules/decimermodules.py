import os
import cv2
from decimer_segmentation import segment_chemical_structures_from_file
from DECIMER import predict_SMILES

def getPredictedSegments(path:str):
    """Takes an image filepath and returns a set of paths of segmented images
    Args:
        input_path (str): path of an image
    
    Returns:
        segment_paths (list): a list of paths of segmented images. 
    """
    smiles_predicted = []
    image_name = os.path.split(path)[1]
    segments = segment_chemical_structures_from_file(path)
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
        return '.'.join(smiles_predicted)
