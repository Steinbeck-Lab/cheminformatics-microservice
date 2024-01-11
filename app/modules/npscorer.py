from __future__ import annotations

import gzip
import math
import os
import pickle

import pystow
from rdkit.Chem import rdMolDescriptors

# Set path
default_path = pystow.join("NP_model")

# Model download location
model_url = "https://github.com/rdkit/rdkit/blob/master/Contrib/NP_Score/publicnp.model.gz?raw=true"
model_path = str(default_path) + "/publicnp.model.gz"

# Download models to a default location if not already present
if not os.path.exists(model_path):
    pystow.ensure("NP_model", url=model_url)

fscore = pickle.load(gzip.open(model_path))


def get_np_model(model_path) -> dict:
    """
    Load the NP model from a pickle file.

    Parameters:
        model_path (str): Path to the pickled model file.

    Returns:
        dict: The NP model.
    """
    fscore = pickle.load(gzip.open(model_path.as_posix()))
    return fscore


def score_mol_with_confidence(molecule) -> dict:
    """
    Calculate NP-likeness score and confidence for a molecule.

    Args:
        molecule (rdkit.Chem.rdchem.Mol): The input molecule.

    Returns:
        dict: A dictionary containing NP-likeness score and confidence.
            - 'nplikeness' (float): The NP-likeness score.
            - 'confidence' (float): The confidence in the score.
    """
    if molecule is None:
        raise ValueError("Invalid molecule")
    fp = rdMolDescriptors.GetMorganFingerprint(molecule, 2)
    bits = fp.GetNonzeroElements()

    # Calculating the score
    score = 0.0
    bits_found = 0
    for bit in bits:
        if bit in fscore:
            bits_found += 1
            score += fscore[bit]

    score /= float(molecule.GetNumAtoms())
    confidence = float(bits_found / len(bits))

    # Preventing score explosion for exotic molecules
    if score > 4:
        score = 4.0 + math.log10(score - 4.0 + 1.0)
    elif score < -4:
        score = -4.0 - math.log10(-4.0 - score + 1.0)
    result = {"nplikeness": score, "confidence": confidence}
    return result


def score_mol(molecule) -> float:
    """
    Calculate the Natural Product Likeness score for a given molecule.

    Parameters:
        molecule (rdkit.Chem.Mol): RDKit molecule object.

    Returns:
        float: NP-Likeness score in the range -5 to 5.
    """
    score = score_mol_with_confidence(molecule)["nplikeness"]
    return score


def get_np_score(molecule: any) -> str:
    """
    Convert SMILES string to RDKit molecule object and generate the NP Score.

    Parameters:
        molecule (Chem.Mol): RDKit molecule object.

    Returns:
        str: NP Score as a formatted string or "invalid" if conversion fails.
    """
    if molecule:
        npscore = "%.2f" % score_mol(molecule)
    else:
        npscore = "invalid"

    return npscore
