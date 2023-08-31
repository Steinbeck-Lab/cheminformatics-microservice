import gzip
import math
import os
import pickle
import pystow
from collections import namedtuple
from rdkit import Chem
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


def getnp_model(model_path) -> dict:
    """
    Load the NP model from a pickle file.

    Parameters:
        model_path (str): Path to the pickled model file.

    Returns:
        dict: The NP model.
    """
    fscore = pickle.load(gzip.open(model_path.as_posix()))
    return fscore


def scoremol_with_confidence(mol) -> dict:
    """
    Calculate NP-likeness score and confidence for a molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The input molecule.

    Returns:
        dict: A dictionary containing NP-likeness score and confidence.
            - 'nplikeness' (float): The NP-likeness score.
            - 'confidence' (float): The confidence in the score.
    """
    if mol is None:
        raise ValueError("Invalid molecule")
    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    bits = fp.GetNonzeroElements()

    # Calculating the score
    score = 0.0
    bits_found = 0
    for bit in bits:
        if bit in fscore:
            bits_found += 1
            score += fscore[bit]

    score /= float(mol.GetNumAtoms())
    confidence = float(bits_found / len(bits))

    # Preventing score explosion for exotic molecules
    if score > 4:
        score = 4.0 + math.log10(score - 4.0 + 1.0)
    elif score < -4:
        score = -4.0 - math.log10(-4.0 - score + 1.0)
    result = {"nplikeness": score, "confidence": confidence}
    return result


def score_mol(mol) -> float:
    """
    Calculate the Natural Product Likeness score for a given molecule.

    Parameters:
        mol (rdkit.Chem.Mol): RDKit molecule object.

    Returns:
        float: NP-Likeness score in the range -5 to 5.
    """
    score = scoremol_with_confidence(mol)["nplikeness"]
    return score


def getNPScore(smiles) -> str:
    """
    Convert SMILES string to RDKit molecule object and generate the NP Score.

    Parameters:
        smiles (str): SMILES representation of a molecule.

    Returns:
        str: NP Score as a formatted string or "invalid" if conversion fails.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        npscore = "%.2f" % score_mol(mol)
    else:
        npscore = "invalid"

    return npscore
