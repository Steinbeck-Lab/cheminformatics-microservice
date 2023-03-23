# calculation of natural product-likeness as described in:
#
# Natural Product-likeness Score and Its Application for
# Prioritization of Compound Libraries
# Peter Ertl, Silvio Roggo, and Ansgar Schuffenhauer
# Journal of Chemical Information and Modeling, 48, 68-74 (2008)
# http://pubs.acs.org/doi/abs/10.1021/ci700286x
# Model was contributed by Peter Ertl :
# https://github.com/rdkit/rdkit/tree/master/Contrib/NP_Score
#

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

# model download location
model_url = "https://github.com/rdkit/rdkit/blob/master/Contrib/NP_Score/publicnp.model.gz?raw=true"
model_path = str(default_path) + "/publicnp.model.gz"

# download models to a default location
if not os.path.exists(model_path):
    pystow.ensure("NP_model", url=model_url)

fscore = pickle.load(gzip.open(model_path))


def getnp_model(model_path):
    fscore = pickle.load(gzip.open(model_path.as_posix()))
    return fscore


def scoremol_with_confidence(mol):
    """Next to the NP Likeness Score, this function outputs a confidence value
    between 0..1 that descibes how many fragments of the tested molecule
    were found in the model data set (1: all fragments were found).

    Returns namedtuple NPLikeness(nplikeness, confidence)"""

    if mol is None:
        raise ValueError("invalid molecule")
    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    bits = fp.GetNonzeroElements()

    # calculating the score
    score = 0.0
    bits_found = 0
    for bit in bits:
        if bit in fscore:
            bits_found += 1
            score += fscore[bit]

    score /= float(mol.GetNumAtoms())
    confidence = float(bits_found / len(bits))

    # preventing score explosion for exotic molecules
    if score > 4:
        score = 4.0 + math.log10(score - 4.0 + 1.0)
    elif score < -4:
        score = -4.0 - math.log10(-4.0 - score + 1.0)
    NPLikeness = namedtuple("NPLikeness", "nplikeness,confidence")
    return NPLikeness(score, confidence)


def score_mol(mol):
    """Calculates the Natural Product Likeness of a molecule.

    Returns the score as float in the range -5..5."""
    return scoremol_with_confidence(mol).nplikeness


def getNPScore(smiles):
    """Converts SMILES to RDKit mol object and generates the NP Score

    Returns the NP Score"""

    mol = Chem.MolFromSmiles(smiles)

    if mol:
        npscore = "%.2f" % score_mol(mol)
    else:
        npscore = "invalid"

    return npscore
