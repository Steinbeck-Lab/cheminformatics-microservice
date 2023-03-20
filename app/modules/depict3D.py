import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

def get3Dconformers(smiles):
    '''Convert SMILES to Mol with 3D coordinates
    Args (str): SMILES string.
    Returns (rdkil.mol): A mol object with 3D coodinates
    optimized with MMFF94 forcefield.

    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        AllChem.Compute2DCoords(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        return mol
    else:
        return None

def get3DDepiction(smiles, size=(512, 512), style="stick", surface=False, opacity=0.5):
    """Draw molecule in 3D
    
    Args:
    ----
        smiles (str): SMILES string.
        size (tuple): canvas size.
        style (str): Type of drawing molecule
        	('line', 'stick', 'sphere', 'carton').
        surface (bool): display SAS.
        opacity (float): opacity of surface, range 0.0-1.0.
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js view.
    """
    assert style in ('line', 'stick', 'sphere', 'carton')
    conformers = get3Dconformers(smiles)
    mblock = Chem.MolToMolBlock(conformers)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, 'mol')
    viewer.setStyle({style:{}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    return viewer
