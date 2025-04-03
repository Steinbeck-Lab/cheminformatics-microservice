from __future__ import annotations

from pydantic import BaseModel
from pydantic import Field


class MoleculeHash(BaseModel):
    """Represents a molecule's hash information.

    Attributes:
        Formula (str): Chemical formula of the molecule.
        Isomeric_SMILES (str): Isomeric Simplified Molecular Input Line Entry System (SMILES) representation.
        Canonical_SMILES (str): Canonical SMILES representation.
    """

    Formula: str
    Isomeric_SMILES: str
    Canonical_SMILES: str


class Representations(BaseModel):
    """Represents different representations of a molecule.

    Attributes:
        field_2D_mol (str): The 2D molecular structure of the parent molecule.
        field_3D_mol (str): The 3D molecular structure of the parent molecule.
        InChI (str): IUPAC International Chemical Identifier (InChI) representation.
        InChI_Key (str): InChI key.
        Murcko (str): Murcko scaffold representation.
    """

    field_2D_mol: str = Field(..., alias="2D_mol")
    field_3D_mol: str = Field(..., alias="3D_mol")
    InChI: str
    InChI_Key: str
    Murcko: str


class Descriptors(BaseModel):
    """Represents a collection of molecular descriptors.

    This class provides a structure for storing various molecular descriptors
    that characterize a molecule's chemical properties.

    Attributes:
        atom_count (int): Total count of atoms in the molecule.
        heavy_atom_count (int): Count of heavy (non-hydrogen) atoms.
        molecular_weight (float): Molecular weight of the molecule.
        exact_molecular_weight (float): Exact molecular weight.
        alogp (float): Calculated ALogP (partition coefficient) value.
        rotatable_bond_count (int): Count of rotatable bonds in the molecule.
        topological_polar_surface_area (int): Topological polar surface area.
        hydrogen_bond_acceptors (int): Count of hydrogen bond acceptors.
        hydrogen_bond_donors (int): Count of hydrogen bond donors.
        hydrogen_bond_acceptors_lipinski (int): Lipinski's count of hydrogen bond acceptors.
        hydrogen_bond_donors_lipinski (int): Lipinski's count of hydrogen bond donors.
        lipinski_rule_of_five_violations (int): Count of Lipinski's rule-of-five violations.
        aromatic_rings_count (int): Count of aromatic rings in the molecule.
        qed_drug_likeliness (float): Quantitative Estimate of Drug-likeness (QED) value.
        formal_charge (int): Formal charge of the molecule.
        fraction_csp3 (int): Fraction of carbon atoms that are sp3 hybridized.
        number_of_minimal_rings (int): Count of minimal rings in the molecule.
        van_der_waals_volume (str): Van der Waals volume description.
        linear_sugars (bool): True if linear sugars are present, False otherwise.
        circular_sugars (bool): True if circular sugars are present, False otherwise.
        Murcko_framework (str): Murcko framework description.
        nplikeness (float): NP-likeness value.
    """

    atom_count: int
    heavy_atom_count: int
    molecular_weight: float
    exact_molecular_weight: float
    alogp: float
    rotatable_bond_count: int
    topological_polar_surface_area: int
    hydrogen_bond_acceptors: int
    hydrogen_bond_donors: int
    hydrogen_bond_acceptors_lipinski: int
    hydrogen_bond_donors_lipinski: int
    lipinski_rule_of_five_violations: int
    aromatic_rings_count: int
    qed_drug_likeliness: float
    formal_charge: int
    fraction_csp3: int
    number_of_minimal_rings: int
    van_der_waals_volume: str
    linear_sugars: bool
    circular_sugars: bool
    Murcko_framework: str
    nplikeness: float


class Original(BaseModel):
    """Represents the parent molecule with various properties.

    Attributes:
        representations (Representations): Molecular representations.
        has_stereo (bool): Indicates presence of stereochemical variants.
        descriptors (Descriptors): Molecular descriptors.
        errors (dict): Information on the errors found in the given molecule
    """

    representations: Representations
    has_stereo: bool
    descriptors: Descriptors
    errors: dict


class Standardized(BaseModel):
    """Represents the parent molecule with various properties.

    Attributes:
        representations (Representations): Molecular representations.
        has_stereo (bool): Indicates presence of stereochemical variants.
        descriptors (Descriptors): Molecular descriptors.
        errors (dict): Information on the errors found in the given molecule
    """

    representations: Representations
    has_stereo: bool
    descriptors: Descriptors
    errors: dict


class Parent(BaseModel):
    """Represents the parent molecule with various properties.

    Attributes:
        representations (Representations): Molecular representations.
        has_stereo (bool): Indicates presence of stereochemical variants.
        descriptors (Descriptors): Molecular descriptors.
        errors (dict): Information on the errors found in the given molecule
    """

    representations: Representations
    has_stereo: bool
    descriptors: Descriptors
    errors: dict


class COCONUTPreprocessingModel(BaseModel):
    """Represents a molecule after CocoNut preprocessing.

    Attributes:
        original_mol (Original): Original molecule information.
        standardised_mol (Standardized): Standardized molecule information.
        parent (Parent): Parent molecule information.
    """

    original_mol: Original
    standardised_mol: Standardized
    parent: Parent
