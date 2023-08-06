from __future__ import annotations

from pydantic import BaseModel, Field


class MoleculeHash(BaseModel):
    Formula: str
    Isomeric_SMILES: str
    Canonical_SMILES: str


class Representations(BaseModel):
    InChI: str
    InChI_Key: str
    Murko: str


class Descriptors(BaseModel):
    atom_count: int
    heavy_atom_count: int
    molecular_weight: float
    exactmolecular_weight: float
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
    fractioncsp3: int
    number_of_minimal_rings: int
    van_der_walls_volume: str
    linear_sugars: bool
    circular_sugars: bool
    murko_framework: str
    nplikeness: float


class Parent(BaseModel):
    field_2D_mol: str = Field(..., alias="2D_mol")
    field_3D_mol: str = Field(..., alias="3D_mol")
    v3000: str
    representations: Representations
    descriptors: Descriptors


class CococnutPreprocessingModel(BaseModel):
    original_mol: str
    standardised_mol: str
    standardised_SMILES: str
    molecule_hash: MoleculeHash
    parent: Parent
    stereochemical_variants: bool
