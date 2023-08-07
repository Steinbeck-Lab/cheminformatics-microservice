from __future__ import annotations

from typing import Any, List

from pydantic import BaseModel, Field


class ClassyFireJob(BaseModel):
    id: int
    label: str
    finished_at: Any
    created_at: str
    updated_at: str
    query_errors: Any
    finished_processing_at: Any
    query_type: str
    fstruc_file_name: Any
    fstruc_content_type: Any
    fstruc_file_size: Any
    fstruc_updated_at: Any
    query_input: str
    tag_list: List[str]


class Kingdom(BaseModel):
    name: str
    description: str
    chemont_id: str
    url: str


class Superclass(BaseModel):
    name: str
    description: str
    chemont_id: str
    url: str


class Class(BaseModel):
    name: str
    description: str
    chemont_id: str
    url: str


class Subclass(BaseModel):
    name: str
    description: str
    chemont_id: str
    url: str


class DirectParent(BaseModel):
    name: str
    description: str
    chemont_id: str
    url: str


class ExternalDescriptor(BaseModel):
    source: str
    source_id: str
    annotations: List[str]


class Entity(BaseModel):
    identifier: str
    smiles: str
    inchikey: str
    kingdom: Kingdom
    superclass: Superclass
    class_: Class = Field(..., alias="class")
    subclass: Subclass
    intermediate_nodes: List
    direct_parent: DirectParent
    alternative_parents: List
    molecular_framework: str
    substituents: List[str]
    description: str
    external_descriptors: List[ExternalDescriptor]
    ancestors: List[str]
    predicted_chebi_terms: List[str]
    predicted_lipidmaps_terms: List
    classification_version: str


class ClassyFireResult(BaseModel):
    id: int
    label: str
    classification_status: str
    number_of_elements: int
    number_of_pages: int
    invalid_entities: List
    entities: List[Entity]
