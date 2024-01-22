from __future__ import annotations

from typing import Any
from typing import List

from pydantic import BaseModel
from pydantic import Field


class ClassyFireJob(BaseModel):
    """Represents a ClassyFire job.

    Attributes:
        id (int): The ID of the job.
        label (str): The label of the job.
        finished_at (Any): The timestamp when the job finished.
        created_at (str): The timestamp when the job was created.
        updated_at (str): The timestamp when the job was last updated.
        query_errors (Any): Any errors related to the query.
        finished_processing_at (Any): The timestamp when processing finished.
        query_type (str): The type of query.
        fstruc_file_name (Any): The name of the file containing structural data.
        fstruc_content_type (Any): The content type of the structural file.
        fstruc_file_size (Any): The size of the structural file.
        fstruc_updated_at (Any): The timestamp when the structural file was updated.
        query_input (str): The input for the query.
        tag_list (List[str]): List of tags associated with the job.
    """

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
    """Represents a Kingdom in ClassyFire taxonomy.

    Attributes:
        name (str): The name of the Kingdom.
        description (str): The description of the Kingdom.
        chemont_id (str): The ChemOnt ID of the Kingdom.
        url (str): The URL associated with the Kingdom.
    """

    name: str
    description: str
    chemont_id: str
    url: str


class Superclass(BaseModel):
    """Represents a Superclass in ClassyFire taxonomy.

    Attributes:
        name (str): The name of the Superclass.
        description (str): The description of the Superclass.
        chemont_id (str): The ChemOnt ID of the Superclass.
        url (str): The URL associated with the Superclass.
    """

    name: str
    description: str
    chemont_id: str
    url: str


class Class(BaseModel):
    """Represents a Class in ClassyFire taxonomy.

    Attributes:
        name (str): The name of the Class.
        description (str): The description of the Class.
        chemont_id (str): The ChemOnt ID of the Class.
        url (str): The URL associated with the Class.
    """

    name: str
    description: str
    chemont_id: str
    url: str


class Subclass(BaseModel):
    """Represents a Subclass within a classification.

    Attributes:
        name (str): The name of the subclass.
        description (str): Description of the subclass.
        chemont_id (str): Identifier for the subclass.
        url (str): URL associated with the subclass.
    """

    name: str
    description: str
    chemont_id: str
    url: str


class DirectParent(BaseModel):
    """Represents a Direct Parent within a classification.

    Attributes:
        name (str): The name of the direct parent.
        description (str): Description of the direct parent.
        chemont_id (str): Identifier for the direct parent.
        url (str): URL associated with the direct parent.
    """

    name: str
    description: str
    chemont_id: str
    url: str


class ExternalDescriptor(BaseModel):
    """Represents an External Descriptor associated with an Entity.

    Attributes:
        source (str): The source of the external descriptor.
        source_id (str): Identifier for the external descriptor.
        annotations (List[str]): List of annotations for the descriptor.
    """

    source: str
    source_id: str
    annotations: List[str]


class Entity(BaseModel):
    """Represents a chemical Entity with classification information.

    Attributes:
        identifier (str): Unique identifier for the entity.
        smiles (str): SMILES representation of the entity.
        inchikey (str): InChIKey of the entity.
        kingdom (Kingdom): Kingdom classification of the entity.
        superclass (Superclass): Superclass classification of the entity.
        class_ (Class): Class classification of the entity.
        subclass (Subclass): Subclass classification of the entity.
        intermediate_nodes (List): List of intermediate nodes.
        direct_parent (DirectParent): Direct parent classification of the entity.
        alternative_parents (List): List of alternative parent classifications.
        molecular_framework (str): The molecular framework of the entity.
        substituents (List[str]): List of substituents in the entity.
        description (str): Description of the entity.
        external_descriptors (List[ExternalDescriptor]): List of external descriptors.
        ancestors (List[str]): List of ancestor classifications.
        predicted_chebi_terms (List[str]): List of predicted ChEBI terms.
        predicted_lipidmaps_terms (List): List of predicted LipidMaps terms.
        classification_version (str): Version of the classification.
    """

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
    """Represents a classification result from ClassyFire.

    Attributes:
        id (int): Identifier for the result.
        label (str): Label for the result.
        classification_status (str): Status of the classification.
        number_of_elements (int): Number of elements in the classification.
        number_of_pages (int): Number of pages.
        invalid_entities (List): List of invalid entities.
        entities (List[Entity]): List of entities in the classification.
    """

    id: int
    label: str
    classification_status: str
    number_of_elements: int
    number_of_pages: int
    invalid_entities: List
    entities: List[Entity]
