from __future__ import annotations

import io
import logging
from typing import Annotated
from typing import Literal
from typing import Optional
from typing import Union

import pandas as pd
from chembl_structure_pipeline import checker
from chembl_structure_pipeline import standardizer
from fastapi import APIRouter
from fastapi import Body
from fastapi import HTTPException
from fastapi import Query
from fastapi import status
from fastapi.responses import JSONResponse
from fastapi.responses import Response
from fastapi.templating import Jinja2Templates
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
)

from app.modules.all_descriptors import get_tanimoto_similarity
from app.modules.classyfire import classify
from app.modules.classyfire import result
from app.modules.coconut.descriptors import get_COCONUT_descriptors
from app.modules.coconut.preprocess import get_COCONUT_preprocessing
from app.modules.npscorer import get_np_score
from app.modules.toolkits.cdk_wrapper import get_CDK_HOSE_codes
from app.modules.toolkits.cdk_wrapper import get_tanimoto_similarity_CDK
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.rdkit_wrapper import check_RO5_violations
from app.modules.toolkits.rdkit_wrapper import check_RO5_violations_detailed
from app.modules.toolkits.rdkit_wrapper import get_ertl_functional_groups
from app.modules.toolkits.rdkit_wrapper import get_GhoseFilter
from app.modules.toolkits.rdkit_wrapper import get_GhoseFilter_detailed
from app.modules.toolkits.rdkit_wrapper import get_PAINS
from app.modules.toolkits.rdkit_wrapper import get_PAINS_detailed
from app.modules.toolkits.rdkit_wrapper import get_properties
from app.modules.toolkits.rdkit_wrapper import get_rdkit_HOSE_codes
from app.modules.toolkits.rdkit_wrapper import get_REOSFilter
from app.modules.toolkits.rdkit_wrapper import get_REOSFilter_detailed
from app.modules.toolkits.rdkit_wrapper import get_RuleofThree
from app.modules.toolkits.rdkit_wrapper import get_RuleofThree_detailed
from app.modules.toolkits.rdkit_wrapper import get_sas_score
from app.modules.toolkits.rdkit_wrapper import get_tanimoto_similarity_rdkit
from app.modules.toolkits.rdkit_wrapper import get_VeberFilter
from app.modules.toolkits.rdkit_wrapper import get_VeberFilter_detailed
from app.modules.toolkits.rdkit_wrapper import get_standardized_tautomer
from app.modules.toolkits.rdkit_wrapper import QED
from app.modules.pubchem_retrieve import PubChemClient
from app.schemas import HealthCheck
from app.schemas.chem_schema import FilteredMoleculesResponse
from app.schemas.chem_schema import GenerateDescriptorsResponse
from app.schemas.chem_schema import GenerateFunctionalGroupResponse
from app.schemas.chem_schema import GenerateHOSECodeResponse
from app.schemas.chem_schema import GenerateMultipleDescriptorsResponse
from app.schemas.chem_schema import GenerateStandardizeResponse
from app.schemas.chem_schema import GenerateStereoisomersResponse
from app.schemas.chem_schema import NPlikelinessScoreResponse
from app.schemas.chem_schema import TanimotoMatrixResponse
from app.schemas.chem_schema import TanimotoSimilarityResponse
from app.schemas.chem_schema import StandarizedTautomerResponse
from app.schemas.chemblstandardizer import SMILESStandardizedResult
from app.schemas.chemblstandardizer import SMILESValidationResult
from app.schemas.classyfire import ClassyFireJob
from app.schemas.classyfire import ClassyFireResult
from app.schemas.pubchem_schema import PubChemResponse
from app.schemas.coconut import COCONUTPreprocessingModel
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse
from app.schemas.error import NotFoundModel


router = APIRouter(
    prefix="/chem",
    tags=["chem"],
    dependencies=[],
    responses={
        200: {"description": "OK"},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)

templates = Jinja2Templates(directory="app/templates")

pubchem_client = PubChemClient()
logger = logging.getLogger(__name__)


@router.get("/", include_in_schema=False)
@router.get(
    "/health",
    tags=["healthcheck"],
    summary="Perform a Health Check on Chem Module",
    response_description="Return HTTP Status Code 200 (OK)",
    status_code=status.HTTP_200_OK,
    include_in_schema=False,
    response_model=HealthCheck,
)
def get_health() -> HealthCheck:
    """## Perform a Health Check.

    Endpoint to perform a health check on. This endpoint can primarily be used by Docker
    to ensure a robust container orchestration and management are in place. Other
    services that rely on the proper functioning of the API service will not deploy if this
    endpoint returns any other HTTP status code except 200 (OK).
    Returns:
        HealthCheck: Returns a JSON response with the health status
    """
    return HealthCheck(status="OK")


@router.get(
    "/stereoisomers",
    summary="Enumerate all possible stereoisomers",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateStereoisomersResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def get_stereoisomers(
    smiles: str = Query(
        title="SMILES",
        description="SMILES string to be enumerated",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
):
    """For a given SMILES string this function enumerates all possible.

    stereoisomers.

    Parameters:
    - **SMILES**: required (query parameter): The SMILES string to be enumerated.

    Returns:
    - List[str]: A list of stereoisomer SMILES strings if successful, otherwise returns an error message.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    mol = parse_input(smiles, "rdkit", False)
    if mol:
        isomers = tuple(EnumerateStereoisomers(mol))
        smilesArray = []
        for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
            smilesArray.append(smi)
        return smilesArray


@router.get(
    "/descriptors",
    summary="Generates descriptors for the input molecule",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateDescriptorsResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def get_descriptors(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    format: Literal["json", "html"] = Query(
        default="json",
        description="Desired display format",
    ),
    toolkit: Literal["cdk", "rdkit", "all"] = Query(
        default="rdkit",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Generates standard descriptors for the input molecule (SMILES).

    Parameters:
    - **SMILES**: required (query): The SMILES representation of the molecule.
    - **format**: optional (query): The desired format for the output.
        - Supported values: "html" / "json" (default), "json".
    - **toolkit**: optional (query): The chemical toolkit to use for descriptor calculation. The default is "rdkit". Allowed "all", "cdk".
        - Supported values: "cdk"/ "rdkit" / "all" (default), "rdkit".

    Returns:
    - If the format is "html", return an HTML response containing a table of the descriptors and their values.
    - Return the descriptors and their values in the specified format if the format is not "html" (default: JSON).

    Raises:
    - None
    """
    data = get_COCONUT_descriptors(smiles, toolkit)
    if format == "html":
        if toolkit == "all":
            headers = [
                "Descriptor name",
                "RDKit Descriptors",
                "CDK Descriptors",
            ]
            df = pd.DataFrame.from_dict(
                data,
                orient="index",
                columns=headers[1:],
                dtype=object,
            )
            df.insert(0, headers[0], df.index)
        else:
            headers = ["Descriptor name", "Values"]
            df = pd.DataFrame.from_dict(
                data,
                orient="index",
                columns=headers[1:],
            )
            df.insert(0, headers[0], df.index)

        with open("app/templates/style.css", "r") as file:
            css_style = file.read()

        html_table = df.to_html(index=False)
        return Response(content=css_style + html_table, media_type="text/html")
    else:
        return JSONResponse(content=data)


@router.get(
    "/descriptors/multiple",
    summary="Generates descriptors for the input molecules",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateMultipleDescriptorsResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def get_multiple_descriptors(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecules",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine, Topiramate-13C6",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C, CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6, Ethane",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1, CC",
            },
        },
    ),
    toolkit: Literal["cdk", "rdkit"] = Query(
        default="rdkit",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Retrieve multiple descriptors for a list of SMILES strings.

    Parameters:
    - **SMILES**: required (query): Comma-separated list of SMILES strings.
    - **toolkit**: optional (query): Toolkit to use for descriptor calculation.
        - Supported values: "rdkit" / "cdk" (default), "rdkit".

    Returns:
    - Union[Dict[str, Any], str]: If multiple SMILES are provided, return a dictionary with each SMILES as the key and the corresponding descriptors as the value. If only one SMILES is provided, returns an error message.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    Example:
    - Request: GET /descriptors/multiple?smiles=CCO,CCN&toolkit=rdkit
      Response: {"CCO": {"descriptor1": value1, "descriptor2": value2}, "CCN": {"descriptor1": value3, "descriptor2": value4}}

    - Request: GET /descriptors/multiple?smiles=CCC
      Response: "Error invalid SMILES"
    """
    molecules = [m.strip() for m in smiles.split(",")]

    if len(molecules) < 2:
        raise HTTPException(
            status_code=422,
            detail="At least two molecules are required.",
        )

    descriptors_dict = {}

    for molecule in molecules:
        descriptors = get_COCONUT_descriptors(molecule, toolkit)
        # print(type(descriptors))
        descriptors_dict[molecule] = descriptors

    return JSONResponse(content=descriptors_dict)


@router.get(
    "/HOSEcode",
    summary="Generates HOSE codes for the input molecules",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateHOSECodeResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def hose_codes(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    spheres: int = Query(
        title="spheres",
        description="Number of spheres to use for generating HOSE codes",
        examples=[
            "2",
            "1",
        ],
    ),
    toolkit: Literal["cdk", "rdkit"] = Query(
        default="rdkit",
        description="Cheminformatics toolkit used in the backend",
    ),
    ringsize: bool = Query(
        False,
        title="ringsize",
        description="Determines whether to include information about ring sizes",
    ),
):
    """Generates HOSE codes for a given SMILES string.

    Parameters:
    - **SMILES**: required (query): The SMILES string represents the chemical compound.
    - **spheres**: required (query): The number of spheres to use for generating HOSE codes.
    - **toolkit**: optional (default: cdk): The chemical toolkit to use for generating HOSE codes.
            Supported values: "cdk" (default), "rdkit".
    - **ringsize**: optional (default: False): Determines whether to include information about ring sizes
            in the HOSE codes. The default is False.

    Returns:
    - List[str]: A list of HOSE codes if successful, indicating the HOSE codes
        for each atom in the molecule. Otherwise, returns an error message.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    if toolkit == "cdk":
        mol = parse_input(smiles, "cdk", False)
        hose_codes = await get_CDK_HOSE_codes(mol, spheres, ringsize)
    elif toolkit == "rdkit":
        mol = parse_input(smiles, "rdkit", False)
        hose_codes = await get_rdkit_HOSE_codes(mol, spheres)

    if hose_codes:
        return hose_codes
    else:
        raise HTTPException(
            status_code=500,
            detail="Failed to generate HOSE codes, check settings.",
        )


@router.post(
    "/standardize",
    summary="Standardize molblock using the ChEMBL curation pipeline",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateStandardizeResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def standardize_mol(
    data: Annotated[
        str,
        Body(
            embed=False,
            media_type="text/plain",
            openapi_examples={
                "example1": {
                    "summary": "Example: C",
                    "value": """
  CDK     09012310592D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END""",
                },
            },
        ),
    ],
):
    """Standardize molblock using the ChEMBL curation pipeline.

    and return the standardized molecule, SMILES, InChI, and InCHI-Key.

    Parameters:
    - **molblock**: The request body containing the "molblock" string representing the molecule to be standardized.

    Returns:
    - dict: A dictionary containing the following keys:
        - "standardized_mol" (str): The standardized molblock of the molecule.
        - "canonical_smiles" (str): The canonical SMILES representation of the molecule.
        - "inchi" (str): The InChI representation of the molecule.
        - "inchikey" (str): The InChI-Key of the molecule.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    try:
        if data:
            suppl = Chem.SDMolSupplier()
            suppl.SetData(data.encode("utf-8"))
            if len(suppl) == 1 and suppl[0]:
                mol_data = suppl[0]
            mol = Chem.MolToMolBlock(mol_data)
            standardized_mol = standardizer.standardize_molblock(mol)
            rdkit_mol = Chem.MolFromMolBlock(standardized_mol)

        else:
            raise HTTPException(
                status_code=422,
                detail="Invalid or missing molblock",
            )

        if rdkit_mol:
            smiles = Chem.MolToSmiles(rdkit_mol, kekuleSmiles=True)
            response = dict(
                standardized_mol=standardized_mol,
                canonical_smiles=smiles,
                inchi=Chem.inchi.MolToInchi(rdkit_mol),
                inchikey=Chem.inchi.MolToInchiKey(rdkit_mol),
            )
            original_properties = get_properties(data)
            response.update(original_properties)
            return response

    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/errors",
    summary="Check a given SMILES string and the represented structure for issues and standardize it",
    responses={
        200: {
            "description": "Successful response",
            "model": Union[SMILESStandardizedResult, SMILESValidationResult],
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def check_errors(
    smiles: str = Query(
        title="SMILES",
        description="The SMILES string to check and standardize.",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    fix: bool = Query(
        False,
        title="Fix",
        description="Flag indicating whether to fix the issues by standardizing the SMILES.",
    ),
):
    """Check a given SMILES string and the represented structure for issues and.

    standardize it using the ChEMBL curation pipeline.

    Parameters:
    - **SMILES**: required (str, query) The SMILES string to check and standardize.
    - **fix**: optional (bool): Flag indicating whether to fix the issues by standardizing the SMILES. Defaults to False.

    Returns:
    - If the fix is False:
        - If issues are found in the SMILES string, return a list of issues.
        - If no issues are found, return the string "No Errors Found".

    - If the fix is True:
        - If issues are found in the SMILES string, return a dictionary containing the original SMILES, original issues,
          standardized SMILES, and new issues after standardization.
        - If no issues are found after standardization, return a dictionary with the original SMILES and "No Errors Found".

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    Notes:
    - If the SMILES string contains spaces, they will be replaced with "+" characters before processing.
    - If the SMILES string cannot be read, the function returns the string "Error reading SMILES string, check again."
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol:
        mol_block = Chem.MolToMolBlock(mol)
        issues = checker.check_molblock(mol_block)

        if issues:
            if fix:
                standardized_mol = standardizer.standardize_molblock(mol_block)
                issues_new = checker.check_molblock(standardized_mol)
                rdkit_mol = Chem.MolFromMolBlock(standardized_mol)
                standardized_smiles = Chem.MolToSmiles(rdkit_mol)
                result = SMILESStandardizedResult(
                    original=SMILESValidationResult(
                        smi=smiles,
                        messages=issues,
                    ),
                    standardized=SMILESValidationResult(
                        smi=standardized_smiles,
                        messages=issues_new if issues_new else ("No Errors Found",),
                    ),
                )
                return result
            else:
                return SMILESValidationResult(smi=smiles, messages=issues)
        else:
            return SMILESValidationResult(smi=smiles, messages=("No Errors Found",))
    else:
        raise HTTPException(
            status_code=422,
            detail="Error reading SMILES string, please check again.",
        )


@router.get(
    "/nplikeness/score",
    summary="Generates descriptors for the input molecules",
    responses={
        200: {
            "description": "Successful response",
            "model": NPlikelinessScoreResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def np_likeness_score(
    smiles: str = Query(
        title="SMILES",
        description="The SMILES string to calculate the natural product likeness score",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
):
    """Calculates the natural product likeness score based on the RDKit.

    implementation.

    Parameters:
    - **SMILES**: required (query): The SMILES representation of the molecule.

    Returns:
    - np_score (float): The natural product likeness score.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    mol = parse_input(smiles, "rdkit", False)
    try:
        np_score = get_np_score(mol)
        if np_score:
            return float(np_score)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/tanimoto",
    summary="Generates the Tanimoto similarity index for a given pair of SMILES strings",
    responses={
        200: {
            "description": "Successful response",
            "model": Union[TanimotoSimilarityResponse, TanimotoMatrixResponse],
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def tanimoto_similarity(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecules",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine, Topiramate-13C6",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C,CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6, Ethane",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1,CC",
            },
        },
    ),
    toolkit: Literal["cdk", "rdkit"] = Query(
        default="rdkit",
        description="Cheminformatics toolkit used in the backend",
    ),
    fingerprinter: Literal[
        "RDKit", "Atompairs", "MACCS", "PubChem", "ECFP", "MAPC"
    ] = Query(
        "ECFP",
        description="Molecule fingerprint generation algorithm to use for RDKit and CDK",
    ),
    nBits: Optional[int] = Query(
        "2048",
        title="nBits size",
        description="The number of bits for fingerprint vectors in RDKit and CDK. Ignored for MACCS and PubChem keys.",
    ),
    radius: Optional[int] = Query(
        "2",
        title="radius size - ECFP",
        description="ECFP 2/4/6 are allowed for using CDK Circular fingerprinter and for MAPC default is radius 2 and permutations 2048. The default is 6",
    ),
):
    """Calculate the Tanimoto similarity index for a pair of SMILES strings.

    using specified parameters.

    Args:
        smiles (str): A comma-separated pair of SMILES strings for the molecules to compare.
        toolkit (Literal["cdk", "rdkit"]): The cheminformatics toolkit used in the backend.
        fingerprinter (Literal["RDKit", "Atompairs", "MACCS","Pubchem","ECFP"]): The molecule fingerprint generation algorithm to use for RDKit (supports: "RDKit", "Atompairs", "MACCS" and "ECFP") and CDK (supports: "Pubchem" and"ECFP").
        nBits (int, optional): The number of bits for fingerprint vectors in RDKit. Ignored for MACCS keys. Defaults to 2048.
        radius (int, optional): The radius size for ECFP (Circular Fingerprints) when using CDK. Defaults to 6.

    Returns:
        The Tanimoto similarity index as a floating-point value.

    Raises:
        HTTPException: Raised when there is an error reading SMILES strings or invalid input.
    """

    molecules = [m.strip() for m in smiles.split(",")]

    if len(molecules) < 2:
        raise HTTPException(
            status_code=422,
            detail="At least two molecules are required.",
        )

    elif len(smiles.split(",")) == 2:
        try:
            smiles1, smiles2 = smiles.split(",")
            if toolkit == "rdkit":
                mol1 = parse_input(smiles1, "rdkit", False)
                mol2 = parse_input(smiles2, "rdkit", False)
                Tanimoto = get_tanimoto_similarity_rdkit(
                    mol1,
                    mol2,
                    fingerprinter,
                    radius,
                    nBits,
                )
            else:
                mol1 = parse_input(smiles1, "cdk", False)
                mol2 = parse_input(smiles2, "cdk", False)
                Tanimoto = get_tanimoto_similarity_CDK(
                    mol1,
                    mol2,
                    fingerprinter,
                    radius,
                    nBits,
                )
            return float(Tanimoto)
        except Exception:
            raise HTTPException(
                status_code=422,
                detail="Error reading SMILES string, please check again.",
            )

    elif len(smiles.split(",")) > 2:
        try:
            matrix = get_tanimoto_similarity(smiles, toolkit)
            return Response(content=matrix, media_type="text/html")
        except Exception:
            raise HTTPException(
                status_code=422,
                detail="Error reading SMILES string, please check again.",
            )
    else:
        raise HTTPException(
            status_code=422,
            detail="Error reading SMILES string, please check again.",
        )


@router.get(
    "/coconut/pre-processing",
    summary="Generates an Input JSON file with information for COCONUT database",
    responses={
        200: {
            "description": "Successful response",
            "model": COCONUTPreprocessingModel,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def coconut_preprocessing(
    smiles: str = Query(
        ...,
        title="SMILES",
        description="SMILES string representing a chemical compound",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    _3d_mol: bool = Query(
        False,
        title="3D_mol",
        description="Flag indicating whether to generate 3D coordinates for a given molecule",
    ),
    descriptors: bool = Query(
        False,
        title="descriptors",
        description="Flag indicating whether to generate COCONUT descriptors for a given molecule",
    ),
):
    """Generates an Input JSON file with information on different molecular.

    representations and descriptors suitable for submission to the COCONUT
    database.

    Parameters:
    - **SMILES**: required (query): The SMILES string represents a chemical compound.

    Returns:
    - JSONResponse: The generated Input JSON file for COCONUT.

    Raises:
    - HTTPException: If there is an error reading the SMILES string.
    """
    try:
        data = get_COCONUT_preprocessing(smiles, _3d_mol, descriptors)
        if data:
            return data
        else:
            raise HTTPException(
                status_code=422,
                detail="Error reading SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(
            status_code=422,
            detail="Error processing request: " + str(e),
        )


@router.get(
    "/classyfire/classify",
    summary="Generate ClassyFire-based classifications using SMILES as input",
    responses={
        200: {
            "description": "Successful response",
            "model": ClassyFireJob,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def classyfire_classify(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the compound to be classified",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
):
    """Generate ClassyFire-based classifications using SMILES as input.

    Parameters:
    - **SMILES**: required (query): The SMILES representation of the compound to be classified.

    Returns:
    - The classification data generated by ClassyFire.

    Raises:
    - HTTPException: If the SMILES string is not provided or if an error occurs during the classification process.

    Note:
    - ClassyFire is a chemical taxonomy classification tool that predicts the chemical class and subclass of a compound based on its structural features.
    - This service pings the http://classyfire.wishartlab.com server for information retrieval.
    """
    mol = parse_input(smiles, "rdkit", False)
    if mol:
        classification_data = await classify(smiles)
        if classification_data:
            return classification_data


@router.get(
    "/classyfire/{jobid}/result",
    summary="Retrieve the ClassyFire classification results based on the provided Job ID",
    responses={
        200: {
            "description": "Successful response",
            "model": ClassyFireResult,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def classyfire_result(jobid: str):
    """Retrieve the ClassyFire classification results based on the provided Job.

    ID.

    To obtain the results from ClassyFire, please initiate a new request and obtain a unique job ID.
    Once you have the job ID, you need to submit another request using this ID in order to retrieve the desired outcome.

    Parameters:
    - **jobid**: required (query): The Job ID used to query the ClassyFire classification results.

    Raises:
    - HTTPException 422: If the Job ID is not provided.
    - HTTPException 500: If an error occurs during the classification process.

    Returns:
    - The ClassyFire classification results as JSON.
    """

    if jobid:
        try:
            # Replace with your function to retrieve the result
            data = await result(jobid)
            return data
        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail="Error processing request: " + str(e),
            )
    else:
        raise HTTPException(status_code=422, detail="Job ID is required.")


@router.post(
    "/all_filters",
    summary="Filters a given list of molecules using selected filters",
    responses={
        200: {
            "description": "Successful response",
            "model": FilteredMoleculesResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def all_filter_molecules(
    smiles_list: str = Body(
        embed=False,
        media_type="text/plain",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    pains: bool = Query(
        True,
        title="PAINS filter",
        description="Calculate PAINS filter (true or false)",
    ),
    lipinski: bool = Query(
        True,
        title="Lipinski Rule of 5",
        description="Calculate Lipinski Rule of 5 (true or false)",
    ),
    veber: bool = Query(
        True,
        title="Veber filter",
        description="Calculate Veber filter (true or false)",
    ),
    reos: bool = Query(
        True,
        title="REOS filter",
        description="Calculate REOS filter (true or false)",
    ),
    ghose: bool = Query(
        True,
        title="Ghose filter",
        description="Calculate Ghose filter (true or false)",
    ),
    ruleofthree: bool = Query(
        True,
        title="Rule of 3",
        description="Calculate Rule of 3 filter (true or false)",
    ),
    qedscore: str = Query(
        "0-10",
        title="QED Druglikeliness",
        description="Calculate QED Score in the range (e.g., 0-10)",
    ),
    sascore: str = Query(
        "0-10",
        title="SAScore",
        description="Calculate SAScore in the range (e.g., 0-10)",
    ),
    nplikeness: str = Query(
        "0-10",
        title="NPlikenessScore",
        description="Calculate NPlikenessScore in the range (e.g., 0-10)",
    ),
    filterOperator: Literal["AND", "OR"] = Query(
        "OR",
        title="Filter Operator",
        description="Logic for combining filter results: AND requires all selected filters to pass, OR requires at least one filter to pass",
    ),
):
    all_smiles = []
    for item in io.StringIO(smiles_list):
        smiles = item.strip()
        mol = parse_input(smiles, "rdkit", False)
        results = []
        results.append(smiles + ":")
        if mol:
            if pains:
                pains_status = str(get_PAINS(mol))
                # PAINS logic: Finding a PAINS is BAD (False), not finding is GOOD (True)
                if "family" in pains_status:
                    results.append("False")  # Found PAINS = BAD = False
                else:
                    results.append("True")  # No PAINS = GOOD = True
            if lipinski:
                lipinski_violations = check_RO5_violations(mol)
                if lipinski_violations == 0:
                    results.append("True")
                else:
                    results.append("False")

            if veber:
                if get_VeberFilter(mol) == "True":
                    results.append("True")
                else:
                    results.append("False")
            if reos:
                if get_REOSFilter(mol) == "True":
                    results.append("True")
                else:
                    results.append("False")

            if ghose:
                if get_GhoseFilter(mol) == "True":
                    results.append("True")
                else:
                    results.append("False")

            if ruleofthree:
                if get_RuleofThree(mol) == "True":
                    results.append("True")
                else:
                    results.append("False")

            if len(qedscore.split("-")) == 2:
                range_start, range_end = float(qedscore.split("-")[0]), float(
                    qedscore.split("-")[1],
                )
                qed_value = QED.qed(mol)
                if range_start <= qed_value <= range_end:
                    results.append("True")
                else:
                    results.append("False")

            if len(sascore.split("-")) == 2:
                range_start, range_end = float(sascore.split("-")[0]), float(
                    sascore.split("-")[1],
                )
                sascore_value = get_sas_score(mol)
                if range_start <= sascore_value <= range_end:
                    results.append("True")
                else:
                    results.append("False")

            if len(nplikeness.split("-")) == 2:
                range_start, range_end = float(nplikeness.split("-")[0]), float(
                    nplikeness.split("-")[1],
                )
                nplikeness_value = get_np_score(mol)
                if range_start <= float(nplikeness_value) <= range_end:
                    results.append("True")
                else:
                    results.append("False")

            # Filter results based on AND/OR operator
            # First element is the SMILES, so we start checking from index 1
            filter_results = (
                [r == "True" for r in results[1:]] if len(results) > 1 else []
            )

            # Apply filter operator logic
            passes_filters = False
            if filter_results:
                if filterOperator == "AND":
                    passes_filters = all(filter_results)  # All filters must pass
                else:  # OR
                    passes_filters = any(
                        filter_results
                    )  # At least one filter must pass

            # Only include SMILES that pass the filter operator logic
            if not filter_results or passes_filters:
                output_list = [
                    x if x != "True" and x != "False" else ("T" if x == "True" else "F")
                    for x in results
                ]
                final_results = ", ".join(output_list).replace(":,", " : ")
                all_smiles.append(final_results)

    return all_smiles


@router.post(
    "/all_filters_detailed",
    summary="Filters molecules with detailed information about violations and substructure matches",
    responses={
        200: {
            "description": "Successful response with detailed filter information",
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def all_filter_molecules_detailed(
    smiles_list: str = Body(
        embed=False,
        media_type="text/plain",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine and problematic compound",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
    pains: bool = Query(
        True,
        title="PAINS filter",
        description="Check for Pan-Assay Interference compounds (true or false)",
    ),
    lipinski: bool = Query(
        True,
        title="Lipinski Rule of 5",
        description="Check Lipinski Rule of 5 violations (true or false)",
    ),
    veber: bool = Query(
        True,
        title="Veber filter",
        description="Check Veber filter criteria (true or false)",
    ),
    reos: bool = Query(
        True,
        title="REOS filter",
        description="Check REOS filter criteria (true or false)",
    ),
    ghose: bool = Query(
        True,
        title="Ghose filter",
        description="Check Ghose filter criteria (true or false)",
    ),
    ruleofthree: bool = Query(
        True,
        title="Rule of 3",
        description="Check Rule of 3 criteria (true or false)",
    ),
    qedscore: str = Query(
        "0-1",
        title="QED Druglikeliness",
        description="QED Score range (e.g., 0-1)",
    ),
    sascore: str = Query(
        "0-10",
        title="SAScore",
        description="Synthetic Accessibility Score range (e.g., 0-10)",
    ),
    nplikeness: str = Query(
        "-5-5",
        title="NPlikenessScore",
        description="Natural Product likeness Score range (e.g., -5-5)",
    ),
    filterOperator: Literal["AND", "OR"] = Query(
        "OR",
        title="Filter Operator",
        description="Logic for combining filter results: AND requires all selected filters to pass, OR requires at least one filter to pass",
    ),
):
    """
    Apply multiple chemical filters with detailed violation information.

    This endpoint provides comprehensive information about why molecules pass or fail
    specific filters, including:
    - PAINS substructure family and description when matched
    - Specific property values that cause Lipinski, Veber, REOS, Ghose, or Rule of 3 violations
    - Clear indication of whether finding a match is good or bad for drug-likeness

    Returns detailed results for each molecule with pass/fail status and violation details.
    """
    results = []

    for item in io.StringIO(smiles_list):
        smiles = item.strip()
        if not smiles:
            continue

        mol = parse_input(smiles, "rdkit", False)
        if not mol:
            results.append(
                {
                    "smiles": smiles,
                    "valid": False,
                    "error": "Invalid SMILES string",
                    "filters": {},
                }
            )
            continue

        molecule_result = {
            "smiles": smiles,
            "valid": True,
            "filters": {},
            "overall_pass": False,
            "filter_passes": [],
        }

        # Apply each requested filter with detailed information
        if pains:
            pains_result = get_PAINS_detailed(mol)
            molecule_result["filters"]["pains"] = pains_result
            molecule_result["filter_passes"].append(pains_result["passes"])

        if lipinski:
            lipinski_result = check_RO5_violations_detailed(mol)
            molecule_result["filters"]["lipinski"] = lipinski_result
            molecule_result["filter_passes"].append(lipinski_result["passes"])

        if veber:
            veber_result = get_VeberFilter_detailed(mol)
            molecule_result["filters"]["veber"] = veber_result
            molecule_result["filter_passes"].append(veber_result["passes"])

        if reos:
            reos_result = get_REOSFilter_detailed(mol)
            molecule_result["filters"]["reos"] = reos_result
            molecule_result["filter_passes"].append(reos_result["passes"])

        if ghose:
            ghose_result = get_GhoseFilter_detailed(mol)
            molecule_result["filters"]["ghose"] = ghose_result
            molecule_result["filter_passes"].append(ghose_result["passes"])

        if ruleofthree:
            ro3_result = get_RuleofThree_detailed(mol)
            molecule_result["filters"]["rule_of_three"] = ro3_result
            molecule_result["filter_passes"].append(ro3_result["passes"])

        # Handle numeric score filters
        if len(qedscore.split("-")) == 2:
            range_start, range_end = float(qedscore.split("-")[0]), float(
                qedscore.split("-")[1]
            )
            qed_value = QED.qed(mol)
            qed_passes = range_start <= qed_value <= range_end
            molecule_result["filters"]["qed"] = {
                "passes": qed_passes,
                "value": round(qed_value, 3),
                "range": f"{range_start}-{range_end}",
                "details": f"QED = {qed_value:.3f} ({'within' if qed_passes else 'outside'} range {range_start}-{range_end})",
            }
            molecule_result["filter_passes"].append(qed_passes)

        if len(sascore.split("-")) == 2:
            range_start, range_end = float(sascore.split("-")[0]), float(
                sascore.split("-")[1]
            )
            sas_value = get_sas_score(mol)
            sas_passes = range_start <= sas_value <= range_end
            molecule_result["filters"]["sa_score"] = {
                "passes": sas_passes,
                "value": sas_value,
                "range": f"{range_start}-{range_end}",
                "details": f"SA Score = {sas_value} ({'within' if sas_passes else 'outside'} range {range_start}-{range_end})",
            }
            molecule_result["filter_passes"].append(sas_passes)

        if len(nplikeness.split("-")) == 2:
            range_start, range_end = float(nplikeness.split("-")[0]), float(
                nplikeness.split("-")[1]
            )
            np_value = get_np_score(mol)
            np_passes = range_start <= float(np_value) <= range_end
            molecule_result["filters"]["np_likeness"] = {
                "passes": np_passes,
                "value": float(np_value),
                "range": f"{range_start}-{range_end}",
                "details": f"NP-likeness = {np_value} ({'within' if np_passes else 'outside'} range {range_start}-{range_end})",
            }
            molecule_result["filter_passes"].append(np_passes)

        # Determine overall pass based on filter operator
        if molecule_result["filter_passes"]:
            if filterOperator == "AND":
                molecule_result["overall_pass"] = all(molecule_result["filter_passes"])
            else:  # OR
                molecule_result["overall_pass"] = any(molecule_result["filter_passes"])

        results.append(molecule_result)

    # Filter results based on overall pass status
    passing_results = [
        r for r in results if r.get("overall_pass", False) or not r.get("valid", True)
    ]

    return {
        "total_molecules": len(results),
        "passing_molecules": len(passing_results),
        "filter_operator": filterOperator,
        "results": passing_results,
    }


@router.get(
    "/ertlfunctionalgroup",
    summary="using the algorithm proposed by Peter Ertl to identify functional groups",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateFunctionalGroupResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def get_functional_groups(
    smiles: str = Query(
        title="SMILES",
        description="SMILES string to be enumerated",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
):
    """For a given SMILES string this function generates a list of identified.

    functional groups.

    Parameters:
    - **SMILES**: required (query parameter): The SMILES string to be checked for functional groups.

    Returns:
    - List[str]: A list of identified functional groups, otherwise returns an error message.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    mol = parse_input(smiles, "rdkit", False)
    if mol:
        try:
            f_groups = get_ertl_functional_groups(mol)
            return f_groups
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))
    else:
        raise HTTPException(
            status_code=422,
            detail="Error reading SMILES string, please check again.",
        )


@router.get(
    "/standarizedTautomer",
    summary="Standardize tautomeric SMILES using RDKit EnumerateStereoisomers module",
    responses={
        200: {
            "description": "Successful response",
            "model": StandarizedTautomerResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def get_standardized_tautomer_smiles(
    smiles: str = Query(
        title="SMILES",
        description="SMILES string to be standardized",
        openapi_examples={
            "example1": {
                "summary": "Example: Caffeine",
                "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            },
            "example2": {
                "summary": "Example: Topiramate-13C6",
                "value": "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1",
            },
        },
    ),
):
    mol = parse_input(smiles, "rdkit", False)
    if mol:
        standardized_smiles = get_standardized_tautomer(mol)
        return standardized_smiles


@router.get(
    "/pubchem/smiles",
    summary="Retrieve canonical SMILES from PubChem for a given chemical identifier",
    responses={
        200: {"description": "Successful response", "model": PubChemResponse},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def get_pubchem_smiles(
    identifier: str = Query(
        title="Chemical Identifier",
        description="Chemical identifier (name, CID, InChI, InChIKey, SMILES, formula, CAS number)",
        openapi_examples={
            "example1": {
                "summary": "Example: Aspirin by name",
                "value": "aspirin",
            },
            "example2": {
                "summary": "Example: Ethanol by SMILES",
                "value": "CCO",
            },
            "example3": {
                "summary": "Example: Glucose by formula",
                "value": "C6H12O6",
            },
            "example4": {
                "summary": "Example: Water by CAS",
                "value": "7732-18-5",
            },
        },
    ),
):
    """
    Retrieve the canonical SMILES for a molecule from PubChem via the PUG REST API.

    The endpoint automatically detects the input type and handles:
      - CID: a string of digits
      - InChI: a string starting with "InChI="
      - InChIKey: matching the standard format (e.g., "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
      - CAS number: if the input matches a pattern like "7732-18-5"
      - Molecular formula: if the input matches a formula pattern (e.g., "C6H12O6")
      - SMILES: if the input contains no spaces, is short, and is composed of SMILES characters
      - Chemical name: default option (covers IUPAC names, synonyms, trivial names, etc.)

    Parameters:
        identifier (str): The chemical identifier to look up

    Returns:
        PubChemResponse: Object containing the input, canonical SMILES, detected input type, CIDs, PubChem links, and success status

    Raises:
        HTTPException: If the identifier is invalid or no results are found
    """
    if not identifier:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail="Chemical identifier is required",
        )

    try:
        # Get comprehensive compound information from PubChem
        result = pubchem_client.get_compound_info(identifier)

        # Create response with all available information
        return PubChemResponse(
            input=result["input"],
            canonical_smiles=result["canonical_smiles"],
            input_type=result["input_type"],
            success=result["success"],
            cids=result["cids"],
            pubchem_links=result["pubchem_links"],
        )

    except Exception as e:
        logger.error(f"Error processing PubChem request for '{identifier}': {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Error processing request: {str(e)}",
        )
