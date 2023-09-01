from fastapi import APIRouter, Query, status, HTTPException, Body
from typing import Optional, Literal, Annotated
from rdkit import Chem
from typing import Union
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
)
from chembl_structure_pipeline import standardizer, checker
from fastapi.responses import Response, JSONResponse
from app.modules.npscorer import getNPScore
from app.modules.classyfire import classify, result
from app.modules.toolkits.cdk_wrapper import (
    getTanimotoSimilarityCDK,
    getCDKHOSECodes,
)
from app.modules.toolkits.rdkit_wrapper import (
    getTanimotoSimilarityRDKit,
    getRDKitHOSECodes,
    getProperties,
)
from app.modules.coconut.descriptors import getCOCONUTDescriptors
from app.modules.all_descriptors import getTanimotoSimilarity
from app.modules.coconut.preprocess import getCOCONUTpreprocessing
import pandas as pd
from fastapi.templating import Jinja2Templates
from app.schemas import HealthCheck
from app.schemas.chemblstandardizer import (
    SMILESValidationResult,
    SMILESStandardizedResult,
    StandardizedResult,
)
from app.schemas.classyfire import ClassyFireJob, ClassyFireResult
from app.schemas.coconut import COCONUTPreprocessingModel
from app.schemas.error import ErrorResponse, BadRequestModel, NotFoundModel
from app.schemas.chem_schema import (
    GenerateStereoisomersResponse,
    GenerateDescriptorsResponse,
    GenerateMultipleDescriptorsResponse,
    GenerateHOSECodeResponse,
    GenerateStandardizeResponse,
    NPlikelinessScoreResponse,
    TanimotoSimilarityResponse,
)

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
    """
    ## Perform a Health Check
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
    )
):
    """
    For a given SMILES string this function enumerates all possible stereoisomers

    Parameters:
    - **SMILES**: required (query parameter): The SMILES string to be enumerated.

    Returns:
    - List[str]: A list of stereoisomer SMILES strings if successful, otherwise returns an error message.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        isomers = tuple(EnumerateStereoisomers(mol))
        smilesArray = []
        for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
            smilesArray.append(smi)
        return smilesArray
    else:
        raise HTTPException(
            status_code=400, detail="Error reading SMILES string, please check again."
        )


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
        default="json", description="Desired display format"
    ),
    toolkit: Literal["cdk", "rdkit", "all"] = Query(
        default="rdkit", description="Cheminformatics toolkit used in the backend"
    ),
):
    """
    Generates standard descriptors for the input molecule (SMILES).

    Parameters:
    - **SMILES**: required (query): The SMILES representation of the molecule.
    - **format**: optional (query): The desired format for the output.
        - Supported values: "html" / "json" (default), "json".
    - **toolkit**: optional (query): The chemical toolkit to use for descriptor calculation. Default is "rdkit". Allowed "all", "cdk".
        - Supported values: "cdk"/ "rdkit" / "all" (default), "rdkit".

    Returns:
    - If format is "html", returns an HTML response containing a table of descriptors and their values.
    - If format is not "html", returns the descriptors and their values in the specified format (default: JSON).

    Raises:
    - None

    """
    data = getCOCONUTDescriptors(smiles, toolkit)
    if format == "html":
        if toolkit == "all":
            headers = ["Descriptor name", "RDKit Descriptors", "CDK Descriptors"]
            df = pd.DataFrame.from_dict(
                data, orient="index", columns=headers[1:], dtype=object
            )
            df.insert(0, headers[0], df.index)
        else:
            headers = ["Descriptor name", "Values"]
            df = pd.DataFrame.from_dict(data, orient="index", columns=headers[1:])
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
        default="rdkit", description="Cheminformatics toolkit used in the backend"
    ),
):
    """
    Retrieve multiple descriptors for a list of SMILES strings.

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
            status_code=400, detail="At least two molecules are required."
        )

    descriptors_dict = {}

    for molecule in molecules:
        descriptors = getCOCONUTDescriptors(molecule, toolkit)
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
        default="rdkit", description="Cheminformatics toolkit used in the backend"
    ),
    ringsize: Optional[bool] = Query(
        False,
        title="ringsize",
        description="Determines whether to include information about ring sizes",
    ),
):
    """
    Generates HOSE codes for a given SMILES string.

    Parameters:
    - **SMILES**: required (query): The SMILES string represents the chemical compound.
    - **spheres**: required (query): The number of spheres to use for generating HOSE codes.
    - **toolkit**: optional (default:cdk): The chemical toolkit to use for generating HOSE codes.
            Supported values: "cdk" (default), "rdkit".
    - **ringsize**: optional (default:False): Determines whether to include information about ring sizes
            in the HOSE codes. Default is False.

    Returns:
    - List[str]: A list of HOSE codes if successful, indicating the HOSE codes
        for each atom in the molecule. Otherwise, returns an error message.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    """
    if toolkit == "cdk":
        hose_codes = await getCDKHOSECodes(smiles, spheres, ringsize)
    elif toolkit == "rdkit":
        hose_codes = await getRDKitHOSECodes(smiles, spheres)

    if hose_codes:
        return hose_codes
    else:
        raise HTTPException(
            status_code=500, detail="Failed to generate HOSE codes, check settings."
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
    ]
):
    """
    Standardize molblock using the ChEMBL curation pipeline
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
            raise HTTPException(status_code=400, detail="Invalid or missing molblock")

        if rdkit_mol:
            smiles = Chem.MolToSmiles(rdkit_mol, kekuleSmiles=True)
            response = dict(
                standardized_mol=standardized_mol,
                canonical_smiles=smiles,
                inchi=Chem.inchi.MolToInchi(rdkit_mol),
                inchikey=Chem.inchi.MolToInchiKey(rdkit_mol),
            )
            original_properties = getProperties(data)
            response.update(original_properties)
            return response

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


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
    fix: Optional[bool] = Query(
        False,
        title="Fix",
        description="Flag indicating whether to fix the issues by standardizing the SMILES.",
    ),
):
    """
    Check a given SMILES string and the represented structure for issues and standardize it using the ChEMBL curation pipeline.

    Parameters:
    - **SMILES**: required (str, query) The SMILES string to check and standardize.
    - **fix**: optional (bool): Flag indicating whether to fix the issues by standardizing the SMILES. Defaults to False.

    Returns:
    - If fix is False:
        - If issues are found in the SMILES string, return a list of issues.
        - If no issues are found, return the string "No Errors Found".

    - If fix is True:
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

                result = StandardizedResult(
                    original=SMILESValidationResult(smi=smiles, messages=issues),
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
            status_code=400, detail="Error reading SMILES string, please check again."
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
    )
):
    """
    Calculates the natural product likeness score based on the RDKit implementation.

    Parameters:
    - **SMILES**: required (query): The SMILES representation of the molecule.

    Returns:
    - np_score (float): The natural product likeness score.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    """
    try:
        np_score = getNPScore(smiles)
        if np_score:
            return float(np_score)
        else:
            raise HTTPException(
                status_code=400,
                detail="Error reading SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get(
    "/tanimoto",
    summary="Generates the Tanimoto similarity index for a given pair of SMILES strings",
    responses={
        200: {
            "description": "Successful response",
            "model": TanimotoSimilarityResponse,
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
        default="rdkit", description="Cheminformatics toolkit used in the backend"
    ),
):
    """
    Generates the Tanimoto similarity index for a given pair of SMILES strings.

    The Tanimoto similarity index is calculated using different toolkits:

    - When using the **'cdk'** toolkit (default), the Tanimoto similarity is calculated
      using the PubchemFingerprinter from the CDK library.
      More information: https://cdk.github.io/cdk/2.8/docs/api/org/openscience/cdk/fingerprint/PubchemFingerprinter.html

    - When using the **'rdkit'** toolkit, the Tanimoto similarity is calculated
      using Morgan fingerprints with a radius of 2 and nBits=1024.
      Additional modifications can be found in the rdkit_wrapper module.

    Usage:
    - To calculate the Tanimoto similarity for a pair of SMILES strings, provide them as a comma-separated parameter:
      Example: api.naturalproducts.net/latest/chem/tanimoto?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C,CN1C=NC2=C1C(=O)NC(=O)N2C

    Parameters:
    - **SMILES**: required (query): The SMILES strings for which the Tanimoto similarity will be calculated. Two SMILES strings should be provided as a comma-separated parameter.
    - **toolkit**: optional (defaults: cdk): The toolkit to use for Tanimoto calculation.
        - Supported values: "rdkit" / "cdk" (default), "rdkit".

    Returns:
    - The Tanimoto similarity index as a floating-point value.

    Raises:
    - If an error occurs during the Tanimoto similarity calculation or if the input SMILES strings are invalid, an appropriate error message will be returned.

    """
    molecules = [m.strip() for m in smiles.split(",")]

    if len(molecules) < 2:
        raise HTTPException(
            status_code=400, detail="At least two molecules are required."
        )

    elif len(smiles.split(",")) == 2:
        try:
            smiles1, smiles2 = smiles.split(",")
            if toolkit == "rdkit":
                Tanimoto = getTanimotoSimilarityRDKit(smiles1, smiles2)
            else:
                Tanimoto = getTanimotoSimilarityCDK(smiles1, smiles2)
            return Tanimoto
        except Exception:
            raise HTTPException(
                status_code=400,
                detail="Error reading SMILES string, please check again.",
            )

    elif len(smiles.split(",")) > 2:
        try:
            matrix = getTanimotoSimilarity(smiles, toolkit)
            return Response(content=matrix, media_type="text/html")
        except Exception:
            raise HTTPException(
                status_code=400,
                detail="Error reading SMILES string, please check again.",
            )
    else:
        raise HTTPException(
            status_code=400, detail="Error reading SMILES string, please check again."
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
    )
):
    """
    Generates an Input JSON file with information on different molecular representations and descriptors suitable for submission to the COCONUT database.

    Parameters:
    - **SMILES**: required (query): The SMILES string represents a chemical compound.

    Returns:
    - JSONResponse: The generated Input JSON file for COCONUT.

    Raises:
    - HTTPException: If there is an error reading the SMILES string.

    """
    try:
        data = getCOCONUTpreprocessing(smiles)
        if data:
            return JSONResponse(content=data)
        else:
            raise HTTPException(
                status_code=400,
                detail="Error reading SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(
            status_code=400, detail="Error processing request: " + str(e)
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
    )
):
    """
    Generate ClassyFire-based classifications using SMILES as input.

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
    try:
        classification_data = await classify(smiles)
        if classification_data:
            return classification_data
        else:
            raise HTTPException(
                status_code=400,
                detail="Error reading SMILES string, please check again.",
            )
    except Exception as e:
        raise HTTPException(
            status_code=500, detail="Error processing request: " + str(e)
        )


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
async def classyFire_result(jobid: str):
    """
    Retrieve the ClassyFire classification results based on the provided Job ID.
    To obtain the results from ClassyFire, please initiate a new request and obtain a unique job ID.
    Once you have the job ID, you need to submit another request using this ID in order to retrieve the desired outcome.

    Parameters:
    - **jobid**: required (query): The Job ID used to query the ClassyFire classification results.

    Raises:
    - HTTPException 400: If the Job ID is not provided.
    - HTTPException 500: If an error occurs during the classification process.

    Returns:
    - The ClassyFire classification results as JSON.

    """

    if jobid:
        try:
            # Replace with your function to retrieve result
            data = await result(jobid)
            return data
        except Exception as e:
            raise HTTPException(
                status_code=500, detail="Error processing request: " + str(e)
            )
    else:
        raise HTTPException(status_code=400, detail="Job ID is required.")
