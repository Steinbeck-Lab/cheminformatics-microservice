from fastapi import APIRouter, Query, status, HTTPException
from typing import Optional, Literal, List, Union, Dict
from rdkit import Chem
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
)
from app.modules.coconut.descriptors import getCOCONUTDescriptors
from app.modules.alldescriptors import getTanimotoSimilarity
from app.modules.coconut.preprocess import getCOCONUTpreprocessing
import pandas as pd
from fastapi.templating import Jinja2Templates
from app.schemas import HealthCheck
from app.schemas.pydanticmodels import (
    ErrorResponse,
    StandardizeRequest,
    StandardizeResponse,
    SMILESValidationResult,
    StandardizedResult,
    NPlikelinessScoreResponse,
)

router = APIRouter(
    prefix="/chem",
    tags=["chem"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
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
    Endpoint to perform a healthcheck on. This endpoint can primarily be used Docker
    to ensure a robust container orchestration and management is in place. Other
    services which rely on proper functioning of the API service will not deploy if this
    endpoint returns any other HTTP status code except 200 (OK).
    Returns:
        HealthCheck: Returns a JSON response with the health status
    """
    return HealthCheck(status="OK")


@router.get(
    "/stereoisomers",
    response_model=Union[List[str], str],
    summary="Enumerate all possible stereoisomers",
    responses={400: {"model": ErrorResponse}},
)
async def get_stereoisomers(
    smiles: str = Query(
        title="SMILES",
        description="SMILES string to be converted",
        examples=[
            "CCO",
            "C=O",
        ],
    )
):
    """
    For a given SMILES string this function enumerates all possible stereoisomers

    Parameters:
    - **SMILES**: required (query parameter): The SMILES string to be converted.

    Returns:
    - Union[List[str], str]: A list of stereo isomer SMILES strings if successful, otherwise returns an error message.

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
    responses={400: {"model": ErrorResponse}},
)
async def get_descriptors(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecule",
        examples=[
            "CCO",
            "C=O",
        ],
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
        return data


@router.get(
    "/descriptors/multiple",
    response_model=Union[Dict, Dict],
    summary="Generates descriptors for the input molecules",
    responses={400: {"model": ErrorResponse}},
)
async def calculate_descriptors(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecules",
        examples=[
            "CCO,CCC",
            "C=O,CC=O",
        ],
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
    - Union[Dict[str, Any], str]: If multiple SMILES are provided, returns a dictionary with each SMILES as the key and the corresponding descriptors as the value. If only one SMILES is provided, returns an error message.

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
        print(type(descriptors))
        descriptors_dict[molecule] = descriptors

    return descriptors_dict


@router.get(
    "/HOSEcode",
    response_model=List,
    summary="Generates HOSE codes for the input molecules",
    responses={400: {"model": ErrorResponse}},
)
async def HOSE_Codes(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecule",
        examples=[
            "CCO",
            "C=O",
        ],
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
    - **SMILES**: required (query): The SMILES string representing the chemical compound.
    - **spheres**: required (query): The number of spheres to use for generating HOSE codes.
    - **toolkit**: optional (default:cdk): The chemical toolkit to use for generating HOSE codes.
            Supported values: "cdk" (default), "rdkit".
    - **ringsize**: optional (default:False): Determines whether to include information about ring sizes
            in the HOSE codes. Default is False.

    Returns:
    - Union[List[str], str]: A list of HOSE codes if successful, indicating the HOSE codes
        for each atom in the molecule. Otherwise, returns an error message.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    """
    if smiles:
        if toolkit == "cdk":
            hose_codes = await getCDKHOSECodes(smiles, spheres, ringsize)
        elif toolkit == "rdkit":
            hose_codes = await getRDKitHOSECodes(smiles, spheres)
    else:
        raise HTTPException(
            status_code=400, detail="Error reading SMILES string, please check again."
        )

    if hose_codes:
        return hose_codes
    else:
        raise HTTPException(
            status_code=500, detail="Failed to generate HOSE codes, check settings."
        )


@router.post(
    "/standardize",
    response_model=StandardizeResponse,
    summary="Standardize molblock using the ChEMBL curation pipeline",
    responses={400: {"model": ErrorResponse}},
)
async def standardize_mol(request_data: StandardizeRequest):
    """
    Standardize molblock using the ChEMBL curation pipeline
    and return the standardized molecule, SMILES, InChI, and InCHI-Key.

    Parameters:
    - **request_data**: The request body containing the "molblock" string representing the molecule to be standardized.

    Returns:
    - dict: A dictionary containing the following keys:
        - "standardized_mol" (str): The standardized molblock of the molecule.
        - "canonical_smiles" (str): The canonical SMILES representation of the molecule.
        - "inchi" (str): The InChI representation of the molecule.
        - "inchikey" (str): The InChI-Key of the molecule.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.

    """
    molblock = request_data.molblock

    if molblock:
        standardized_mol = standardizer.standardize_molblock(molblock)
        rdkit_mol = Chem.MolFromMolBlock(standardized_mol)

    else:
        raise HTTPException(status_code=400, detail="Invalid or missing molblock")

    if rdkit_mol:
        smiles = Chem.MolToSmiles(rdkit_mol, kekuleSmiles=True)
        response = StandardizeResponse(
            standardized_mol=standardized_mol,
            canonical_smiles=smiles,
            inchi=Chem.inchi.MolToInchi(rdkit_mol),
            inchikey=Chem.inchi.MolToInchiKey(rdkit_mol),
        )
        return response
    else:
        raise HTTPException(status_code=400, detail="Invalid molblock format")


@router.get(
    "/errors",
    summary="Check a given SMILES string and the represented structure for issues and standardizes it",
    responses={400: {"model": ErrorResponse}},
)
async def check_errors(
    smiles: str = Query(
        title="SMILES",
        description="The SMILES string to check and standardize.",
        examples=[
            "CCO",
            "C=O",
        ],
    ),
    fix: Optional[bool] = Query(
        False,
        title="Fix",
        description="Flag indicating whether to fix the issues by standardizing the SMILES.",
    ),
):
    """
    Check a given SMILES string and the represented structure for issues and standardizes it using the ChEMBL curation pipeline.

    Parameters:
    - **SMILES**: required (str,query) The SMILES string to check and standardize.
    - **fix**: optional (bool): Flag indicating whether to fix the issues by standardizing the SMILES. Defaults to False.

    Returns:
    - If fix is False:
        - If issues are found in the SMILES string, returns a list of issues.
        - If no issues are found, returns the string "No Errors Found".

    - If fix is True:
        - If issues are found in the SMILES string, returns a dictionary containing the original SMILES, original issues,
          standardized SMILES, and new issues after standardization.
        - If no issues are found after standardization, returns a dictionary with the original SMILES and "No Errors Found".

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
    response_model=NPlikelinessScoreResponse,
    summary="Generates descriptors for the input molecules",
    responses={400: {"model": ErrorResponse}},
)
async def NPlikeliness_Score(
    smiles: str = Query(
        title="SMILES",
        description="The SMILES string to calculate the natural product likeness score",
        examples=[
            "CCO",
            "C=O",
        ],
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
            return {"np_score": float(np_score)}
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
    responses={400: {"model": ErrorResponse}},
)
async def Tanimoto_Similarity(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the molecules",
        examples=[
            "CCO,CCC",
            "C=O,CC=O",
        ],
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
    response_model=Dict,
    summary="Generates an Input JSON file with information for COCONUT database",
    responses={400: {"model": ErrorResponse}},
)
async def COCONUT_Preprocessing(
    smiles: str = Query(
        title="SMILES",
        description="SMILES string representing a chemical compound",
        examples=[
            "CCO",
            "C=O",
        ],
    )
):
    """
    Generates an Input JSON file with information of different molecular representations and descriptors suitable for submission to COCONUT database.

    Parameters:
    - **SMILES**: required (query): The SMILES string representing a chemical compound.

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
    response_model=Dict,
    summary="Generate ClassyFire-based classifications using SMILES as input",
    responses={400: {"model": ErrorResponse}},
)
async def ClassyFire_Classify(
    smiles: str = Query(
        title="SMILES",
        description="SMILES representation of the compound to be classified",
        examples=[
            "CCO",
            "C=O",
        ],
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
        print(e)
        raise HTTPException(
            status_code=500, detail="Error during classification process"
        )


@router.get(
    "/classyfire/{jobid}/result",
    response_model=Dict,
    summary="Retrieve the ClassyFire classification results based on the provided Job ID",
    responses={400: {"model": ErrorResponse}},
)
async def ClassyFire_result(jobid: str):
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
            data = await result(jobid)  # Replace with your function to retrieve result
            return data
        except Exception as e:
            print(e)
            raise HTTPException(
                status_code=500, detail="Error during classification process."
            )
    else:
        raise HTTPException(status_code=400, detail="Job ID is required.")
