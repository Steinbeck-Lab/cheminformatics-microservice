from __future__ import annotations

from typing import Literal

import selfies as sf
from fastapi import FastAPI
from fastapi import APIRouter
from fastapi import HTTPException
from fastapi import Query
from fastapi import status
from fastapi import Request
from fastapi import Body
from fastapi.responses import Response
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded
from rdkit import Chem

# Schema imports
from app.schemas import HealthCheck
from app.schemas.error import BadRequestModel
from app.schemas.error import ErrorResponse
from app.schemas.error import NotFoundModel
from app.schemas.converters_schema import GenerateCanonicalResponse
from app.schemas.converters_schema import GenerateCXSMILESResponse
from app.schemas.converters_schema import GenerateFormatsResponse
from app.schemas.converters_schema import GenerateInChIKeyResponse
from app.schemas.converters_schema import GenerateInChIResponse
from app.schemas.converters_schema import GenerateSELFIESResponse
from app.schemas.converters_schema import GenerateSMILESResponse
from app.schemas.converters_schema import ThreeDCoordinatesResponse
from app.schemas.converters_schema import TwoDCoordinatesResponse
from app.schemas.converters_schema import GenerateSMARTSResponse

# Module imports
from app.modules.toolkits.cdk_wrapper import get_canonical_SMILES
from app.modules.toolkits.cdk_wrapper import get_CDK_SDG_mol
from app.modules.toolkits.cdk_wrapper import get_CXSMILES
from app.modules.toolkits.cdk_wrapper import get_InChI
from app.modules.toolkits.cdk_wrapper import get_smiles_opsin
from app.modules.toolkits.helpers import parse_input
from app.modules.toolkits.openbabel_wrapper import get_ob_canonical_SMILES
from app.modules.toolkits.openbabel_wrapper import get_ob_InChI
from app.modules.toolkits.openbabel_wrapper import get_ob_mol
from app.modules.toolkits.rdkit_wrapper import get_2d_mol
from app.modules.toolkits.rdkit_wrapper import get_3d_conformers
from app.modules.toolkits.rdkit_wrapper import get_rdkit_CXSMILES

# Create the Limiter instance
limiter = Limiter(key_func=get_remote_address)

# Initialize FastAPI app
app = FastAPI()

# Add the middleware to handle rate limit exceeded errors
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

router = APIRouter(
    prefix="/convert",
    tags=["convert"],
    dependencies=[],
    responses={
        200: {"description": "OK"},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)


@router.get("/", include_in_schema=False)
@router.get(
    "/health",
    tags=["healthcheck"],
    summary="Perform a Health Check on Converters Module",
    response_description="Return HTTP Status Code 200 (OK)",
    status_code=status.HTTP_200_OK,
    response_model=HealthCheck,
    include_in_schema=False,
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
    "/mol2D",
    summary="Generates 2D Coordinates for the input molecules",
    responses={
        200: {
            "description": "Successful response",
            "model": TwoDCoordinatesResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
@limiter.limit("20/minute")
async def create2d_coordinates(
    request: Request,
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
    toolkit: Literal["cdk", "rdkit", "openbabel"] = Query(
        default="cdk",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Generates 2D Coordinates using the CDK Structure diagram.

    generator/RDKit/Open Babel and returns the mol block.

    Parameters:
    - **SMILES**: required (str): The SMILES string.
    - **toolkit** (str, optional): The toolkit to use for generating 2D coordinates.
        - Supported values: "cdk" (default), "rdkit", "openbabel".

    Returns:
    - molblock (str): The generated mol block with 2D coordinates as a plain text response.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """
    if toolkit == "cdk":
        mol = parse_input(smiles, "cdk", False)
        return Response(
            content=get_CDK_SDG_mol(mol).replace("$$$$\n", ""),
            media_type="text/plain",
        )
    elif toolkit == "rdkit":
        mol = parse_input(smiles, "rdkit", False)
        return Response(
            content=get_2d_mol(mol),
            media_type="text/plain",
        )
    else:
        mol = parse_input(smiles, "rdkit", False)
        if mol:
            return Response(
                content=get_ob_mol(smiles),
                media_type="text/plain",
            )


@router.get(
    "/mol3D",
    summary="Generates 3D Coordinates for the input molecules",
    responses={
        200: {
            "description": "Successful response",
            "model": ThreeDCoordinatesResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
@limiter.limit("20/minute")
async def create3d_coordinates(
    request: Request,
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
    toolkit: Literal["rdkit", "openbabel"] = Query(
        default="openbabel",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Generates a random 3D conformer from SMILES using the specified molecule.

    toolkit.

    Parameters:
    - **SMILES**: required (str): The SMILES representation of the molecule.
    - **toolkit**: optional (str): The molecule toolkit to use.
        - Supported values: "rdkit"  & "openbabel" (default).

    Returns:
    - molblock (str): The generated mol block with 3D coordinates as a plain text response.

    Raises:
    - ValueError: If the SMILES string is not provided or is invalid.
    """

    if toolkit == "rdkit":
        mol = parse_input(smiles, "rdkit", False)
        return Response(
            content=get_3d_conformers(mol, depict=False),
            media_type="text/plain",
        )
    elif toolkit == "openbabel":
        mol = parse_input(smiles, "rdkit", False)
        if mol:
            return Response(
                content=get_ob_mol(smiles, threeD=True),
                media_type="text/plain",
            )


@router.get(
    "/smiles",
    summary="Generate SMILES from a given input",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateSMILESResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
@limiter.limit("10/minute")
async def iupac_name_or_selfies_to_smiles(
    request: Request,
    input_text: str = Query(
        title="Input IUPAC name or SELFIES",
        description="IUPAC name or SELFIES representation of the molecule",
        openapi_examples={
            "example1": {
                "summary": "Example: IUPAC name",
                "value": "1,3,7-trimethylpurine-2,6-dione",
            },
            "example2": {
                "summary": "Example: SELFIES",
                "value": "[C][N][C][=N][C][=C][Ring1][Branch1][C][=Branch1][C][=O][N][Branch1][=Branch2][C][=Branch1][C][=O][N][Ring1][Branch2][C][C]",
            },
        },
    ),
    representation: Literal["iupac", "selfies"] = Query(
        default="iupac",
        description="Required type of format conversion",
    ),
    converter: Literal["opsin"] = Query(
        default="opsin",
        description="Required type of converter for IUPAC",
    ),
):
    """Generate SMILES from a given IUPAC name or a SELFIES representation.

    Parameters:
    - **input_text**: required (str): The input text containing either the IUPAC name or SELFIES representation.
    - **representation**: optional (str): The representation type of the input text.
        - Supported values: "iupac" (default) & "selfies".

    Returns:
    - If representation is "iupac": The generated SMILES string corresponding to the given IUPAC name.
    - If the representation is "selfies": The generated SMILES string corresponds to the given SELFIES representation.

    Notes:
    - The IUPAC name should follow the standard IUPAC naming conventions for organic compounds.
    - SELFIES (Self-Referencing Embedded Strings) is a concise yet expressive chemical string notation.

    Example Usage:
    - To generate SMILES from an IUPAC name: /smiles?input_text=benzene&representation=iupac
    - To generate SMILES from a SELFIES representation: /smiles?input_text=[C][C][C]&representation=selfies
    """
    try:
        if representation == "iupac":
            if converter == "opsin":
                iupac_name = get_smiles_opsin(input_text)
            if iupac_name:
                return str(iupac_name)
        elif representation == "selfies":
            selfies_out = sf.decoder(input_text)
            if selfies_out:
                return str(selfies_out)
        else:
            raise HTTPException(
                status_code=422,
                detail="Error reading input text, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))


@router.get(
    "/canonicalsmiles",
    summary="Generate CanonicalSMILES from a given SMILES",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateCanonicalResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def smiles_canonicalise(
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
    toolkit: Literal["cdk", "rdkit", "openbabel"] = Query(
        default="cdk",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Canonicalizes a given SMILES string according to the allowed toolkits.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to be canonicalized.
    - **toolkit**: optional (str): The toolkit to use for canonicalization.
        - Supported values: "cdk" (default), "rdkit" & "openbabel".

    Returns:
    - SMILES (str): The canonicalized SMILES string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.
    """
    if toolkit == "cdk":
        mol = parse_input(smiles, "cdk", False)
        return str(get_canonical_SMILES(mol))
    elif toolkit == "rdkit":
        mol = parse_input(smiles, "rdkit", False)
        return str(Chem.MolToSmiles(mol, kekuleSmiles=True))
    elif toolkit == "openbabel":
        smiles = get_ob_canonical_SMILES(smiles)
        return smiles


@router.get(
    "/cxsmiles",
    summary="Generate CXSMILES from a given SMILES",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateCXSMILESResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def smiles_to_cxsmiles(
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
    toolkit: Literal["cdk", "rdkit"] = Query(
        default="cdk",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Convert SMILES to CXSMILES.

    For more information:
    - https://docs.chemaxon.com/display/docs/chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.md

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default) & "rdkit".

    Returns:
    - CXSMILES (str): The converted CXSMILES string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.

    Note:
    - CXSMILES is a Chemaxon Extended SMILES which is used for storing special features of the molecules after the SMILES string.
    """
    if toolkit == "cdk":
        mol = parse_input(smiles, "cdk", False)
        cxsmiles = get_CXSMILES(mol)
        if cxsmiles:
            return str(cxsmiles)
    else:
        mol = parse_input(smiles, "rdkit", False)
        cxsmiles = get_rdkit_CXSMILES(mol)
        if cxsmiles:
            return str(cxsmiles)


@router.get(
    "/inchi",
    summary="Generate InChI from a given SMILES",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateInChIResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def smiles_to_inchi(
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
    toolkit: Literal["cdk", "rdkit", "openbabel"] = Query(
        default="cdk",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Convert SMILES to InChI.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default), "openbabel" & "rdkit".

    Returns:
    - InChI (str): The resulting InChI string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.
    """
    if toolkit == "cdk":
        mol = parse_input(smiles, "cdk", False)
        inchi = get_InChI(mol)
        if inchi:
            return str(inchi)
    elif toolkit == "rdkit":
        mol = parse_input(smiles, "rdkit", False)
        if mol:
            inchi = Chem.inchi.MolToInchi(mol)
            if inchi:
                return str(inchi)
    elif toolkit == "openbabel":
        inchi = get_ob_InChI(smiles)
        if inchi:
            return str(inchi)


@router.get(
    "/inchikey",
    summary="Generate InChI-Key from a given SMILES",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateInChIKeyResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def smiles_to_inchikey(
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
    toolkit: Literal["cdk", "rdkit", "openbabel"] = Query(
        default="cdk",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Convert SMILES to InChI-Key.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default), "openbabel" & "rdkit".

    Returns:
    - InChI-Key (str): The resulting InChI-Key string.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.
    """
    if toolkit == "cdk":
        mol = parse_input(smiles, "cdk", False)
        inchikey = get_InChI(mol, InChIKey=True)
        if inchikey:
            return str(inchikey)

    elif toolkit == "rdkit":
        mol = parse_input(smiles, "rdkit", False)
        if mol:
            inchikey = Chem.inchi.MolToInchiKey(mol)
            if inchikey:
                return str(inchikey)
    elif toolkit == "openbabel":
        inchikey = get_ob_InChI(smiles, InChIKey=True)
        if inchikey:
            return str(inchikey)


@router.get(
    "/selfies",
    summary="Generates SELFIES string for a given SMILES string",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateSELFIESResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def encode_selfies(
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
):
    """Generates SELFIES string for a given SMILES string.

    For more information:
    - Krenn et al, SELFIES and the future of molecular string representations, Patterns, https://doi.org/10.1016/j.patter.2022.100588.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.

    Returns:
    - SELFIES (str): The resulting SELFIES of the chemical compound.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    """
    try:
        selfies_e = sf.encoder(smiles)
        if selfies_e:
            return str(selfies_e)
        else:
            raise HTTPException(
                status_code=400,
                detail="Error reading input text, please check again.",
            )
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get(
    "/formats",
    summary="Convert SMILES to various molecular formats using different toolkits",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateFormatsResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def smiles_convert_to_formats(
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
    toolkit: Literal["cdk", "rdkit", "openbabel"] = Query(
        default="cdk",
        description="Cheminformatics toolkit used in the backend",
    ),
):
    """Convert SMILES to various molecular formats using different toolkits.

    Parameters:
    - **SMILES**: required (str): The input SMILES string to convert.
    - **toolkit**: optional (str): The toolkit to use for conversion.
        - Supported values: "cdk" (default), "openbabel" & "rdkit".

    Returns:
    - dict: A dictionary containing the converted data in various formats. The dictionary has the following keys:
        - "mol" (str): The generated 2D mol block of the molecule.
        - "canonicalsmiles" (str): The canonical SMILES representation of the molecule.
        - "inchi" (str): The InChI representation of the molecule.
        - "inchikey" (str): The InChIKey representation of the molecule.

    Note:
        - The returned dictionary may contain empty strings if conversion fails or the input SMILES string is invalid.

    Raises:
    - ValueError: If the SMILES string is empty or contains invalid characters.
    - ValueError: If an unsupported toolkit option is provided.
    """
    try:
        if toolkit == "cdk":
            response = {}
            mol = parse_input(smiles, "cdk", False)
            response["mol"] = get_CDK_SDG_mol(mol).replace("$$$$\n", "")
            response["canonicalsmiles"] = str(get_canonical_SMILES(mol))
            response["inchi"] = str(get_InChI(mol))
            response["inchikey"] = str(get_InChI(mol, InChIKey=True))
            return response

        elif toolkit == "rdkit":
            mol = parse_input(smiles, "rdkit", False)
            if mol:
                response = {}
                response["mol"] = Chem.MolToMolBlock(mol)
                response["canonicalsmiles"] = Chem.MolToSmiles(
                    mol,
                    kekuleSmiles=True,
                )
                response["inchi"] = Chem.inchi.MolToInchi(mol)
                response["inchikey"] = Chem.inchi.MolToInchiKey(mol)
                return response
        elif toolkit == "openbabel":
            response = {}
            response["mol"] = get_ob_mol(smiles)
            response["canonicalsmiles"] = get_ob_canonical_SMILES(smiles)
            response["inchi"] = get_ob_InChI(smiles)
            response["inchikey"] = get_ob_InChI(smiles, InChIKey=True)
            return response
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
    "/smarts",
    summary="Generate SMARTS from a given SMILES",
    responses={
        200: {
            "description": "Successful response",
            "model": GenerateSMARTSResponse,
        },
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
async def smiles_to_smarts(
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
    toolkit: Literal["rdkit"] = Query(
        default="rdkit",
        description="Cheminformatics toolkit used in the backend",
    ),
):

    if toolkit == "rdkit":
        mol = parse_input(smiles, "rdkit", False)
        if mol:
            smarts = Chem.MolToSmarts(mol)
            if smarts:
                return str(smarts)


@router.post(
    "/batch",
    summary="Batch convert chemical structures to various formats",
    responses={
        200: {"description": "Successful response"},
        400: {"description": "Bad Request", "model": BadRequestModel},
        404: {"description": "Not Found", "model": NotFoundModel},
        422: {"description": "Unprocessable Entity", "model": ErrorResponse},
    },
)
@limiter.limit("10/minute")
async def batch_convert(
    request: Request,
    body: dict = Body(...),
    output_format: str = Query(
        default="smiles",
        description="Format to convert to (smiles, canonicalsmiles, inchi, inchikey, selfies, cxsmiles, smarts, mol2d, mol3d)",
    ),
    toolkit: Literal["cdk", "rdkit", "openbabel"] = Query(
        default="cdk",
        description="Cheminformatics toolkit to use for conversion",
    ),
):
    """Batch convert chemical structures to various formats.

    This endpoint accepts a list of inputs with different formats and converts them
    to the specified output format using the selected toolkit.

    Parameters:
    - **body**: required (dict): JSON object with a list of inputs to convert.
        - Structure:
            ```json
            {
                "inputs": [
                    {
                        "value": "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
                        "input_format": "smiles"
                    }
                ]
            }
            ```
    - **output_format**: optional (str): Format to convert to.
        - Supported values: "smiles", "canonicalsmiles", "inchi", "inchikey", "selfies", "cxsmiles", "smarts", "mol2d", "mol3d".
    - **toolkit**: optional (str): Toolkit to use for conversion.
        - Supported values: "cdk" (default), "rdkit", "openbabel".

    Returns:
    - JSON object containing conversion results and summary.

    Note:
    - Some conversion combinations may not be supported by all toolkits.
    - Failed conversions will be included in the response with error messages.
    """
    results = []
    success_count = 0
    failure_count = 0

    # Convert dict to our expected format
    inputs = []
    if "inputs" in body:
        inputs = body["inputs"]

    for input_item in inputs:
        try:
            # Extract values from the input dictionary
            value = input_item.get("value", "")
            input_format = input_item.get("input_format", "")

            if not value or not input_format:
                raise ValueError("Missing required fields: value or input_format")

            # First convert input to SMILES if it's not already in SMILES format
            smiles = value

            if input_format.lower() == "iupac":
                smiles = get_smiles_opsin(value)
                if not smiles:
                    raise ValueError(
                        f"Failed to convert IUPAC name '{value}' to SMILES"
                    )
            elif input_format.lower() == "selfies":
                smiles = sf.decoder(value)
                if not smiles:
                    raise ValueError(f"Failed to decode SELFIES '{value}' to SMILES")
            elif input_format.lower() == "inchi":
                # Use RDKit to convert InChI to SMILES
                mol = Chem.inchi.MolFromInchi(value)
                if not mol:
                    raise ValueError(f"Failed to convert InChI '{value}' to molecule")
                smiles = Chem.MolToSmiles(mol)
            elif input_format.lower() != "smiles":
                raise ValueError(f"Unsupported input format: {input_format}")

            # Now convert SMILES to the desired output format
            output_value = ""

            if output_format.lower() == "smiles":
                output_value = smiles

            elif output_format.lower() == "canonicalsmiles":
                if toolkit == "cdk":
                    mol = parse_input(smiles, "cdk", False)
                    output_value = str(get_canonical_SMILES(mol))
                elif toolkit == "rdkit":
                    mol = parse_input(smiles, "rdkit", False)
                    output_value = str(Chem.MolToSmiles(mol, kekuleSmiles=True))
                elif toolkit == "openbabel":
                    output_value = get_ob_canonical_SMILES(smiles)

            elif output_format.lower() == "inchi":
                if toolkit == "cdk":
                    mol = parse_input(smiles, "cdk", False)
                    output_value = str(get_InChI(mol))
                elif toolkit == "rdkit":
                    mol = parse_input(smiles, "rdkit", False)
                    output_value = str(Chem.inchi.MolToInchi(mol))
                elif toolkit == "openbabel":
                    output_value = get_ob_InChI(smiles)

            elif output_format.lower() == "inchikey":
                if toolkit == "cdk":
                    mol = parse_input(smiles, "cdk", False)
                    output_value = str(get_InChI(mol, InChIKey=True))
                elif toolkit == "rdkit":
                    mol = parse_input(smiles, "rdkit", False)
                    output_value = str(Chem.inchi.MolToInchiKey(mol))
                elif toolkit == "openbabel":
                    output_value = get_ob_InChI(smiles, InChIKey=True)

            elif output_format.lower() == "selfies":
                output_value = str(sf.encoder(smiles))

            elif output_format.lower() == "cxsmiles":
                if toolkit == "cdk":
                    mol = parse_input(smiles, "cdk", False)
                    output_value = str(get_CXSMILES(mol))
                elif toolkit == "rdkit":
                    mol = parse_input(smiles, "rdkit", False)
                    output_value = str(get_rdkit_CXSMILES(mol))
                else:
                    raise ValueError(
                        f"CXSMILES conversion not supported by toolkit: {toolkit}"
                    )

            elif output_format.lower() == "smarts":
                if toolkit == "rdkit":
                    mol = parse_input(smiles, "rdkit", False)
                    output_value = str(Chem.MolToSmarts(mol))
                else:
                    raise ValueError(
                        f"SMARTS conversion not supported by toolkit: {toolkit}"
                    )

            elif output_format.lower() == "mol2d":
                if toolkit == "cdk":
                    mol = parse_input(smiles, "cdk", False)
                    output_value = get_CDK_SDG_mol(mol).replace("$$$$\n", "")
                elif toolkit == "rdkit":
                    mol = parse_input(smiles, "rdkit", False)
                    output_value = get_2d_mol(mol)
                elif toolkit == "openbabel":
                    output_value = get_ob_mol(smiles)

            elif output_format.lower() == "mol3d":
                if toolkit == "rdkit":
                    mol = parse_input(smiles, "rdkit", False)
                    output_value = get_3d_conformers(mol, depict=False)
                elif toolkit == "openbabel":
                    output_value = get_ob_mol(smiles, threeD=True)
                else:
                    raise ValueError(
                        f"3D coordinates generation not supported by toolkit: {toolkit}"
                    )

            else:
                raise ValueError(f"Unsupported output format: {output_format}")

            # Create a result dictionary
            results.append(
                {
                    "input": {"value": value, "input_format": input_format},
                    "output": output_value,
                    "success": True,
                    "error": "",
                }
            )
            success_count += 1

        except Exception as e:
            # Create an error result dictionary
            results.append(
                {
                    "input": {
                        "value": input_item.get("value", ""),
                        "input_format": input_item.get("input_format", ""),
                    },
                    "output": "",
                    "success": False,
                    "error": str(e),
                }
            )
            failure_count += 1

    summary = {
        "total": len(inputs),
        "successful": success_count,
        "failed": failure_count,
    }

    # Return the response as a dictionary
    return {"results": results, "summary": summary}
