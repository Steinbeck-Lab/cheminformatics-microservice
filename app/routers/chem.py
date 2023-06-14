from fastapi import Body, APIRouter
from typing import Optional
from typing_extensions import Annotated
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
from app.modules.coconut.preprocess import COCONUTpreprocessing
import pandas as pd
from fastapi.templating import Jinja2Templates

router = APIRouter(
    prefix="/chem",
    tags=["chem"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)

templates = Jinja2Templates(directory="app/templates")


@router.get("/")
async def chem_index():
    return {"module": "chem/test", "message": "Successful", "status": 200}


@router.get("/stereoisomers")
async def SMILES_to_Stereo_Isomers(smiles: str):
    """
    Enumerate all possible stereoisomers based on the chiral centres in the given SMILES:

    - **SMILES**: required (query parameter)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        isomers = tuple(EnumerateStereoisomers(mol))
        smilesArray = []
        for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
            smilesArray.append(smi)
        return smilesArray
    else:
        return "Error reading SMILES string, check again."


@router.get("/descriptors")
async def SMILES_Descriptors(
    smiles: str, format: Optional[str] = "json", toolkit: Optional[str] = "rdkit"
):
    """
    Generate standard descriptors for the input molecules (SMILES):

    - **SMILES**: required (query)
    """
    if smiles:
        if format == "html":
            data = getCOCONUTDescriptors(smiles, toolkit)
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
            return getCOCONUTDescriptors(smiles, toolkit)


@router.get("/descriptors/multiple")
async def SMILES_Descriptors_multiple(smiles: str, toolkit: Optional[str] = "rdkit"):
    """
    Generate standard descriptors for the input molecules (SMILES):

    - **SMILES**: required (query)
    - **toolkit**: Optional (default:rdkit, allowed:cdk)
    """
    if len(smiles.split(",")) > 1:
        combinedDict = {
            entry: getCOCONUTDescriptors(entry, toolkit) for entry in smiles.split(",")
        }
        return combinedDict
    else:
        return "Error invalid SMILES"


@router.get("/HOSEcode")
async def HOSE_Codes(
    smiles: str,
    spheres: int,
    toolkit: Optional[str] = "cdk",
    ringsize: Optional[bool] = False,
):
    """
    Generates HOSE Codes using CDK/RDKit.

    - **SMILES**: required (query)
    - **spheres**: required (query)
    - **toolkit**: Optional (default:cdk)
    - **ringsize**: Optional (default:False)
    """
    if smiles:
        if toolkit == "cdk":
            return await getCDKHOSECodes(smiles, spheres, ringsize)
        elif toolkit == "rdkit":
            return await getRDKitHOSECodes(smiles, spheres)
    else:
        return "Error reading SMILES string, check again."


@router.post("/standardize")
async def Standardize_Mol(mol: Annotated[str, Body(embed=True)]):
    """
    Standardize molblock using the ChEMBL curation pipeline routine
    and return the Standardized molecule, SMILES, InChI and InCHI-Key:

    - **mol**: required
    """
    if mol:
        standardized_mol = standardizer.standardize_molblock(mol)
        rdkit_mol = Chem.MolFromMolBlock(standardized_mol)
        smiles = Chem.MolToSmiles(rdkit_mol, kekuleSmiles=True)
        response = {}
        response["standardized_mol"] = standardized_mol
        response["cannonical_smiles"] = smiles
        response["inchi"] = Chem.inchi.MolToInchi(rdkit_mol)
        response["inchikey"] = Chem.inchi.MolToInchiKey(rdkit_mol)
        return response
    else:
        return "Error reading SMILES string, check again."


@router.get("/errors")
async def Check_Errors(smiles: str, fix: Optional[bool] = False):
    """
    Check issues for a given SMILES string and standardize it using the ChEMBL curation pipeline.

    - **SMILES**: required (query)
    - **fix**: optional (defaults: False)
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    if smiles:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            mol_block = Chem.MolToMolBlock(mol)
            if len(checker.check_molblock(mol_block)) == 0:
                return "No Errors Found"
            else:
                issues = checker.check_molblock(mol_block)
                if fix:
                    issues = checker.check_molblock(mol_block)
                    standardized_mol = standardizer.standardize_molblock(mol_block)
                    issues_new = checker.check_molblock(standardized_mol)
                    rdkit_mol = Chem.MolFromMolBlock(standardized_mol)
                    standardizedsmiles = Chem.MolToSmiles(rdkit_mol)
                    if len(issues_new) == 0:
                        issues_new = "No Errors Found"

                    parsed_data = {
                        "source": {
                            "SMILES": smiles,
                            "messages": issues,
                        },
                        "standardized": {
                            "SMILES": standardizedsmiles,
                            "messages": issues_new,
                        },
                    }
                    return parsed_data
                else:
                    return issues
        else:
            return "Error reading SMILES string, check again."
    else:
        return "Error reading SMILES string, check again."


@router.get("/nplikeness/score")
async def NPlikeliness_Score(smiles: str):
    """
    Generate natural product likeliness score based on RDKit implementation

    - **SMILES**: required (query)
    """
    if smiles:
        np_score = getNPScore(smiles)
        return np_score


@router.get("/tanimoto")
async def Tanimoto_Similarity(smiles: str, toolkit: Optional[str] = "cdk"):
    """
    Generate the Tanimoto similarity index for a given pair of SMILES strings.

    - **SMILES**: required (query)
    - **toolkit**: optional (defaults: cdk)
    """
    if len(smiles.split(",")) == 2:
        try:
            smiles1, smiles2 = smiles.split(",")
            if toolkit == "rdkit":
                Tanimoto = getTanimotoSimilarityRDKit(smiles1, smiles2)
            else:
                Tanimoto = getTanimotoSimilarityCDK(smiles1, smiles2)
            return Tanimoto
        except ValueError:
            return 'Please give a SMILES pair with "," separated. (Example: api.naturalproducts.net/chem/tanimoto?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C,CN1C=NC2=C1C(=O)NC(=O)N2C)'
    elif len(smiles.split(",")) > 2:
        try:
            matrix = getTanimotoSimilarity(smiles, toolkit)
            return Response(content=matrix, media_type="text/html")
        except ValueError:
            return 'Please give a SMILES pair with "," separated. (Example: api.naturalproducts.net/chem/tanimoto?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C,CN1C=NC2=C1C(=O)NC(=O)N2C)'
    else:
        return 'Please give a SMILES pair with "," separated. (Example: api.naturalproducts.net/chem/tanimoto?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C,CN1C=NC2=C1C(=O)NC(=O)N2C)'


@router.get("/coconut/pre-processing")
async def COCONUT_Preprocessing(smiles: str):
    """
    Generate Input JSON file for COCONUT.

    - **SMILES**: required (query)
    """
    if smiles:
        data = COCONUTpreprocessing(smiles)
        return JSONResponse(content=data)
    else:
        return "Error reading SMILES string, check again."


@router.get("/classyfire/classify")
async def ClassyFire_Classify(smiles: str):
    """
    Generate ClassyFire-based classifications using SMILES as input.

    - **SMILES**: required (query)
    """
    if smiles:
        data = await classify(smiles)
        return data


@router.get("/classyfire/{id}/result")
async def ClassyFire_result(id: str):
    """
    Get the ClassyFire classification results using ID.

    - **ID**: required (query)
    """
    if id:
        data = await result(id)
        return data
