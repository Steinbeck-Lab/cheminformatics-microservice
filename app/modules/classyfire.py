import requests
import json


async def classify(smiles: str):
    """
    This fucntion takes a smiles string and returns a json response
    from classyfire API.
    Args (str): SMILES string.
    Returns (json): classyfire results.
    """
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")
    url = "http://classyfire.wishartlab.com/queries/?format=json"

    payload = json.dumps(
        {"label": "curl_test", "query_input": smiles, "query_type": "STRUCTURE"}
    )
    headers = {"Content-Type": "application/json"}
    response = requests.request("POST", url, headers=headers, data=payload)
    return response.json()


async def result(id):
    """
    This fucntion takes a ID and returns a json response
    from classyfire API.
    Args : ID.
    Returns (json): classyfire results.
    """
    url = "http://classyfire.wishartlab.com/queries/" + str(id) + "?format=json"

    headers = {"Content-Type": "application/json"}
    response = requests.request("GET", url, headers=headers, data={})
    return response.json()
