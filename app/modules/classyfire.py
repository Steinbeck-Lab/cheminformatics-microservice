import requests
import json


async def classify(smiles: str) -> dict:
    """
    This function queries the ClassyFire API to classify a chemical compound
    represented by a SMILES string.

    Args:
        smiles (str): A SMILES string representing the chemical compound.

    Returns:
        dict: A dictionary containing the response from the ClassyFire API.

    Raises:
        requests.RequestException: If there's an issue with the API request.
    """

    # Replace any spaces in the SMILES string with '+'
    if any(char.isspace() for char in smiles):
        smiles = smiles.replace(" ", "+")

    # ClassyFire API endpoint URL
    url = "http://classyfire.wishartlab.com/queries/?format=json"

    # Prepare payload for the API request
    payload = json.dumps(
        {"label": "query", "query_input": smiles, "query_type": "STRUCTURE"}
    )

    # Set headers for the API request
    headers = {"Content-Type": "application/json"}

    try:
        # Make a POST request to the API
        response = requests.post(url, headers=headers, data=payload)
        response.raise_for_status()  # Raise exception for HTTP errors
        return response.json()
    except requests.RequestException as e:
        # Handle request-related errors
        raise e


async def result(id: int) -> dict:
    """
    Fetches JSON response from the ClassyFire API for a given ID.

    This function takes an ID and retrieves the corresponding chemical classification
    information from the ClassyFire API in JSON format.

    Args:
        id (int): The ID associated with the chemical compound.

    Returns:
        dict: A dictionary containing ClassyFire classification results.
              The structure of the dictionary includes various classification
              details of the chemical compound, such as class, superclass, direct
              parent, etc.

    Raises:
        requests.exceptions.RequestException: If there is an issue with the HTTP request
            to the ClassyFire API.
    """
    url = f"http://classyfire.wishartlab.com/queries/{id}?format=json"

    headers = {"Content-Type": "application/json"}

    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise an exception for 4xx/5xx status codes
        return response.json()
    except requests.exceptions.RequestException as e:
        # Handle connection errors, timeouts, etc.
        raise e
