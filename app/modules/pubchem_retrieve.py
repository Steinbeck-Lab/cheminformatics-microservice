import re
import logging
from functools import lru_cache
from typing import Optional
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from requests.exceptions import RequestException
from urllib3.util.retry import Retry


# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("pubchem_client")


class PubChemClient:
    """Client for interacting with the PubChem PUG REST API to retrieve chemical information."""

    # API Constants
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    DEFAULT_TIMEOUT = 10  # seconds
    MAX_RETRIES = 3
    BACKOFF_FACTOR = 0.5

    def __init__(
        self, timeout: int = None, max_retries: int = None, cache_size: int = 128
    ):
        """
        Initialize the PubChem client.

        Args:
            timeout (int, optional): Request timeout in seconds. Defaults to DEFAULT_TIMEOUT.
            max_retries (int, optional): Maximum number of retries for failed requests. Defaults to MAX_RETRIES.
            cache_size (int, optional): Size of the LRU cache. Defaults to 128.
        """
        self.timeout = timeout or self.DEFAULT_TIMEOUT
        self.max_retries = max_retries or self.MAX_RETRIES
        self.cache_size = cache_size

        # Set up a session with retry capabilities
        self.session = requests.Session()
        retry_strategy = Retry(
            total=self.max_retries,
            backoff_factor=self.BACKOFF_FACTOR,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    @lru_cache(maxsize=128)
    def _query_by_cid(self, cid: str) -> Optional[str]:
        """
        Retrieve the canonical SMILES for a compound by its PubChem CID.

        Args:
            cid (str): The PubChem Compound ID.

        Returns:
            Optional[str]: The canonical SMILES string if found, None otherwise.
        """
        if not cid.isdigit():
            logger.error(f"Invalid CID format: {cid}")
            return None

        url = f"{self.BASE_URL}/cid/{cid}/property/IsomericSMILES/JSON"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
        except (RequestException, KeyError, IndexError) as e:
            logger.error(f"Error querying by CID {cid}: {str(e)}")
            return None

    def _query_by_inchi(self, inchi: str) -> Optional[str]:
        """
        Retrieve the canonical SMILES for a compound by its InChI string.

        Args:
            inchi (str): The InChI string.

        Returns:
            Optional[str]: The canonical SMILES string if found, None otherwise.
        """
        url = f"{self.BASE_URL}/inchi/cids/txt"

        try:
            response = self.session.post(
                url, data={"inchi": inchi}, timeout=self.timeout
            )
            response.raise_for_status()
            cid = response.text.strip().splitlines()[0]
            return self._query_by_cid(cid)
        except (RequestException, IndexError) as e:
            logger.error(f"Error querying by InChI: {str(e)}")
            return None

    def _query_by_inchikey(self, inchikey: str) -> Optional[str]:
        """
        Retrieve the canonical SMILES for a compound by its InChIKey.

        Args:
            inchikey (str): The InChIKey.

        Returns:
            Optional[str]: The canonical SMILES string if found, None otherwise.
        """
        url = f"{self.BASE_URL}/inchikey/{quote(inchikey, safe='')}/property/IsomericSMILES/JSON"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
        except (RequestException, KeyError, IndexError) as e:
            logger.error(f"Error querying by InChIKey: {str(e)}")
            return None

    def _query_by_formula(self, formula: str) -> Optional[str]:
        """
        Retrieve the canonical SMILES for a compound by its molecular formula.

        Args:
            formula (str): The molecular formula.

        Returns:
            Optional[str]: The canonical SMILES string if found, None otherwise.
        """
        url = f"{self.BASE_URL}/fastformula/{quote(formula, safe='')}/cids/txt"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            cid = response.text.strip().splitlines()[0]
            return self._query_by_cid(cid)
        except (RequestException, IndexError) as e:
            logger.error(f"Error querying by formula: {str(e)}")
            return None

    def _query_by_smiles(self, smiles: str) -> Optional[str]:
        """
        Retrieve the canonical SMILES for a compound by its SMILES string.

        Args:
            smiles (str): The SMILES string.

        Returns:
            Optional[str]: The canonical SMILES string if found, None otherwise.
        """
        url = f"{self.BASE_URL}/smiles/cids/txt"

        try:
            response = self.session.post(
                url, data={"smiles": smiles}, timeout=self.timeout
            )
            response.raise_for_status()
            cid = response.text.strip().splitlines()[0]
            return self._query_by_cid(cid)
        except (RequestException, IndexError) as e:
            logger.error(f"Error querying by SMILES: {str(e)}")
            return None

    def _query_by_name(self, name: str) -> Optional[str]:
        """
        Retrieve the canonical SMILES for a compound by its name.

        Args:
            name (str): The chemical name (IUPAC, synonym, trivial name, etc.) or CAS number.

        Returns:
            Optional[str]: The canonical SMILES string if found, None otherwise.
        """
        encoded = quote(name, safe="")
        url = f"{self.BASE_URL}/name/{encoded}/property/IsomericSMILES/JSON"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
        except (RequestException, KeyError, IndexError) as e:
            logger.error(f"Error querying by name: {str(e)}")
            return None

    def cache_with_dynamic_size(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            cache = lru_cache(maxsize=self.cache_size)(func)
            return cache(self, *args, **kwargs)
        return wrapper

    @cache_with_dynamic_size
    def get_smiles(self, user_input: str) -> Optional[str]:
        """
        Retrieve the canonical SMILES for a molecule from PubChem via the PUG REST API.

        The function supports multiple input types:
          - CID: a string of digits.
          - InChI: a string starting with "InChI=".
          - InChIKey: matching the standard format (e.g., "LFQSCWFLJHTTHZ-UHFFFAOYSA-N").
          - CAS number: if the input matches a pattern like "7732-18-5" (handled as a name).
          - Molecular formula: if the input matches a formula pattern (e.g., "C6H12O6").
          - SMILES: if the input contains no spaces, is short (≤10 chars), and is composed solely of characters typically found in SMILES.
          - Chemical name: default option (covers IUPAC names, synonyms, trivial names, etc.).

        Args:
            user_input (str): The chemical identifier.

        Returns:
            Optional[str]: The canonical SMILES string if found, None otherwise.

        Examples:
            >>> client = PubChemClient()
            >>> client.get_smiles("2244")  # Aspirin by CID
            'CC(=O)OC1=CC=CC=C1C(=O)O'
            >>> client.get_smiles("aspirin")  # Aspirin by name
            'CC(=O)OC1=CC=CC=C1C(=O)O'
            >>> client.get_smiles("C6H12O6")  # Glucose by formula
            'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'
        """
        if not user_input or not isinstance(user_input, str):
            logger.error(f"Invalid input: {user_input}")
            return None

        # Trim whitespace
        user_input = user_input.strip()

        # 1. If input is a CID (all digits)
        if user_input.isdigit():
            return self._query_by_cid(user_input)

        # 2. If input is an InChI (starts with "InChI=")
        if user_input.startswith("InChI="):
            return self._query_by_inchi(user_input)

        # 3. If input is an InChIKey
        if re.match(r"^[A-Z0-9]{14}-[A-Z0-9]{10}-[A-Z0-9]$", user_input):
            return self._query_by_inchikey(user_input)

        # 4. If input is a CAS number (e.g., "7732-18-5")
        if re.match(r"^\d{2,7}-\d{2}-\d$", user_input):
            return self._query_by_name(user_input)  # Handle CAS as name

        # 5. If input is a molecular formula (e.g., "C6H12O6")
        if re.match(r"^(?:[A-Z][a-z]?\d+)+$", user_input):
            return self._query_by_formula(user_input)

        # 6. If input is a SMILES string
        if (
            " " not in user_input
            and len(user_input) <= 10
            and re.match(r"^[CHONSPFIClBr0-9@+\-\\/]+$", user_input)
        ):
            return self._query_by_smiles(user_input)

        # 7. Default: treat as a chemical name
        return self._query_by_name(user_input)
