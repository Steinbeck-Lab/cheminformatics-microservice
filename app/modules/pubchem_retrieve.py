import re
import logging
from functools import lru_cache
from typing import Optional, Dict, Any, List
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

        # Create a cached version of the internal _get_smiles method
        self._get_smiles_cached = lru_cache(maxsize=self.cache_size)(self._get_smiles)

    def detect_input_type(self, user_input: str) -> str:
        """
        Detect the type of chemical identifier from the input string.

        Args:
            user_input (str): The input string to analyze

        Returns:
            str: The detected input type
        """
        if not user_input or not isinstance(user_input, str):
            return "unknown"

        user_input = user_input.strip()

        # 1. CID: all digits
        if user_input.isdigit():
            return "CID"

        # 2. InChI: starts with "InChI="
        if user_input.startswith("InChI="):
            return "InChI"

        # 3. InChIKey: specific format
        if re.match(r"^[A-Z0-9]{14}-[A-Z0-9]{10}-[A-Z0-9]$", user_input):
            return "InChIKey"

        # 4. CAS number: pattern like "50-78-2" or "7732-18-5"
        if re.match(r"^\d{2,7}-\d{2}-\d$", user_input):
            return "CAS"

        # 5. SMILES: check before molecular formula since some short formulas look like SMILES
        # SMILES typically contain organic chemistry characters and structural notation
        smiles_pattern = r"^[A-Za-z0-9\(\)\[\]=#+\-\\/\\@\.%:]*$"

        if (
            " " not in user_input  # SMILES shouldn't contain spaces
            and len(user_input) >= 1  # Allow very short SMILES
            and len(user_input) <= 500  # Reasonable maximum length
            and re.match(smiles_pattern, user_input)
            # Improved heuristic: SMILES usually contain structural elements OR basic organic patterns
            and (
                any(char in user_input for char in "()=[]#@\\/")  # Structural notation
                or re.match(
                    r"^[CNOPS]+(Cl|Br|[cnops]|\d)*$", user_input
                )  # Simple organic patterns
                or (
                    len(user_input) <= 10
                    and re.match(r"^[CNOSPcnops]+\d*$", user_input)
                )  # Short organic
                or (
                    len(user_input) >= 3
                    and re.match(r"^[CNOSHPFcnoshpf]+$", user_input)
                )  # Common atoms
            )
        ):
            return "SMILES"

        # 6. Molecular formula: pattern like "C9H8O4" (check after SMILES)
        if re.match(r"^[A-Z][a-z]?(\d*[A-Z][a-z]?\d*)*$", user_input):
            return "formula"

        # 7. Default: chemical name
        return "name"

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

        url = f"{self.BASE_URL}/cid/{cid}/property/SMILES/JSON"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data["PropertyTable"]["Properties"][0]["SMILES"]
        except (RequestException, KeyError, IndexError) as e:
            logger.error(f"Error querying by CID {cid}: {str(e)}")
            return None

    def _get_cids_by_identifier(
        self, identifier: str, input_type: str
    ) -> Optional[List[str]]:
        """
        Get CIDs for a given identifier based on its type.

        Args:
            identifier (str): The chemical identifier
            input_type (str): The type of identifier

        Returns:
            Optional[List[str]]: List of CIDs if found, None otherwise
        """
        try:
            if input_type == "CID":
                return [identifier] if identifier.isdigit() else None

            elif input_type == "InChI":
                url = f"{self.BASE_URL}/inchi/cids/txt"
                response = self.session.post(
                    url, data={"inchi": identifier}, timeout=self.timeout
                )
                response.raise_for_status()
                cids = response.text.strip().splitlines()
                return cids[:10]  # Limit to first 10 CIDs

            elif input_type == "InChIKey":
                url = f"{self.BASE_URL}/inchikey/{quote(identifier, safe='')}/cids/txt"
                response = self.session.get(url, timeout=self.timeout)
                response.raise_for_status()
                cids = response.text.strip().splitlines()
                return cids[:10]  # Limit to first 10 CIDs

            elif input_type == "formula":
                url = (
                    f"{self.BASE_URL}/fastformula/{quote(identifier, safe='')}/cids/txt"
                )
                response = self.session.get(url, timeout=self.timeout)
                response.raise_for_status()
                cids = response.text.strip().splitlines()
                return cids[:10]  # Limit to first 10 CIDs

            elif input_type == "SMILES":
                url = f"{self.BASE_URL}/smiles/cids/txt"
                response = self.session.post(
                    url, data={"smiles": identifier}, timeout=self.timeout
                )
                response.raise_for_status()
                cids = response.text.strip().splitlines()
                return cids[:10]  # Limit to first 10 CIDs

            elif input_type in ["name", "CAS"]:
                encoded = quote(identifier, safe="")
                url = f"{self.BASE_URL}/name/{encoded}/cids/txt"
                response = self.session.get(url, timeout=self.timeout)
                response.raise_for_status()
                cids = response.text.strip().splitlines()
                return cids[:10]  # Limit to first 10 CIDs

        except (RequestException, IndexError) as e:
            logger.error(f"Error getting CIDs for {identifier}: {str(e)}")
            return None

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
        url = (
            f"{self.BASE_URL}/inchikey/{quote(inchikey, safe='')}/property/SMILES/JSON"
        )

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data["PropertyTable"]["Properties"][0]["SMILES"]
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
        url = f"{self.BASE_URL}/name/{encoded}/property/SMILES/JSON"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            return data["PropertyTable"]["Properties"][0]["SMILES"]
        except (RequestException, KeyError, IndexError) as e:
            logger.error(f"Error querying by name: {str(e)}")
            return None

    def _get_smiles(self, user_input: str) -> Optional[str]:
        """
        Internal implementation of get_smiles without caching.
        This will be wrapped with LRU cache in the constructor.
        """
        if not user_input or not isinstance(user_input, str):
            logger.error(f"Invalid input: {user_input}")
            return None

        # Trim whitespace
        user_input = user_input.strip()
        input_type = self.detect_input_type(user_input)

        # Route to appropriate query method based on detected type
        if input_type == "CID":
            return self._query_by_cid(user_input)
        elif input_type == "InChI":
            return self._query_by_inchi(user_input)
        elif input_type == "InChIKey":
            return self._query_by_inchikey(user_input)
        elif input_type == "formula":
            return self._query_by_formula(user_input)
        elif input_type == "SMILES":
            return self._query_by_smiles(user_input)
        else:  # name or CAS
            return self._query_by_name(user_input)

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
        return self._get_smiles_cached(user_input)

    def get_compound_info(self, user_input: str) -> Dict[str, Any]:
        """
        Retrieve comprehensive compound information from PubChem.

        Args:
            user_input (str): The chemical identifier

        Returns:
            Dict[str, Any]: Dictionary containing:
                - input: original input
                - input_type: detected input type
                - canonical_smiles: canonical SMILES if found
                - cids: list of PubChem CIDs
                - pubchem_links: list of PubChem compound page URLs
                - success: boolean indicating if compound was found
        """
        if not user_input or not isinstance(user_input, str):
            return {
                "input": user_input,
                "input_type": "unknown",
                "canonical_smiles": None,
                "cids": None,
                "pubchem_links": None,
                "success": False,
            }

        user_input = user_input.strip()
        input_type = self.detect_input_type(user_input)

        # Get CIDs
        cids = self._get_cids_by_identifier(user_input, input_type)

        if not cids:
            return {
                "input": user_input,
                "input_type": input_type,
                "canonical_smiles": None,
                "cids": None,
                "pubchem_links": None,
                "success": False,
            }

        # Get canonical SMILES from first CID
        canonical_smiles = self._query_by_cid(cids[0])

        # Generate PubChem links
        pubchem_links = [
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" for cid in cids
        ]

        return {
            "input": user_input,
            "input_type": input_type,
            "canonical_smiles": canonical_smiles,
            "cids": cids,
            "pubchem_links": pubchem_links,
            "success": canonical_smiles is not None,
        }
