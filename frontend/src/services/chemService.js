// This file contains functions to interact with the chemical structure API
import api from './api';

const CHEM_URL = '/chem';

/**
 * Generate all possible stereoisomers for a given SMILES string
 * @param {string} smiles - SMILES string
 * @returns {Promise<Array<string>>} - Array of stereoisomer SMILES strings
 */
const generateStereoisomers = async (smiles) => {
  try {
    const response = await api.get(`${CHEM_URL}/stereoisomers`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate stereoisomers: ${error.message}`);
  }
};

/**
 * Calculate molecular descriptors for a given SMILES string
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (rdkit, cdk, all)
 * @param {string} format - Output format (json, html)
 * @returns {Promise<Object>} - Object containing calculated descriptors
 */
const calculateDescriptors = async (smiles, toolkit = 'rdkit', format = 'json') => {
  try {
    const response = await api.get(`${CHEM_URL}/descriptors`, {
      params: { smiles, toolkit, format }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to calculate descriptors: ${error.message}`);
  }
};

/**
 * Calculate descriptors for multiple SMILES strings
 * @param {string} smilesList - Comma-separated list of SMILES strings
 * @param {string} toolkit - Toolkit to use (rdkit, cdk)
 * @returns {Promise<Object>} - Object with SMILES as keys and descriptors as values
 */
const calculateMultipleDescriptors = async (smilesList, toolkit = 'rdkit') => {
  try {
    const response = await api.get(`${CHEM_URL}/descriptors/multiple`, {
      params: { smiles: smilesList, toolkit }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to calculate multiple descriptors: ${error.message}`);
  }
};

/**
 * Generate HOSE codes for a given SMILES string
 * @param {string} smiles - SMILES string
 * @param {number} spheres - Number of spheres for HOSE code generation
 * @param {string} toolkit - Toolkit to use (cdk, rdkit)
 * @param {boolean} ringsize - Whether to include ringsize information
 * @returns {Promise<Array<string>>} - Array of HOSE codes
 */
const generateHOSECodes = async (smiles, spheres = 2, toolkit = 'rdkit', ringsize = false) => {
  try {
    const response = await api.get(`${CHEM_URL}/HOSEcode`, {
      params: { smiles, spheres, toolkit, ringsize }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate HOSE codes: ${error.message}`);
  }
};

/**
 * Check a SMILES string for errors and optionally standardize it
 * @param {string} smiles - SMILES string
 * @param {boolean} fix - Whether to fix issues by standardizing
 * @returns {Promise<Object>} - Object containing validation results
 */
const checkStructureErrors = async (smiles, fix = false) => {
  try {
    const response = await api.get(`${CHEM_URL}/errors`, {
      params: { smiles, fix }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to check structure errors: ${error.message}`);
  }
};

/**
 * Calculate natural product (NP) likeness score
 * @param {string} smiles - SMILES string
 * @returns {Promise<number>} - NP likeness score
 */
const calculateNPLikeness = async (smiles) => {
  try {
    const response = await api.get(`${CHEM_URL}/nplikeness/score`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to calculate NP likeness score: ${error.message}`);
  }
};

/**
 * Calculate Tanimoto similarity between two or more molecules
 * @param {string} smiles - Comma-separated list of SMILES strings
 * @param {string} toolkit - Toolkit to use (cdk, rdkit)
 * @param {string} fingerprinter - Fingerprint type (RDKit, Atompairs, MACCS, PubChem, ECFP, MAPC)
 * @param {number} nBits - Number of bits for fingerprint vectors
 * @param {number} radius - Radius for ECFP fingerprints
 * @returns {Promise<number|string>} - Similarity score or HTML table for multiple molecules
 */
const calculateTanimotoSimilarity = async (
  smiles,
  toolkit = 'rdkit',
  fingerprinter = 'ECFP',
  nBits = 2048,
  radius = 2
) => {
  try {
    const response = await api.get(`${CHEM_URL}/tanimoto`, {
      params: { smiles, toolkit, fingerprinter, nBits, radius }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to calculate Tanimoto similarity: ${error.message}`);
  }
};

/**
 * Preprocess a molecule for COCONUT database submissions
 * @param {string} smiles - SMILES string
 * @param {boolean} _3d_mol - Whether to generate 3D coordinates
 * @param {boolean} descriptors - Whether to calculate descriptors
 * @returns {Promise<Object>} - Object containing preprocessed data
 */
const coconutPreprocessing = async (smiles, _3d_mol = false, descriptors = false) => {
  try {
    const response = await api.get(`${CHEM_URL}/coconut/pre-processing`, {
      params: { smiles, _3d_mol, descriptors }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to preprocess for COCONUT: ${error.message}`);
  }
};

/**
 * Generate functional groups for a molecule using Ertl's algorithm
 * @param {string} smiles - SMILES string
 * @returns {Promise<Array>} - Array of functional groups
 */
const generateFunctionalGroups = async (smiles) => {
  try {
    const response = await api.get(`${CHEM_URL}/ertlfunctionalgroup`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate functional groups: ${error.message}`);
  }
};

/**
 * Generate standardized tautomer for a molecule
 * @param {string} smiles - SMILES string
 * @returns {Promise<string>} - Standardized tautomer SMILES
 */
const generateStandardizedTautomer = async (smiles) => {
  try {
    const response = await api.get(`${CHEM_URL}/standarizedTautomer`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate standardized tautomer: ${error.message}`);
  }
};


/**
 * Find chemical structure by name, formula, or identifier using the PubChem database
 * @param {string} identifier - Chemical name, formula, or identifier (CAS, InChI, SMILES, etc.)
 * @returns {Promise<Object>} - Object containing canonical_smiles, input, input_type, and success status
 */
const lookupPubChem = async (identifier) => {
  try {
    const response = await api.get(`${CHEM_URL}/pubchem/smiles`, {
      params: { identifier }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to look up structure in PubChem: ${error.message}`);
  }
};

// Create an alias for calculateDescriptors as getDescriptors for backward compatibility
const getDescriptors = calculateDescriptors;


// Export individual functions
export {
  generateStereoisomers,
  calculateDescriptors,
  calculateMultipleDescriptors,
  generateHOSECodes,
  checkStructureErrors,
  calculateNPLikeness,
  calculateTanimotoSimilarity,
  coconutPreprocessing,
  generateFunctionalGroups,
  generateStandardizedTautomer,
  getDescriptors,
  lookupPubChem
};

// Assemble all functions into a service object
const chemService = {
  generateStereoisomers,
  calculateDescriptors,
  getDescriptors,  // Include the alias
  calculateMultipleDescriptors,
  generateHOSECodes,
  checkStructureErrors,
  calculateNPLikeness,
  calculateTanimotoSimilarity,
  coconutPreprocessing,
  generateFunctionalGroups,
  generateStandardizedTautomer,
  lookupPubChem
};

export default chemService;