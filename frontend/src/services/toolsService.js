// This file contains functions to interact with the tools API for chemical structure generation, sugar detection, and filtering.
import api from './api';

const TOOLS_URL = '/tools';

/**
 * Generate chemical structures based on a molecular formula
 * @param {string} molecularFormula - Molecular formula (e.g., "C6H6")
 * @returns {Promise<Array<string>>} - Array of generated SMILES strings
 */
export const generateStructures = async (molecularFormula) => {
  try {
    const response = await api.get(`${TOOLS_URL}/generate-structures`, {
      params: { molecular_formula: molecularFormula }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate structures: ${error.message}`);
  }
};

/**
 * Get information about sugar moieties in a molecule
 * @param {string} smiles - SMILES string
 * @returns {Promise<string>} - Sugar information message
 */
export const getSugarInfo = async (smiles) => {
  try {
    const response = await api.get(`${TOOLS_URL}/sugars-info`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to get sugar information: ${error.message}`);
  }
};

/**
 * Detect sugar moieties in a molecule (alias for getSugarInfo)
 * @param {string} smiles - SMILES string
 * @returns {Promise<string>} - Sugar information message
 */
export const detectSugars = getSugarInfo;

/**
 * Remove linear sugars from a molecule
 * @param {string} smiles - SMILES string
 * @returns {Promise<string>} - SMILES with linear sugars removed
 */
export const removeLinearSugars = async (smiles) => {
  try {
    const response = await api.get(`${TOOLS_URL}/remove-linear-sugars`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to remove linear sugars: ${error.message}`);
  }
};

/**
 * Remove circular sugars from a molecule
 * @param {string} smiles - SMILES string
 * @returns {Promise<string>} - SMILES with circular sugars removed
 */
export const removeCircularSugars = async (smiles) => {
  try {
    const response = await api.get(`${TOOLS_URL}/remove-circular-sugars`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to remove circular sugars: ${error.message}`);
  }
};

/**
 * Remove both linear and circular sugars from a molecule
 * @param {string} smiles - SMILES string
 * @returns {Promise<string>} - SMILES with all sugars removed
 */
export const removeAllSugars = async (smiles) => {
  try {
    const response = await api.get(`${TOOLS_URL}/remove-sugars`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to remove all sugars: ${error.message}`);
  }
};

/**
 * General sugar removal function (alias for removeAllSugars)
 * @param {string} smiles - SMILES string
 * @param {string} type - Type of sugars to remove ('all', 'linear', 'circular')
 * @returns {Promise<string>} - SMILES with sugars removed
 */
export const removeSugars = async (smiles, type = 'all') => {
  try {
    switch (type) {
      case 'linear':
        return await removeLinearSugars(smiles);
      case 'circular':
        return await removeCircularSugars(smiles);
      case 'all':
      default:
        return await removeAllSugars(smiles);
    }
  } catch (error) {
    throw new Error(`Failed to remove sugars: ${error.message}`);
  }
};

/**
 * Apply multiple chemical filters to a list of molecules
 * @param {string} smilesList - Newline-separated list of SMILES strings
 * @param {Object} filterOptions - Filter options to apply
 * @returns {Promise<Array<string>>} - Filtered results
 */
export const applyChemicalFilters = async (smilesList, filterOptions = {}) => {
  const {
    pains = true,
    lipinski = true,
    veber = true,
    reos = true,
    ghose = true,
    ruleofthree = true,
    qedscore = '0-10',
    sascore = '0-10',
    nplikeness = '0-10'
  } = filterOptions;

  try {
    const response = await api.post(`/chem/all_filters`, smilesList, {
      params: {
        pains,
        lipinski,
        veber,
        reos,
        ghose,
        ruleofthree,
        qedscore,
        sascore,
        nplikeness
      },
      headers: {
        'Content-Type': 'text/plain'
      }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to apply chemical filters: ${error.message}`);
  }
};

/**
 * Apply chemical filters with detailed violation information
 * @param {string} smilesList - Newline-separated list of SMILES strings
 * @param {Object} filterOptions - Filter options to apply
 * @returns {Promise<Object>} - Detailed filter results with violation information
 */
export const applyChemicalFiltersDetailed = async (smilesList, filterOptions = {}) => {
  const {
    pains = true,
    lipinski = true,
    veber = true,
    reos = true,
    ghose = true,
    ruleofthree = true,
    qedscore = '0-1',
    sascore = '0-10',
    nplikeness = '-5-5',
    filterOperator = 'OR'
  } = filterOptions;

  try {
    const response = await api.post(`/chem/all_filters_detailed`, smilesList, {
      params: {
        pains,
        lipinski,
        veber,
        reos,
        ghose,
        ruleofthree,
        qedscore,
        sascore,
        nplikeness,
        filterOperator
      },
      headers: {
        'Content-Type': 'text/plain',
        'Accept': 'application/json'
      }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to apply detailed chemical filters: ${error.message}`);
  }
};

/**
 * Standardize a molecule using the ChEMBL curation pipeline
 * @param {string} molblock - Molecule in molblock format
 * @returns {Promise<Object>} - Standardized molecule data
 */
export const standardizeMolecule = async (molblock) => {
  try {
    const response = await api.post(`/chem/standardize`, molblock, {
      headers: {
        'Content-Type': 'text/plain'
      }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to standardize molecule: ${error.message}`);
  }
};

/**
 * Get ClassyFire classification for a molecule
 * @param {string} smiles - SMILES string
 * @returns {Promise<Object>} - ClassyFire job data
 */
export const classifyMolecule = async (smiles) => {
  try {
    const response = await api.get(`/chem/classyfire/classify`, {
      params: { smiles }
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to classify molecule: ${error.message}`);
  }
};

/**
 * Get ClassyFire classification results
 * @param {string} jobId - ClassyFire job ID
 * @returns {Promise<Object>} - ClassyFire classification results
 */
export const getClassificationResults = async (jobId) => {
  try {
    const response = await api.get(`/chem/classyfire/${jobId}/result`);
    return response.data;
  } catch (error) {
    throw new Error(`Failed to get classification results: ${error.message}`);
  }
};

// Assemble all functions into a service object
const toolsService = {
  generateStructures,
  getSugarInfo,
  detectSugars, // Add alias
  removeLinearSugars,
  removeCircularSugars,
  removeAllSugars,
  removeSugars, // Add new function
  applyChemicalFilters,
  applyChemicalFiltersDetailed, // Add detailed filters function
  standardizeMolecule,
  classifyMolecule,
  getClassificationResults
};

export default toolsService;