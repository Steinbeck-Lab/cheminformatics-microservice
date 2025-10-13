// This module provides functions to convert chemical representations using a REST API.
import api from "./api";

const CONVERT_URL = "/convert";

/**
 * Generate 2D coordinates for a molecule
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (cdk, rdkit, openbabel)
 * @returns {Promise<string>} - Mol block with 2D coordinates
 */
export const generate2DCoordinates = async (smiles, toolkit = "cdk") => {
  try {
    const response = await api.get(`${CONVERT_URL}/mol2D`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate 2D coordinates: ${error.message}`);
  }
};

/**
 * Generate 3D coordinates for a molecule
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (rdkit, openbabel)
 * @returns {Promise<string>} - Mol block with 3D coordinates
 */
export const generate3DCoordinates = async (smiles, toolkit = "openbabel") => {
  try {
    const response = await api.get(`${CONVERT_URL}/mol3D`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate 3D coordinates: ${error.message}`);
  }
};

/**
 * Generate SMILES from IUPAC name or SELFIES
 * @param {string} inputText - IUPAC name or SELFIES
 * @param {string} representation - Input representation type (iupac, selfies)
 * @param {string} converter - Converter to use for IUPAC (opsin)
 * @returns {Promise<string>} - Generated SMILES
 */
export const generateSMILES = async (inputText, representation = "iupac", converter = "opsin") => {
  try {
    const response = await api.get(`${CONVERT_URL}/smiles`, {
      params: { input_text: inputText, representation, converter },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate SMILES: ${error.message}`);
  }
};

/**
 * Generate canonical SMILES
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (cdk, rdkit, openbabel)
 * @returns {Promise<string>} - Canonical SMILES
 */
export const generateCanonicalSMILES = async (smiles, toolkit = "cdk") => {
  try {
    const response = await api.get(`${CONVERT_URL}/canonicalsmiles`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate canonical SMILES: ${error.message}`);
  }
};

/**
 * Generate CXSMILES (ChemAxon Extended SMILES)
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (cdk, rdkit)
 * @returns {Promise<string>} - CXSMILES
 */
export const generateCXSMILES = async (smiles, toolkit = "cdk") => {
  try {
    const response = await api.get(`${CONVERT_URL}/cxsmiles`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate CXSMILES: ${error.message}`);
  }
};

/**
 * Generate InChI (IUPAC International Chemical Identifier)
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (cdk, rdkit, openbabel)
 * @returns {Promise<string>} - InChI
 */
export const generateInChI = async (smiles, toolkit = "cdk") => {
  try {
    const response = await api.get(`${CONVERT_URL}/inchi`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate InChI: ${error.message}`);
  }
};

/**
 * Generate InChI Key
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (cdk, rdkit, openbabel)
 * @returns {Promise<string>} - InChI Key
 */
export const generateInChIKey = async (smiles, toolkit = "cdk") => {
  try {
    const response = await api.get(`${CONVERT_URL}/inchikey`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate InChI Key: ${error.message}`);
  }
};

/**
 * Generate SELFIES (Self-Referencing Embedded Strings)
 * @param {string} smiles - SMILES string
 * @returns {Promise<string>} - SELFIES
 */
export const generateSELFIES = async (smiles) => {
  try {
    const response = await api.get(`${CONVERT_URL}/selfies`, {
      params: { smiles },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate SELFIES: ${error.message}`);
  }
};

/**
 * Generate SMARTS pattern
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (rdkit)
 * @returns {Promise<string>} - SMARTS pattern
 */
export const generateSMARTS = async (smiles, toolkit = "rdkit") => {
  try {
    const response = await api.get(`${CONVERT_URL}/smarts`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate SMARTS: ${error.message}`);
  }
};

/**
 * Convert MOL/SDF block to SMILES
 * @param {string} molblock - MOL or SDF block string
 * @param {string} toolkit - Toolkit to use (cdk, rdkit)
 * @returns {Promise<string>} - SMILES string
 */
export const molblockToSMILES = async (molblock, toolkit = "cdk") => {
  try {
    // Get the base URL from the api instance or use default
    const baseURL = api.defaults.baseURL || "https://dev.api.naturalproducts.net/latest";

    // Use fetch API instead of axios to ensure proper text/plain handling
    const response = await fetch(`${baseURL}${CONVERT_URL}/molblock?toolkit=${toolkit}`, {
      method: "POST",
      headers: {
        "Content-Type": "text/plain; charset=utf-8",
        Accept: "application/json",
      },
      body: molblock, // Send molblock as plain text
    });

    if (!response.ok) {
      let errorMsg = `Error ${response.status}: ${response.statusText}`;
      try {
        const errorData = await response.json();
        errorMsg = errorData.detail || errorMsg;
      } catch (jsonError) {
        // Ignore if response is not JSON
      }
      throw new Error(errorMsg);
    }

    // Get the response as text and clean it up
    let result = await response.text();

    // Remove surrounding quotes if present
    if (result.startsWith('"') && result.endsWith('"')) {
      result = result.substring(1, result.length - 1);
    }

    return result;
  } catch (error) {
    throw new Error(`Failed to convert MOL/SDF to SMILES: ${error.message}`);
  }
};

/**
 * Generate multiple formats at once
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (cdk, rdkit, openbabel)
 * @returns {Promise<Object>} - Object with multiple formats
 */
export const generateMultipleFormats = async (smiles, toolkit = "cdk") => {
  try {
    const response = await api.get(`${CONVERT_URL}/formats`, {
      params: { smiles, toolkit },
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate multiple formats: ${error.message}`);
  }
};

// Assemble all functions into a service object
const convertService = {
  generate2DCoordinates,
  generate3DCoordinates,
  generateSMILES,
  generateCanonicalSMILES,
  generateCXSMILES,
  generateInChI,
  generateInChIKey,
  generateSELFIES,
  generateSMARTS,
  molblockToSMILES,
  generateMultipleFormats,
};

export default convertService;
