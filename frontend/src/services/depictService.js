// This module provides functions to generate 2D and 3D depictions of molecules
import api from './api';

const DEPICT_URL = '/depict';

/**
 * Generate a 2D depiction of a molecule
 * @param {string} smiles - SMILES string
 * @param {Object} options - Depiction options
 * @param {string} options.toolkit - Toolkit to use (cdk, rdkit)
 * @param {number} options.width - Width of the generated image
 * @param {number} options.height - Height of the generated image
 * @param {number} options.rotate - Rotation angle in degrees
 * @param {boolean} options.CIP - Include Cahn-Ingold-Prelog stereochemistry
 * @param {boolean} options.unicolor - Use a single color for the molecule
 * @param {string} options.highlight - SMARTS pattern to highlight atoms/bonds
 * @returns {Promise<string>} - SVG depiction as text
 */
export const generate2DDepiction = async (smiles, options = {}) => {
  const {
    toolkit = 'rdkit',
    width = 512,
    height = 512,
    rotate = 0,
    CIP = false,
    unicolor = false,
    highlight = 'COSN'
  } = options;

  try {
    const response = await api.get(`${DEPICT_URL}/2D`, {
      params: {
        smiles,
        toolkit,
        width,
        height,
        rotate,
        CIP,
        unicolor,
        highlight
      },
      responseType: 'text'
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate 2D depiction: ${error.message}`);
  }
};

/**
 * Generate a 3D depiction of a molecule
 * @param {string} smiles - SMILES string
 * @param {string} toolkit - Toolkit to use (rdkit, openbabel)
 * @returns {Promise<string>} - HTML with embedded 3D viewer
 */
export const generate3DDepiction = async (smiles, toolkit = 'openbabel') => {
  try {
    const response = await api.get(`${DEPICT_URL}/3D`, {
      params: { smiles, toolkit },
      responseType: 'text'
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate 3D depiction: ${error.message}`);
  }
};

/**
 * Get the URL for a 2D depiction image
 * @param {string} smiles - SMILES string
 * @param {Object} options - Depiction options
 * @returns {string} - URL to the depiction image
 */
export const get2DDepictionUrl = (smiles, options = {}) => {
  const {
    width = 512,
    height = 512,
    toolkit = 'rdkit',
    rotate = 0,
    CIP = false,
    unicolor = false,
    highlight = '',
    format = 'svg'
  } = options;

  const baseUrl = api.defaults.baseURL || '';
  let url;

  // Support format parameter in the URL
  if (format && format !== 'svg') {
    url = new URL(`${baseUrl}${DEPICT_URL}/2D/${format}`);
  } else {
    url = new URL(`${baseUrl}${DEPICT_URL}/2D`);
  }

  url.searchParams.append('smiles', smiles);
  url.searchParams.append('width', width);
  url.searchParams.append('height', height);
  url.searchParams.append('toolkit', toolkit);
  url.searchParams.append('rotate', rotate);
  url.searchParams.append('CIP', CIP);
  url.searchParams.append('unicolor', unicolor);

  if (highlight) {
    url.searchParams.append('highlight', highlight);
  }

  return url.toString();
};

/**
 * Generate batch 2D depictions for multiple SMILES strings
 * @param {Array<string>} smilesList - Array of SMILES strings
 * @param {Object} options - Depiction options
 * @returns {Promise<Array<Object>>} - Array of depiction results with SMILES and SVG content
 */
export const generateBatchDepictions = async (smilesList, options = {}) => {
  try {
    const results = [];

    for (const smiles of smilesList) {
      const svg = await generate2DDepiction(smiles, options);
      results.push({
        smiles,
        svg
      });
    }

    return results;
  } catch (error) {
    throw new Error(`Failed to generate batch depictions: ${error.message}`);
  }
};

/**
 * Download multiple depictions as a ZIP file
 * @param {Array<Object>} depictions - Array of depiction objects with smiles, title, and svg properties
 * @param {string} format - File format to download (svg, png)
 * @param {Object} options - Depiction options
 * @returns {Promise<Blob>} - ZIP file as a Blob
 */
export const downloadDepictionsAsZip = async (depictions, format = 'svg', options = {}) => {
  // This function requires JSZip to be installed
  // npm install jszip
  try {
    const JSZip = (await import('jszip')).default;
    const zip = new JSZip();

    for (const depiction of depictions) {
      const { smiles, title, svg } = depiction;
      const filename = `${title.replace(/[^a-z0-9]/gi, '_').toLowerCase()}.${format}`;

      if (format === 'svg') {
        zip.file(filename, svg);
      } else {
        // For other formats, we need to fetch from the API
        const url = get2DDepictionUrl(smiles, { ...options, format });
        const response = await fetch(url);
        const blob = await response.blob();
        zip.file(filename, blob);
      }
    }

    return await zip.generateAsync({ type: 'blob' });
  } catch (error) {
    throw new Error(`Failed to create ZIP file: ${error.message}`);
  }
};

/**
 * Highlight a substructure in a molecule
 * @param {string} smiles - SMILES string of the full molecule
 * @param {string} substructure - SMARTS pattern of the substructure to highlight
 * @param {Object} options - Additional depiction options
 * @returns {Promise<string>} - SVG depiction with highlighted substructure
 */
export const highlightSubstructure = async (smiles, substructure, options = {}) => {
  const depictionOptions = {
    ...options,
    highlight: substructure
  };

  return generate2DDepiction(smiles, depictionOptions);
};

/**
 * Generate a color-coded depiction based on atom properties
 * @param {string} smiles - SMILES string
 * @param {string} colorBy - Property to color by (e.g., 'charge', 'element')
 * @param {Object} options - Additional depiction options
 * @returns {Promise<string>} - SVG depiction with color-coded atoms
 */
export const generateColorCodedDepiction = async (smiles, colorBy = 'element', options = {}) => {
  // This is a placeholder - the actual API might not support this directly
  // but we can implement it as a convenience method that uses other endpoints
  const depictionOptions = {
    ...options,
    unicolor: colorBy === 'none'
  };

  return generate2DDepiction(smiles, depictionOptions);
};

/**
 * Generate multiple depictions with different settings
 * @param {string} smiles - SMILES string
 * @param {Array<Object>} optionsArray - Array of depiction option objects
 * @returns {Promise<Array<string>>} - Array of SVG depictions
 */
export const generateMultipleDepictions = async (smiles, optionsArray = [{}]) => {
  try {
    const promises = optionsArray.map(options =>
      generate2DDepiction(smiles, options)
    );

    return Promise.all(promises);
  } catch (error) {
    throw new Error(`Failed to generate multiple depictions: ${error.message}`);
  }
};

/**
 * Parse SMILES strings from a text input with optional titles
 * @param {string} text - Input text with one SMILES per line, optionally with titles
 * @returns {Array<Object>} - Array of objects with smiles and title properties
 */
export const parseSmilesInput = (text) => {
  return text
    .split('\n')
    .map(line => {
      line = line.trim();
      // Skip empty lines and comments
      if (line.length === 0 || line.startsWith('#')) {
        return null;
      }

      // Extract title if present (after space or tab)
      const spaceIndex = line.search(/[\s\t]/);
      const smiles = spaceIndex > 0 ? line.substring(0, spaceIndex) : line;
      const title = spaceIndex > 0 ? line.substring(spaceIndex + 1).trim() : `Molecule`;

      return { smiles, title };
    })
    .filter(item => item !== null);
};

// Assemble all functions into a service object
const depictService = {
  generate2DDepiction,
  generate3DDepiction,
  get2DDepictionUrl,
  highlightSubstructure,
  generateColorCodedDepiction,
  generateMultipleDepictions,
  generateBatchDepictions,
  downloadDepictionsAsZip,
  parseSmilesInput
};

export default depictService;