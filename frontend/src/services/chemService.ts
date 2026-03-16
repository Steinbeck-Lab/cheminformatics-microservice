// This file contains functions to interact with the chemical structure API
import api from "./api";
import type {
  StructureErrorResult,
  DescriptorResult,
  MultipleDescriptorResult,
  FixRadicalsResult,
  StandardizeResult,
  PubChemLookupResult,
  CoconutPreprocessingResult,
} from "../types/api";

const CHEM_URL = "/chem";

/**
 * Generate all possible stereoisomers for a given SMILES string
 */
const generateStereoisomers = async (smiles: string): Promise<string[]> => {
  try {
    const response = await api.get<string[]>(`${CHEM_URL}/stereoisomers`, {
      params: { smiles },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate stereoisomers: ${message}`);
  }
};

/**
 * Calculate molecular descriptors for a given SMILES string
 */
const calculateDescriptors = async (
  smiles: string,
  toolkit = "rdkit",
  format = "json"
): Promise<DescriptorResult> => {
  try {
    const response = await api.get<DescriptorResult>(`${CHEM_URL}/descriptors`, {
      params: { smiles, toolkit, format },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to calculate descriptors: ${message}`);
  }
};

/**
 * Calculate descriptors for multiple SMILES strings
 */
const calculateMultipleDescriptors = async (
  smilesList: string,
  toolkit = "rdkit"
): Promise<MultipleDescriptorResult> => {
  try {
    const response = await api.get<MultipleDescriptorResult>(`${CHEM_URL}/descriptors/multiple`, {
      params: { smiles: smilesList, toolkit },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to calculate multiple descriptors: ${message}`);
  }
};

/**
 * Generate HOSE codes for a given SMILES string
 */
const generateHOSECodes = async (
  smiles: string,
  spheres = 2,
  toolkit = "rdkit",
  ringsize = false
): Promise<string[]> => {
  try {
    const response = await api.get<string[]>(`${CHEM_URL}/HOSEcode`, {
      params: { smiles, spheres, toolkit, ringsize },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate HOSE codes: ${message}`);
  }
};

/**
 * Check a SMILES string for errors and optionally standardize it
 */
const checkStructureErrors = async (smiles: string, fix = false): Promise<StructureErrorResult> => {
  try {
    const response = await api.get<StructureErrorResult>(`${CHEM_URL}/errors`, {
      params: { smiles, fix },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to check structure errors: ${message}`);
  }
};

/**
 * Calculate natural product (NP) likeness score
 */
const calculateNPLikeness = async (smiles: string): Promise<number> => {
  try {
    const response = await api.get<number>(`${CHEM_URL}/nplikeness/score`, {
      params: { smiles },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to calculate NP likeness score: ${message}`);
  }
};

/**
 * Calculate Tanimoto similarity between two or more molecules
 */
const calculateTanimotoSimilarity = async (
  smiles: string,
  toolkit = "rdkit",
  fingerprinter = "ECFP",
  nBits = 2048,
  radius = 2
): Promise<number | string> => {
  try {
    const response = await api.get<number | string>(`${CHEM_URL}/tanimoto`, {
      params: { smiles, toolkit, fingerprinter, nBits, radius },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to calculate Tanimoto similarity: ${message}`);
  }
};

/**
 * Preprocess a molecule for COCONUT database submissions
 */
const coconutPreprocessing = async (
  smiles: string,
  _3d_mol = false,
  descriptors = false
): Promise<CoconutPreprocessingResult> => {
  try {
    const response = await api.get<CoconutPreprocessingResult>(
      `${CHEM_URL}/coconut/pre-processing`,
      {
        params: { smiles, _3d_mol, descriptors },
      }
    );
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to preprocess for COCONUT: ${message}`);
  }
};

/**
 * Generate functional groups for a molecule using Ertl's algorithm
 */
const generateFunctionalGroups = async (smiles: string): Promise<unknown[]> => {
  try {
    const response = await api.get<unknown[]>(`${CHEM_URL}/ertlfunctionalgroup`, {
      params: { smiles },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate functional groups: ${message}`);
  }
};

/**
 * Generate standardized tautomer for a molecule
 */
const generateStandardizedTautomer = async (smiles: string): Promise<string> => {
  try {
    const response = await api.get<string>(`${CHEM_URL}/standarizedTautomer`, {
      params: { smiles },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate standardized tautomer: ${message}`);
  }
};

/**
 * Fix radicals (single electrons) in molecules
 */
const fixRadicals = async (smiles: string): Promise<FixRadicalsResult> => {
  try {
    const response = await api.get<FixRadicalsResult>(`${CHEM_URL}/fixRadicals`, {
      params: { smiles },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to fix radicals: ${message}`);
  }
};

/**
 * Standardize a molblock using the ChEMBL curation pipeline
 */
const standardizeMolblock = async (molblock: string): Promise<StandardizeResult> => {
  try {
    const response = await api.post<StandardizeResult>(`${CHEM_URL}/standardize`, molblock, {
      headers: {
        "Content-Type": "text/plain",
      },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to standardize molblock: ${message}`);
  }
};

/**
 * Find chemical structure by name, formula, or identifier using the PubChem database
 */
const lookupPubChem = async (identifier: string): Promise<PubChemLookupResult> => {
  try {
    const response = await api.get<PubChemLookupResult>(`${CHEM_URL}/pubchem/smiles`, {
      params: { identifier },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to look up structure in PubChem: ${message}`);
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
  fixRadicals,
  getDescriptors,
  lookupPubChem,
  standardizeMolblock,
};

// Assemble all functions into a service object
const chemService = {
  generateStereoisomers,
  calculateDescriptors,
  getDescriptors, // Include the alias
  calculateMultipleDescriptors,
  generateHOSECodes,
  checkStructureErrors,
  calculateNPLikeness,
  calculateTanimotoSimilarity,
  coconutPreprocessing,
  generateFunctionalGroups,
  generateStandardizedTautomer,
  fixRadicals,
  lookupPubChem,
  standardizeMolblock,
};

export default chemService;
