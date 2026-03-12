// This file contains functions to interact with the tools API for chemical structure generation, sugar detection, and filtering.
import api from "./api";
import type {
  SugarRemovalOptions,
  ExtractionOptions,
  ChemicalFilterOptions,
  ClassyFireResult,
  StandardizeResult,
} from "../types/api";

const TOOLS_URL = "/tools";

/**
 * Generate chemical structures based on a molecular formula
 */
export const generateStructures = async (molecularFormula: string): Promise<string[]> => {
  try {
    const response = await api.get<string[]>(`${TOOLS_URL}/generate-structures`, {
      params: { molecular_formula: molecularFormula },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate structures: ${message}`);
  }
};

/**
 * Get information about sugar moieties in a molecule with full parameter support
 */
export const getSugarInfo = async (
  smiles: string,
  options: SugarRemovalOptions = {}
): Promise<string> => {
  try {
    const params = {
      smiles,
      gly_bond: options.gly_bond ?? false,
      oxygen_atoms: options.oxygen_atoms ?? true,
      oxygen_atoms_threshold: options.oxygen_atoms_threshold ?? 0.5,
      linear_sugars_in_rings: options.linear_sugars_in_rings ?? false,
      linear_sugars_min_size: options.linear_sugars_min_size ?? 4,
      linear_sugars_max_size: options.linear_sugars_max_size ?? 7,
      linear_acidic_sugars: options.linear_acidic_sugars ?? false,
      spiro_sugars: options.spiro_sugars ?? false,
      keto_sugars: options.keto_sugars ?? false,
    };
    const response = await api.get<string>(`${TOOLS_URL}/sugars-info`, { params });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to get sugar information: ${message}`);
  }
};

/**
 * Detect sugar moieties in a molecule (alias for getSugarInfo)
 */
export const detectSugars = getSugarInfo;

/**
 * Remove linear sugars from a molecule with full parameter support
 */
export const removeLinearSugars = async (
  smiles: string,
  options: SugarRemovalOptions = {}
): Promise<string> => {
  try {
    const params = {
      smiles,
      only_terminal: options.only_terminal ?? true,
      preservation_mode: options.preservation_mode ?? 2,
      preservation_threshold: options.preservation_threshold ?? 5,
      linear_sugars_in_rings: options.linear_sugars_in_rings ?? false,
      linear_sugars_min_size: options.linear_sugars_min_size ?? 4,
      linear_sugars_max_size: options.linear_sugars_max_size ?? 7,
      linear_acidic_sugars: options.linear_acidic_sugars ?? false,
      mark_attach_points: options.mark_attach_points ?? false,
    };
    const response = await api.get<string>(`${TOOLS_URL}/remove-linear-sugars`, {
      params,
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to remove linear sugars: ${message}`);
  }
};

/**
 * Remove circular sugars from a molecule with full parameter support
 */
export const removeCircularSugars = async (
  smiles: string,
  options: SugarRemovalOptions = {}
): Promise<string> => {
  try {
    const params = {
      smiles,
      gly_bond: options.gly_bond ?? false,
      only_terminal: options.only_terminal ?? true,
      preservation_mode: options.preservation_mode ?? 2,
      preservation_threshold: options.preservation_threshold ?? 5,
      oxygen_atoms: options.oxygen_atoms ?? true,
      oxygen_atoms_threshold: options.oxygen_atoms_threshold ?? 0.5,
      spiro_sugars: options.spiro_sugars ?? false,
      keto_sugars: options.keto_sugars ?? false,
      mark_attach_points: options.mark_attach_points ?? false,
    };
    const response = await api.get<string>(`${TOOLS_URL}/remove-circular-sugars`, {
      params,
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to remove circular sugars: ${message}`);
  }
};

/**
 * Remove both linear and circular sugars from a molecule with full parameter support
 */
export const removeAllSugars = async (
  smiles: string,
  options: SugarRemovalOptions = {}
): Promise<string> => {
  try {
    const params = {
      smiles,
      gly_bond: options.gly_bond ?? false,
      only_terminal: options.only_terminal ?? true,
      preservation_mode: options.preservation_mode ?? 2,
      preservation_threshold: options.preservation_threshold ?? 5,
      oxygen_atoms: options.oxygen_atoms ?? true,
      oxygen_atoms_threshold: options.oxygen_atoms_threshold ?? 0.5,
      linear_sugars_in_rings: options.linear_sugars_in_rings ?? false,
      linear_sugars_min_size: options.linear_sugars_min_size ?? 4,
      linear_sugars_max_size: options.linear_sugars_max_size ?? 7,
      linear_acidic_sugars: options.linear_acidic_sugars ?? false,
      spiro_sugars: options.spiro_sugars ?? false,
      keto_sugars: options.keto_sugars ?? false,
      mark_attach_points: options.mark_attach_points ?? false,
    };
    const response = await api.get<string>(`${TOOLS_URL}/remove-sugars`, { params });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to remove all sugars: ${message}`);
  }
};

/**
 * Extract aglycone and sugar moieties from a molecule with full parameter support
 */
export const extractAglyconeAndSugars = async (
  smiles: string,
  options: ExtractionOptions = {}
): Promise<string[]> => {
  try {
    const params = {
      smiles,
      extract_circular_sugars: options.extract_circular_sugars ?? true,
      extract_linear_sugars: options.extract_linear_sugars ?? false,
      gly_bond: options.gly_bond ?? false,
      only_terminal: options.only_terminal ?? true,
      preservation_mode: options.preservation_mode ?? 2,
      preservation_threshold: options.preservation_threshold ?? 5,
      oxygen_atoms: options.oxygen_atoms ?? true,
      oxygen_atoms_threshold: options.oxygen_atoms_threshold ?? 0.5,
      linear_sugars_in_rings: options.linear_sugars_in_rings ?? false,
      linear_sugars_min_size: options.linear_sugars_min_size ?? 4,
      linear_sugars_max_size: options.linear_sugars_max_size ?? 7,
      linear_acidic_sugars: options.linear_acidic_sugars ?? false,
      spiro_sugars: options.spiro_sugars ?? false,
      keto_sugars: options.keto_sugars ?? false,
      mark_attach_points: options.mark_attach_points ?? false,
      post_process_sugars: options.post_process_sugars ?? false,
      limit_post_process_by_size: options.limit_post_process_by_size ?? false,
    };
    const response = await api.get<string[]>(`${TOOLS_URL}/extract-aglycone-and-sugars`, {
      params,
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to extract aglycone and sugars: ${message}`);
  }
};

/**
 * General sugar removal function
 */
export const removeSugars = async (
  smiles: string,
  type = "all",
  options: SugarRemovalOptions = {}
): Promise<string> => {
  try {
    switch (type) {
      case "linear":
        return await removeLinearSugars(smiles, options);
      case "circular":
        return await removeCircularSugars(smiles, options);
      case "all":
      default:
        return await removeAllSugars(smiles, options);
    }
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to remove sugars: ${message}`);
  }
};

/**
 * Get atom indices of aglycone and sugar moieties
 */
export const getAglyconeAndSugarIndices = async (
  smiles: string,
  options: ExtractionOptions = {}
): Promise<number[][]> => {
  try {
    const params = {
      smiles,
      extract_circular_sugars: options.extract_circular_sugars ?? true,
      extract_linear_sugars: options.extract_linear_sugars ?? false,
      gly_bond: options.gly_bond ?? false,
      only_terminal: options.only_terminal ?? false, // Set to false to get all sugars for highlighting
      preservation_mode: options.preservation_mode ?? 2,
      preservation_threshold: options.preservation_threshold ?? 5,
      oxygen_atoms: options.oxygen_atoms ?? true,
      oxygen_atoms_threshold: options.oxygen_atoms_threshold ?? 0.5,
      linear_sugars_in_rings: options.linear_sugars_in_rings ?? false,
      linear_sugars_min_size: options.linear_sugars_min_size ?? 4,
      linear_sugars_max_size: options.linear_sugars_max_size ?? 7,
      linear_acidic_sugars: options.linear_acidic_sugars ?? false,
      spiro_sugars: options.spiro_sugars ?? false,
      keto_sugars: options.keto_sugars ?? false,
      mark_attach_points: options.mark_attach_points ?? false,
      post_process_sugars: options.post_process_sugars ?? false,
      limit_post_process_by_size: options.limit_post_process_by_size ?? false,
    };
    const response = await api.get<number[][]>(`${TOOLS_URL}/get-aglycone-and-sugar-indices`, {
      params,
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to get aglycone and sugar indices: ${message}`);
  }
};

/**
 * Apply multiple chemical filters to a list of molecules
 */
export const applyChemicalFilters = async (
  smilesList: string,
  filterOptions: ChemicalFilterOptions = {}
): Promise<string[]> => {
  const {
    pains = true,
    lipinski = true,
    veber = true,
    reos = true,
    ghose = true,
    ruleofthree = true,
    qedscore = "0-10",
    sascore = "0-10",
    nplikeness = "0-10",
  } = filterOptions;

  try {
    const response = await api.post<string[]>(`/chem/all_filters`, smilesList, {
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
      },
      headers: {
        "Content-Type": "text/plain",
      },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to apply chemical filters: ${message}`);
  }
};

/**
 * Apply chemical filters with detailed violation information
 */
export const applyChemicalFiltersDetailed = async (
  smilesList: string,
  filterOptions: ChemicalFilterOptions = {}
): Promise<unknown> => {
  const {
    pains = true,
    lipinski = true,
    veber = true,
    reos = true,
    ghose = true,
    ruleofthree = true,
    qedscore = "0-1",
    sascore = "0-10",
    nplikeness = "-5-5",
    filterOperator = "OR",
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
        filterOperator,
      },
      headers: {
        "Content-Type": "text/plain",
        Accept: "application/json",
      },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to apply detailed chemical filters: ${message}`);
  }
};

/**
 * Standardize a molecule using the ChEMBL curation pipeline
 */
export const standardizeMolecule = async (molblock: string): Promise<StandardizeResult> => {
  try {
    const response = await api.post<StandardizeResult>(`/chem/standardize`, molblock, {
      headers: {
        "Content-Type": "text/plain",
      },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to standardize molecule: ${message}`);
  }
};

/**
 * Get ClassyFire classification for a molecule
 */
export const classifyMolecule = async (smiles: string): Promise<ClassyFireResult> => {
  try {
    const response = await api.get<ClassyFireResult>(`/chem/classyfire/classify`, {
      params: { smiles },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to classify molecule: ${message}`);
  }
};

/**
 * Get ClassyFire classification results
 */
export const getClassificationResults = async (jobId: string): Promise<ClassyFireResult> => {
  try {
    const response = await api.get<ClassyFireResult>(`/chem/classyfire/${jobId}/result`);
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to get classification results: ${message}`);
  }
};

// Assemble all functions into a service object
const toolsService = {
  generateStructures,
  getSugarInfo,
  detectSugars,
  removeLinearSugars,
  removeCircularSugars,
  removeAllSugars,
  removeSugars,
  extractAglyconeAndSugars,
  getAglyconeAndSugarIndices,
  applyChemicalFilters,
  applyChemicalFiltersDetailed,
  standardizeMolecule,
  classifyMolecule,
  getClassificationResults,
};

export default toolsService;
