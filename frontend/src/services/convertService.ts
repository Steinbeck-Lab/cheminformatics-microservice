// This module provides functions to convert chemical representations using a REST API.
import api from "./api";
import { AxiosError } from "axios";
import type {
  MultipleFormatsResult,
  CdxConversionResult,
  XYZBatchConversionResult,
  XYZConversionOptions,
} from "../types/api";

const CONVERT_URL = "/convert";

/**
 * Generate 2D coordinates for a molecule
 */
export const generate2DCoordinates = async (smiles: string, toolkit = "cdk"): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/mol2D`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate 2D coordinates: ${message}`);
  }
};

/**
 * Generate 3D coordinates for a molecule
 */
export const generate3DCoordinates = async (
  smiles: string,
  toolkit = "openbabel"
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/mol3D`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate 3D coordinates: ${message}`);
  }
};

/**
 * Generate SMILES from IUPAC name or SELFIES
 */
export const generateSMILES = async (
  inputText: string,
  representation = "iupac",
  converter = "opsin"
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/smiles`, {
      params: { input_text: inputText, representation, converter },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate SMILES: ${message}`);
  }
};

/**
 * Generate canonical SMILES
 */
export const generateCanonicalSMILES = async (smiles: string, toolkit = "cdk"): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/canonicalsmiles`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate canonical SMILES: ${message}`);
  }
};

/**
 * Generate CXSMILES (ChemAxon Extended SMILES)
 */
export const generateCXSMILES = async (smiles: string, toolkit = "cdk"): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/cxsmiles`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate CXSMILES: ${message}`);
  }
};

/**
 * Generate InChI (IUPAC International Chemical Identifier)
 */
export const generateInChI = async (smiles: string, toolkit = "cdk"): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/inchi`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate InChI: ${message}`);
  }
};

/**
 * Generate InChI Key
 */
export const generateInChIKey = async (smiles: string, toolkit = "cdk"): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/inchikey`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate InChI Key: ${message}`);
  }
};

/**
 * Generate SELFIES (Self-Referencing Embedded Strings)
 */
export const generateSELFIES = async (smiles: string): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/selfies`, {
      params: { smiles },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate SELFIES: ${message}`);
  }
};

/**
 * Generate SMARTS pattern
 */
export const generateSMARTS = async (smiles: string, toolkit = "rdkit"): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/smarts`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate SMARTS: ${message}`);
  }
};

/**
 * Convert MOL/SDF block to SMILES
 */
export const molblockToSMILES = async (molblock: string, toolkit = "cdk"): Promise<string> => {
  try {
    // Get the base URL from the api instance or use default
    const baseURL = api.defaults.baseURL;

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
        const errorData = (await response.json()) as { detail?: string };
        errorMsg = errorData.detail || errorMsg;
      } catch (_jsonError) {
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
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to convert MOL/SDF to SMILES: ${message}`);
  }
};

/**
 * Generate multiple formats at once
 */
export const generateMultipleFormats = async (
  smiles: string,
  toolkit = "cdk"
): Promise<MultipleFormatsResult> => {
  try {
    const response = await api.get<MultipleFormatsResult>(`${CONVERT_URL}/formats`, {
      params: { smiles, toolkit },
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate multiple formats: ${message}`);
  }
};

/**
 * Convert a CDX or CDXML file to a MOL block
 */
export const convertCDXToMol = async (file: File): Promise<string> => {
  const formData = new FormData();
  formData.append("file", file);

  try {
    const response = await api.post<CdxConversionResult>(`${CONVERT_URL}/cdx-to-mol`, formData, {
      headers: { "Content-Type": "multipart/form-data" },
    });
    return response.data.molblock;
  } catch (error) {
    const axiosError = error as AxiosError<{ detail?: string }>;
    const detail = axiosError.response?.data?.detail;
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(detail || `Failed to convert CDX file: ${message}`);
  }
};

/**
 * Convert an XYZ-coordinate block to SMILES, InChI, InChIKey, MOL, and SDF.
 *
 * Uses the RDKit xyz2mol algorithm by default (charge-aware bond perception).
 * Pass ``toolkit: "openbabel"`` for the distance-based fallback (best for
 * neutral species — OpenBabel ignores ``charge`` and ``useHueckel``).
 */
export const convertXYZ = async (
  xyz: string,
  options: XYZConversionOptions = {}
): Promise<XYZBatchConversionResult> => {
  const { charge = 0, useHueckel = false, toolkit = "rdkit" } = options;

  const baseURL = api.defaults.baseURL;
  const params = new URLSearchParams({
    charge: String(charge),
    use_huckel: String(useHueckel),
    toolkit,
  });

  try {
    const response = await fetch(`${baseURL}${CONVERT_URL}/xyz?${params.toString()}`, {
      method: "POST",
      headers: {
        "Content-Type": "text/plain; charset=utf-8",
        Accept: "application/json",
      },
      body: xyz,
    });

    if (!response.ok) {
      let detail = `Error ${response.status}: ${response.statusText}`;
      try {
        const errorData = (await response.json()) as { detail?: string };
        if (errorData.detail) detail = errorData.detail;
      } catch (_jsonError) {
        // Response was not JSON; fall back to status text.
      }
      throw new Error(detail);
    }

    return (await response.json()) as XYZBatchConversionResult;
  } catch (error) {
    if (error instanceof Error) throw error;
    throw new Error("Failed to convert XYZ block: Unknown error");
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
  convertCDXToMol,
  convertXYZ,
};

export default convertService;
