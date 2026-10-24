// This module provides functions to convert chemical representations using a REST API.
import { failed } from "@/lib/error-messages";
import api from "./api";
import { AxiosError } from "axios";
import type {
  MultipleFormatsResult,
  CdxConversionResult,
  XYZBatchConversionResult,
  XYZConversionOptions,
  DetectFormatResult,
  IUPACNameResult,
} from "../types/api";

const CONVERT_URL = "/convert";

/**
 * Detect the chemical structure format of an input string.
 */
export const detectInputFormat = async (input: string): Promise<DetectFormatResult> => {
  const response = await api.get<DetectFormatResult>(`${CONVERT_URL}/detect-format`, {
    params: { input },
  });
  return response.data;
};

type ConvertParams = {
  autoDetect?: boolean;
};

const withAutoDetect = (params: Record<string, unknown>, autoDetect?: boolean) =>
  autoDetect ? { ...params, auto_detect: true } : params;

/**
 * Generate 2D coordinates for a molecule
 */
export const generate2DCoordinates = async (
  smiles: string,
  toolkit = "cdk",
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/mol2D`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate 2D coordinates", error);
  }
};

/**
 * Generate 3D coordinates for a molecule
 */
export const generate3DCoordinates = async (
  smiles: string,
  toolkit = "openbabel",
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/mol3D`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate 3D coordinates", error);
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
    throw failed("Failed to generate SMILES", error);
  }
};

/**
 * Generate canonical SMILES
 */
export const generateCanonicalSMILES = async (
  smiles: string,
  toolkit = "cdk",
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/canonicalsmiles`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate canonical SMILES", error);
  }
};

/**
 * Generate CXSMILES (ChemAxon Extended SMILES)
 */
export const generateCXSMILES = async (
  smiles: string,
  toolkit = "cdk",
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/cxsmiles`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate CXSMILES", error);
  }
};

/**
 * Generate InChI (IUPAC International Chemical Identifier)
 */
export const generateInChI = async (
  smiles: string,
  toolkit = "cdk",
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/inchi`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate InChI", error);
  }
};

/**
 * Generate InChI Key
 */
export const generateInChIKey = async (
  smiles: string,
  toolkit = "cdk",
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/inchikey`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate InChI Key", error);
  }
};

/**
 * Generate SELFIES (Self-Referencing Embedded Strings)
 */
export const generateSELFIES = async (
  smiles: string,
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/selfies`, {
      params: withAutoDetect({ smiles }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate SELFIES", error);
  }
};

/**
 * Generate SMARTS pattern
 */
export const generateSMARTS = async (
  smiles: string,
  toolkit = "rdkit",
  options: ConvertParams = {}
): Promise<string> => {
  try {
    const response = await api.get<string>(`${CONVERT_URL}/smarts`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate SMARTS", error);
  }
};

/**
 * Generate an IUPAC name from SMILES with Orthonym
 */
export const generateIUPACName = async (
  smiles: string,
  tier = "pin",
  options: ConvertParams = {}
): Promise<IUPACNameResult> => {
  try {
    const response = await api.get<IUPACNameResult>(`${CONVERT_URL}/iupac`, {
      params: withAutoDetect({ smiles, tier }, options.autoDetect),
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate IUPAC name", error);
  }
};

/** The `detail` of a JSON error response; undefined when the body is not JSON. */
async function detailOf(response: Response): Promise<string | undefined> {
  try {
    const errorData = (await response.json()) as { detail?: string };
    return errorData.detail;
  } catch {
    return undefined;
  }
}

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
      throw new Error(
        (await detailOf(response)) || `Error ${response.status}: ${response.statusText}`
      );
    }

    // Get the response as text and clean it up
    let result = await response.text();

    // Remove surrounding quotes if present
    if (result.startsWith('"') && result.endsWith('"')) {
      result = result.substring(1, result.length - 1);
    }

    return result;
  } catch (error) {
    throw failed("Failed to convert MOL/SDF to SMILES", error);
  }
};

/**
 * Generate multiple formats at once
 */
export const generateMultipleFormats = async (
  smiles: string,
  toolkit = "cdk",
  options: ConvertParams = {}
): Promise<MultipleFormatsResult> => {
  try {
    const response = await api.get<MultipleFormatsResult>(`${CONVERT_URL}/formats`, {
      params: withAutoDetect({ smiles, toolkit }, options.autoDetect),
    });
    return response.data;
  } catch (error) {
    throw failed("Failed to generate multiple formats", error);
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
    throw detail ? new Error(detail) : failed("Failed to convert CDX file", error);
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
      throw new Error(
        (await detailOf(response)) || `Error ${response.status}: ${response.statusText}`
      );
    }

    return (await response.json()) as XYZBatchConversionResult;
  } catch (error) {
    if (error instanceof Error) throw error;
    throw failed("Failed to convert XYZ block", error);
  }
};

// Assemble all functions into a service object
const convertService = {
  detectInputFormat,
  generate2DCoordinates,
  generate3DCoordinates,
  generateSMILES,
  generateCanonicalSMILES,
  generateCXSMILES,
  generateInChI,
  generateInChIKey,
  generateSELFIES,
  generateSMARTS,
  generateIUPACName,
  molblockToSMILES,
  generateMultipleFormats,
  convertCDXToMol,
  convertXYZ,
};

export default convertService;
