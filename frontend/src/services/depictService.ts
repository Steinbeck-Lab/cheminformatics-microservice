// This module provides functions to generate 2D and 3D depictions of molecules
import api from "./api";
import type {
  DepictionOptions,
  EnhancedDepictionOptions,
  BatchDepictionResult,
  ParsedSmilesEntry,
} from "../types/api";

const DEPICT_URL = "/depict";

/**
 * Generate an enhanced 2D depiction of a molecule with advanced features
 */
export const generate2DDepictionEnhanced = async (
  smiles: string,
  options: EnhancedDepictionOptions = {}
): Promise<string> => {
  const {
    width = 512,
    height = 512,
    rotate = 0,
    CIP = true,
    unicolor = false,
    highlight = "",
    atomIds = null,
    showAtomNumbers = false,
    hydrogen_display = "Smart",
    abbreviate = "off",
    dative = "metals",
    multicenter = "provided",
    annotate = "none",
    style = "cow",
    donuts = false,
    arrow = "",
    alignrxnmap = true,
    showtitle = false,
    title = null,
    bgcolor = null,
    fgcolor = null,
    zoom = 1.0,
    ratio = 1.0,
    flip = false,
    anon = false,
    smalim = 100,
    svgunits = "px",
    perceive_radicals = false,
    apply_mdl_highlighting = true,
  } = options;

  try {
    const params: Record<string, unknown> = {
      smiles,
      width,
      height,
      rotate,
      CIP,
      unicolor,
      showAtomNumbers,
      hydrogen_display,
      abbreviate,
      dative,
      multicenter,
      annotate,
      style,
      donuts,
      arrow,
      alignrxnmap,
      showtitle,
      zoom,
      ratio,
      flip,
      anon,
      smalim,
      svgunits,
      perceive_radicals,
      apply_mdl_highlighting,
    };

    if (highlight) params.highlight = highlight;
    if (title) params.title = title;
    if (bgcolor) params.bgcolor = bgcolor;
    if (fgcolor) params.fgcolor = fgcolor;

    // Convert atomIds to comma-separated string if provided
    if (atomIds) {
      if (Array.isArray(atomIds)) {
        const flatIds = (atomIds as number[][]).flat();
        params.atomIds = flatIds.join(",");
      }
    }

    const response = await api.get<string>(`${DEPICT_URL}/2D_enhanced`, {
      params,
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate enhanced 2D depiction: ${message}`);
  }
};

/**
 * Generate a 2D depiction of a molecule
 */
export const generate2DDepiction = async (
  smiles: string,
  options: DepictionOptions = {}
): Promise<string> => {
  const {
    toolkit = "rdkit",
    width = 512,
    height = 512,
    rotate = 0,
    CIP = false,
    unicolor = false,
    highlight = "COSN",
    atomIds = null,
    showAtomNumbers = false,
    hydrogen_display = "Smart",
  } = options;

  try {
    const params: Record<string, unknown> = {
      smiles,
      toolkit,
      width,
      height,
      rotate,
      CIP,
      unicolor,
      highlight,
      showAtomNumbers,
    };

    // Add hydrogen_display parameter for CDK toolkit
    if (toolkit === "cdk" && hydrogen_display) {
      params.hydrogen_display = hydrogen_display;
    }

    // Convert atomIds to comma-separated string if provided
    if (atomIds) {
      if (Array.isArray(atomIds)) {
        // Flatten if it's an array of arrays
        const flatIds = (atomIds as number[][]).flat();
        params.atomIds = flatIds.join(",");
      }
    }

    const response = await api.get<string>(`${DEPICT_URL}/2D`, {
      params,
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate 2D depiction: ${message}`);
  }
};

/**
 * Generate a 3D depiction of a molecule
 */
export const generate3DDepiction = async (
  smiles: string,
  toolkit = "openbabel"
): Promise<string> => {
  try {
    const response = await api.get<string>(`${DEPICT_URL}/3D`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate 3D depiction: ${message}`);
  }
};

/**
 * Get the URL for an enhanced 2D depiction image with advanced features
 */
export const get2DDepictionUrlEnhanced = (
  smiles: string,
  options: EnhancedDepictionOptions = {}
): string => {
  const {
    width = 512,
    height = 512,
    rotate = 0,
    CIP = true,
    unicolor = false,
    highlight = "",
    atomIds = null,

    format: _format = "svg",
    showAtomNumbers = false,
    hydrogen_display = "Smart",
    abbreviate = "off",
    dative = "metals",
    multicenter = "provided",
    annotate = "none",
    style = "cow",
    donuts = false,
    arrow = "",
    alignrxnmap = true,
    showtitle = false,
    title = null,
    bgcolor = null,
    fgcolor = null,
    zoom = 1.0,
    ratio = 1.0,
    flip = false,
    anon = false,
    smalim = 100,
    svgunits = "px",
    perceive_radicals = false,
    apply_mdl_highlighting = true,
  } = options;

  const baseUrl = api.defaults.baseURL || "";
  const url = new URL(`${baseUrl}${DEPICT_URL}/2D_enhanced`);

  url.searchParams.append("smiles", smiles);
  url.searchParams.append("width", String(width));
  url.searchParams.append("height", String(height));
  url.searchParams.append("rotate", String(rotate));
  url.searchParams.append("CIP", String(CIP));
  url.searchParams.append("unicolor", String(unicolor));
  url.searchParams.append("showAtomNumbers", String(showAtomNumbers));
  url.searchParams.append("hydrogen_display", hydrogen_display);
  url.searchParams.append("abbreviate", abbreviate);
  url.searchParams.append("dative", dative);
  url.searchParams.append("multicenter", multicenter);
  url.searchParams.append("annotate", annotate);
  url.searchParams.append("style", style);
  url.searchParams.append("donuts", String(donuts));
  url.searchParams.append("arrow", arrow);
  url.searchParams.append("alignrxnmap", String(alignrxnmap));
  url.searchParams.append("showtitle", String(showtitle));
  url.searchParams.append("zoom", String(zoom));
  url.searchParams.append("ratio", String(ratio));
  url.searchParams.append("flip", String(flip));
  url.searchParams.append("anon", String(anon));
  url.searchParams.append("smalim", String(smalim));
  url.searchParams.append("svgunits", svgunits);
  url.searchParams.append("perceive_radicals", String(perceive_radicals));
  url.searchParams.append("apply_mdl_highlighting", String(apply_mdl_highlighting));

  if (highlight) {
    url.searchParams.append("highlight", highlight);
  }

  if (title) {
    url.searchParams.append("title", title);
  }

  if (bgcolor) {
    url.searchParams.append("bgcolor", bgcolor);
  }

  if (fgcolor) {
    url.searchParams.append("fgcolor", fgcolor);
  }

  // Add atomIds if provided
  if (atomIds) {
    if (Array.isArray(atomIds)) {
      const flatIds = (atomIds as number[][]).flat();
      url.searchParams.append("atomIds", flatIds.join(","));
    }
  }

  return url.toString();
};

/**
 * Get the URL for a 2D depiction image
 */
export const get2DDepictionUrl = (
  smiles: string,
  options: DepictionOptions & { format?: string } = {}
): string => {
  const {
    width = 512,
    height = 512,
    toolkit = "rdkit",
    rotate = 0,
    CIP = false,
    unicolor = false,
    highlight = "",
    atomIds = null,
    format = "svg",
    showAtomNumbers = false,
    hydrogen_display = "Smart",
  } = options;

  const baseUrl = api.defaults.baseURL || "";
  let url: URL;

  // Support format parameter in the URL
  if (format && format !== "svg") {
    url = new URL(`${baseUrl}${DEPICT_URL}/2D/${format}`);
  } else {
    url = new URL(`${baseUrl}${DEPICT_URL}/2D`);
  }

  url.searchParams.append("smiles", smiles);
  url.searchParams.append("width", String(width));
  url.searchParams.append("height", String(height));
  url.searchParams.append("toolkit", toolkit);
  url.searchParams.append("rotate", String(rotate));
  url.searchParams.append("CIP", String(CIP));
  url.searchParams.append("unicolor", String(unicolor));
  url.searchParams.append("showAtomNumbers", String(showAtomNumbers));

  // Add hydrogen_display parameter for CDK toolkit
  if (toolkit === "cdk" && hydrogen_display) {
    url.searchParams.append("hydrogen_display", hydrogen_display);
  }

  if (highlight) {
    url.searchParams.append("highlight", highlight);
  }

  // Add atomIds if provided
  if (atomIds) {
    if (Array.isArray(atomIds)) {
      // Flatten if it's an array of arrays
      const flatIds = (atomIds as number[][]).flat();
      url.searchParams.append("atomIds", flatIds.join(","));
    }
  }

  return url.toString();
};

/**
 * Get the URL for an enhanced 2D depiction with advanced CDK features
 */
export const get2DDepictionEnhancedUrl = (
  smiles: string,
  options: EnhancedDepictionOptions = {}
): string => {
  const {
    width = 512,
    height = 512,
    rotate = 0,
    CIP = true,
    unicolor = false,
    highlight = "",
    atomIds = null,
    showAtomNumbers = false,
    hydrogen_display = "Smart",
    abbreviate = "off",
    dative = "metals",
    multicenter = "provided",
    annotate = "none",
    style = "cow",
    donuts = false,
    arrow = "",
    alignrxnmap = true,
    showtitle = false,
    title = null,
    bgcolor = null,
    fgcolor = null,
    zoom = 1.0,
    ratio = 1.0,
    flip = false,
    anon = false,
    smalim = 100,
    svgunits = "px",
    perceive_radicals = false,
    apply_mdl_highlighting = true,
  } = options;

  const baseUrl = api.defaults.baseURL || "";
  const url = new URL(`${baseUrl}${DEPICT_URL}/2D_enhanced`);

  url.searchParams.append("smiles", smiles);
  url.searchParams.append("width", String(width));
  url.searchParams.append("height", String(height));
  url.searchParams.append("rotate", String(rotate));
  url.searchParams.append("CIP", String(CIP));
  url.searchParams.append("unicolor", String(unicolor));
  url.searchParams.append("showAtomNumbers", String(showAtomNumbers));
  url.searchParams.append("hydrogen_display", hydrogen_display);
  url.searchParams.append("abbreviate", abbreviate);
  url.searchParams.append("dative", dative);
  url.searchParams.append("multicenter", multicenter);
  url.searchParams.append("annotate", annotate);
  url.searchParams.append("style", style);
  url.searchParams.append("donuts", String(donuts));
  url.searchParams.append("arrow", arrow);
  url.searchParams.append("alignrxnmap", String(alignrxnmap));
  url.searchParams.append("showtitle", String(showtitle));
  url.searchParams.append("zoom", String(zoom));
  url.searchParams.append("ratio", String(ratio));
  url.searchParams.append("flip", String(flip));
  url.searchParams.append("anon", String(anon));
  url.searchParams.append("smalim", String(smalim));
  url.searchParams.append("svgunits", svgunits);
  url.searchParams.append("perceive_radicals", String(perceive_radicals));
  url.searchParams.append("apply_mdl_highlighting", String(apply_mdl_highlighting));

  if (highlight) {
    url.searchParams.append("highlight", highlight);
  }

  if (title) {
    url.searchParams.append("title", title);
  }

  if (bgcolor) {
    url.searchParams.append("bgcolor", bgcolor);
  }

  if (fgcolor) {
    url.searchParams.append("fgcolor", fgcolor);
  }

  // Add atomIds if provided
  if (atomIds) {
    if (Array.isArray(atomIds)) {
      const flatIds = (atomIds as number[][]).flat();
      url.searchParams.append("atomIds", flatIds.join(","));
    }
  }

  return url.toString();
};

/**
 * Generate batch 2D depictions for multiple SMILES strings
 */
export const generateBatchDepictions = async (
  smilesList: string[],
  options: DepictionOptions = {}
): Promise<BatchDepictionResult[]> => {
  try {
    const results: BatchDepictionResult[] = [];

    for (const smiles of smilesList) {
      const svg = await generate2DDepiction(smiles, options);
      results.push({
        smiles,
        svg,
      });
    }

    return results;
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate batch depictions: ${message}`);
  }
};

interface DepictionEntry {
  smiles: string;
  title: string;
  svg: string;
}

/**
 * Download multiple depictions as a ZIP file
 */
export const downloadDepictionsAsZip = async (
  depictions: DepictionEntry[],
  format = "svg",
  options: DepictionOptions = {}
): Promise<Blob> => {
  try {
    const JSZip = (await import("jszip")).default;
    const zip = new JSZip();

    for (const depiction of depictions) {
      const { smiles, title, svg } = depiction;
      const filename = `${title.replace(/[^a-z0-9]/gi, "_").toLowerCase()}.${format}`;

      if (format === "svg") {
        zip.file(filename, svg);
      } else {
        // For other formats, we need to fetch from the API
        const url = get2DDepictionUrl(smiles, { ...options, format });
        const response = await fetch(url);
        const blob = await response.blob();
        zip.file(filename, blob);
      }
    }

    return await zip.generateAsync({ type: "blob" });
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to create ZIP file: ${message}`);
  }
};

/**
 * Highlight a substructure in a molecule
 */
export const highlightSubstructure = async (
  smiles: string,
  substructure: string,
  options: DepictionOptions = {}
): Promise<string> => {
  const depictionOptions = {
    ...options,
    highlight: substructure,
  };

  return generate2DDepiction(smiles, depictionOptions);
};

/**
 * Generate a color-coded depiction based on atom properties
 */
export const generateColorCodedDepiction = async (
  smiles: string,
  colorBy = "element",
  options: DepictionOptions = {}
): Promise<string> => {
  const depictionOptions = {
    ...options,
    unicolor: colorBy === "none",
  };

  return generate2DDepiction(smiles, depictionOptions);
};

/**
 * Generate multiple depictions with different settings
 */
export const generateMultipleDepictions = async (
  smiles: string,
  optionsArray: DepictionOptions[] = [{}]
): Promise<string[]> => {
  try {
    const promises = optionsArray.map((options) => generate2DDepiction(smiles, options));

    return Promise.all(promises);
  } catch (error) {
    const message = error instanceof Error ? error.message : "Unknown error";
    throw new Error(`Failed to generate multiple depictions: ${message}`);
  }
};

/**
 * Parse SMILES strings from a text input with optional titles
 */
export const parseSmilesInput = (text: string): ParsedSmilesEntry[] => {
  return text
    .split("\n")
    .map((rawLine) => {
      const line = rawLine.trim();
      // Skip empty lines and comments
      if (line.length === 0 || line.startsWith("#")) {
        return null;
      }

      // Extract title if present (after space or tab)
      const spaceIndex = line.search(/[\s\t]/);
      const smiles = spaceIndex > 0 ? line.substring(0, spaceIndex) : line;
      const title = spaceIndex > 0 ? line.substring(spaceIndex + 1).trim() : `Molecule`;

      return { smiles, title };
    })
    .filter((item): item is ParsedSmilesEntry => item !== null);
};

// Assemble all functions into a service object
const depictService = {
  generate2DDepiction,
  generate2DDepictionEnhanced,
  generate3DDepiction,
  get2DDepictionUrl,
  get2DDepictionUrlEnhanced,
  highlightSubstructure,
  generateColorCodedDepiction,
  generateMultipleDepictions,
  generateBatchDepictions,
  downloadDepictionsAsZip,
  parseSmilesInput,
};

export default depictService;
