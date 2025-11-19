// This module provides functions to generate 2D and 3D depictions of molecules
import api from "./api";

const DEPICT_URL = "/depict";

/**
 * Generate an enhanced 2D depiction of a molecule with advanced features
 * @param {string} smiles - SMILES or CXSMILES string
 * @param {Object} options - Enhanced depiction options
 * @param {number} options.width - Width of the generated image
 * @param {number} options.height - Height of the generated image
 * @param {number} options.rotate - Rotation angle in degrees
 * @param {boolean} options.CIP - Include Cahn-Ingold-Prelog stereochemistry
 * @param {boolean} options.unicolor - Use a single color for the molecule (deprecated, use style)
 * @param {string} options.highlight - SMARTS pattern to highlight atoms/bonds
 * @param {Array<number>|Array<Array<number>>} options.atomIds - Atom indices to highlight
 * @param {boolean} options.showAtomNumbers - Show atom numbers on the depiction
 * @param {string} options.hydrogen_display - Hydrogen display mode (Provided, Minimal, Explicit, Stereo, Smart)
 * @param {string} options.abbreviate - Chemical abbreviation mode (off, groups, reagents, on)
 * @param {string} options.dative - Dative bond perception (never, metals, always)
 * @param {string} options.multicenter - Multicenter bond display style
 * @param {string} options.annotate - Annotation mode (none, number, bondnumber, mapidx, atomvalue, colmap, cip)
 * @param {string} options.style - Color scheme preset (cow, cob, cot, bow, bot, wob, nob)
 * @param {boolean} options.donuts - Use circle-in-ring for aromatic rings
 * @param {string} options.arrow - Reaction arrow type
 * @param {boolean} options.alignrxnmap - Align reaction mapped atoms
 * @param {boolean} options.showtitle - Display molecule/reaction title
 * @param {string} options.title - Optional title to display
 * @param {string} options.bgcolor - Custom background color (hex)
 * @param {string} options.fgcolor - Custom foreground color (hex)
 * @param {number} options.zoom - Zoom level (0.1 to 5.0)
 * @param {number} options.ratio - Bond thickness ratio (0.5 to 2.0)
 * @param {boolean} options.flip - Horizontally flip the structure
 * @param {boolean} options.anon - Use anonymous atom display
 * @param {number} options.smalim - SMARTS hit limit
 * @param {string} options.svgunits - SVG coordinate units (px, mm, cm, in)
 * @param {boolean} options.perceive_radicals - Detect and mark radicals
 * @param {boolean} options.apply_mdl_highlighting - Apply MDL V3000 highlighting
 * @returns {Promise<string>} - SVG depiction as text
 */
export const generate2DDepictionEnhanced = async (smiles, options = {}) => {
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
    const params = {
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
        const flatIds = atomIds.flat();
        params.atomIds = flatIds.join(",");
      }
    }

    const response = await api.get(`${DEPICT_URL}/2D_enhanced`, {
      params,
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate enhanced 2D depiction: ${error.message}`);
  }
};

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
 * @param {Array<number>|Array<Array<number>>} options.atomIds - Atom indices to highlight (single array or array of arrays for multiple substructures)
 * @param {boolean} options.showAtomNumbers - Show atom numbers on the depiction
 * @param {string} options.hydrogen_display - Hydrogen display mode for CDK (Provided, Minimal, Explicit, Stereo, Smart)
 * @returns {Promise<string>} - SVG depiction as text
 */
export const generate2DDepiction = async (smiles, options = {}) => {
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
    const params = {
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
        const flatIds = atomIds.flat();
        params.atomIds = flatIds.join(",");
      }
    }

    const response = await api.get(`${DEPICT_URL}/2D`, {
      params,
      responseType: "text",
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
export const generate3DDepiction = async (smiles, toolkit = "openbabel") => {
  try {
    const response = await api.get(`${DEPICT_URL}/3D`, {
      params: { smiles, toolkit },
      responseType: "text",
    });
    return response.data;
  } catch (error) {
    throw new Error(`Failed to generate 3D depiction: ${error.message}`);
  }
};

/**
 * Get the URL for an enhanced 2D depiction image with advanced features
 * @param {string} smiles - SMILES or CXSMILES string
 * @param {Object} options - Enhanced depiction options
 * @returns {string} - URL to the enhanced depiction image
 */
export const get2DDepictionUrlEnhanced = (smiles, options = {}) => {
  const {
    width = 512,
    height = 512,
    rotate = 0,
    CIP = true,
    unicolor = false,
    highlight = "",
    atomIds = null,
    format = "svg",
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
  url.searchParams.append("width", width);
  url.searchParams.append("height", height);
  url.searchParams.append("rotate", rotate);
  url.searchParams.append("CIP", CIP);
  url.searchParams.append("unicolor", unicolor);
  url.searchParams.append("showAtomNumbers", showAtomNumbers);
  url.searchParams.append("hydrogen_display", hydrogen_display);
  url.searchParams.append("abbreviate", abbreviate);
  url.searchParams.append("dative", dative);
  url.searchParams.append("multicenter", multicenter);
  url.searchParams.append("annotate", annotate);
  url.searchParams.append("style", style);
  url.searchParams.append("donuts", donuts);
  url.searchParams.append("arrow", arrow);
  url.searchParams.append("alignrxnmap", alignrxnmap);
  url.searchParams.append("showtitle", showtitle);
  url.searchParams.append("zoom", zoom);
  url.searchParams.append("ratio", ratio);
  url.searchParams.append("flip", flip);
  url.searchParams.append("anon", anon);
  url.searchParams.append("smalim", smalim);
  url.searchParams.append("svgunits", svgunits);
  url.searchParams.append("perceive_radicals", perceive_radicals);
  url.searchParams.append("apply_mdl_highlighting", apply_mdl_highlighting);

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
      const flatIds = atomIds.flat();
      url.searchParams.append("atomIds", flatIds.join(","));
    }
  }

  return url.toString();
};

/**
 * Get the URL for a 2D depiction image
 * @param {string} smiles - SMILES string
 * @param {Object} options - Depiction options
 * @param {Array<number>|Array<Array<number>>} options.atomIds - Atom indices to highlight
 * @param {string} options.hydrogen_display - Hydrogen display mode for CDK (Provided, Minimal, Explicit, Stereo, Smart)
 * @returns {string} - URL to the depiction image
 */
export const get2DDepictionUrl = (smiles, options = {}) => {
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
  let url;

  // Support format parameter in the URL
  if (format && format !== "svg") {
    url = new URL(`${baseUrl}${DEPICT_URL}/2D/${format}`);
  } else {
    url = new URL(`${baseUrl}${DEPICT_URL}/2D`);
  }

  url.searchParams.append("smiles", smiles);
  url.searchParams.append("width", width);
  url.searchParams.append("height", height);
  url.searchParams.append("toolkit", toolkit);
  url.searchParams.append("rotate", rotate);
  url.searchParams.append("CIP", CIP);
  url.searchParams.append("unicolor", unicolor);
  url.searchParams.append("showAtomNumbers", showAtomNumbers);

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
      const flatIds = atomIds.flat();
      url.searchParams.append("atomIds", flatIds.join(","));
    }
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
        svg,
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
export const downloadDepictionsAsZip = async (depictions, format = "svg", options = {}) => {
  // This function requires JSZip to be installed
  // npm install jszip
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
    highlight: substructure,
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
export const generateColorCodedDepiction = async (smiles, colorBy = "element", options = {}) => {
  // This is a placeholder - the actual API might not support this directly
  // but we can implement it as a convenience method that uses other endpoints
  const depictionOptions = {
    ...options,
    unicolor: colorBy === "none",
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
    const promises = optionsArray.map((options) => generate2DDepiction(smiles, options));

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
    .split("\n")
    .map((line) => {
      line = line.trim();
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
    .filter((item) => item !== null);
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
