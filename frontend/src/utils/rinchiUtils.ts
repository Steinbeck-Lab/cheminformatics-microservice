/**
 * rinchiUtils.js - Utility functions for RInChI module loading and operations
 *
 * This file provides utilities for loading RInChI WASM modules from remote URLs and performing conversions
 * like reaction to RInChI, RInChI to reaction, and generating RInChIKeys.
 */

import { errorText } from "@/lib/error-messages";
import type {
  RInChIFileResult,
  RInChIKeyResult,
  RInChIResult,
  RInChIVersionConfig,
} from "@/types/api";

// Base URL for the GitHub Pages-hosted WASM files
const RINCHI_BASE_URL = "https://iupac-inchi.github.io/InChI-Web-Demo";

// Configuration for available RInChI modules with remote URLs
export const RINCHI_VERSIONS: Record<string, RInChIVersionConfig> = {
  1.1: {
    label: "RInChI 1.1-dev with Latest InChI",
    scriptSrc: `${RINCHI_BASE_URL}/rinchi/librinchi-1.1.js`,
    moduleName: "rinchiModule11",
    default: true,
  },
};

/** [version, label] pairs for a version picker. */
export const VERSION_OPTIONS: [string, string][] = Object.entries(RINCHI_VERSIONS).map(
  ([value, v]) => [value, v.label]
);

/** The version marked `default: true`. */
export const DEFAULT_VERSION =
  Object.keys(RINCHI_VERSIONS).find((key) => RINCHI_VERSIONS[key]?.default) || "1.1";

// Module references
const moduleInstances: Record<string, InchiWasmModule> = {};
const moduleLoadPromises: Record<string, Promise<InchiWasmModule>> = {};

/**
 * Loads a RInChI module from remote URL
 * @param {string} version - Version identifier ("1.1")
 * @returns {Promise<Object>} - Promise that resolves to the module instance
 */
export const loadRinchiModule = async (version: string): Promise<InchiWasmModule> => {
  // If we already have a loading promise for this version, return it
  if (moduleLoadPromises[version]) {
    return moduleLoadPromises[version];
  }

  // If we already have an instance for this version, return it
  if (moduleInstances[version]) {
    return Promise.resolve(moduleInstances[version]);
  }

  // Create a new loading promise
  moduleLoadPromises[version] = new Promise((resolve, reject) => {
    try {
      const versionConfig = RINCHI_VERSIONS[version];

      if (!versionConfig) {
        reject(new Error(`Unknown RInChI version: ${version}`));
        return;
      }

      // Check if the script is already loaded
      const existingScript = document.getElementById(`rinchi-script-${version}`);
      if (existingScript) {
        // Script exists, try to initialize the module
        const factory = window[versionConfig.moduleName];
        if (typeof factory === "function") {
          factory()
            .then((module) => {
              moduleInstances[version] = module;
              console.log(`RInChI module ${version} initialized from existing script`);
              resolve(module);
            })
            .catch((err) => {
              console.error(`Error initializing existing RInChI module ${version}:`, err);
              reject(err);
            });
        } else {
          reject(
            new Error(
              `Module function ${versionConfig.moduleName} not found after script was loaded`
            )
          );
        }
        return;
      }

      // Load the script
      const script = document.createElement("script");
      script.src = versionConfig.scriptSrc;
      script.id = `rinchi-script-${version}`;
      script.async = true;
      script.crossOrigin = "anonymous"; // Add CORS attribute

      script.onload = () => {
        // Ensure the module function exists
        const factory = window[versionConfig.moduleName];
        if (typeof factory === "function") {
          // Initialize the module
          factory()
            .then((module) => {
              moduleInstances[version] = module;
              resolve(module);
            })
            .catch((err) => {
              console.error(`Error initializing RInChI module ${version}:`, err);
              reject(err);
            });
        } else {
          const error = new Error(
            `Module function ${versionConfig.moduleName} not found in window object`
          );
          console.error(error);
          reject(error);
        }
      };

      script.onerror = (e) => {
        const error = new Error(
          `Failed to load RInChI module ${version} from ${versionConfig.scriptSrc}. Check network connection or CORS configuration.`
        );
        console.error(error, e);
        reject(error);
      };

      document.body.appendChild(script);
    } catch (err) {
      console.error(`Unexpected error loading RInChI module ${version}:`, err);
      reject(err);
    }
  });

  return moduleLoadPromises[version];
};

/**
 * Converts a RXN/RD file to RInChI
 * @param {string} rxnfile - RXN/RD file content
 * @param {boolean} forceEquilibrium - Whether to force equilibrium
 * @param {string} version - RInChI version to use
 * @returns {Promise<Object>} - RInChI result
 */
export const convertRxnfileToRinchi = async (
  rxnfile: string,
  forceEquilibrium: boolean,
  version: string
): Promise<RInChIResult> => {
  if (!rxnfile || !rxnfile.trim()) {
    return Promise.reject(new Error("No RXN/RD file provided"));
  }

  try {
    const module = await loadRinchiModule(version);

    const out_rinchi_stringPtr = module._malloc(4);
    const out_rinchi_auxinfoPtr = module._malloc(4);
    try {
      const res = module.ccall(
        "rinchilib_rinchi_from_file_text",
        "number",
        ["string", "string", "boolean", "number", "number"],
        ["AUTO", rxnfile, forceEquilibrium, out_rinchi_stringPtr, out_rinchi_auxinfoPtr]
      );
      const rinchi = module.UTF8ToString(module.getValue(out_rinchi_stringPtr, "i32"));
      const rauxinfo = module.UTF8ToString(module.getValue(out_rinchi_auxinfoPtr, "i32"));

      let error = "";
      if (res !== 0) {
        error = module.ccall("rinchilib_latest_err_msg", "string", [], []);
      }

      return {
        return_code: res,
        rinchi: rinchi,
        rauxinfo: rauxinfo,
        error: error,
      };
    } finally {
      module._free(out_rinchi_stringPtr);
      module._free(out_rinchi_auxinfoPtr);
    }
  } catch (err) {
    console.error("Error converting RXN/RD file to RInChI:", err);
    throw new Error(`Failed to convert reaction to RInChI: ${errorText(err)}`);
  }
};

/**
 * Generates a RInChIKey from a RInChI string
 * @param {string} rinchi - RInChI string
 * @param {string} keyType - Type of key to generate ("Long", "Short", or "Web")
 * @param {string} version - RInChI version to use
 * @returns {Promise<Object>} - RInChIKey result
 */
export const generateRinchiKey = async (
  rinchi: string,
  keyType: string,
  version: string
): Promise<RInChIKeyResult> => {
  if (!rinchi || !rinchi.trim() || !rinchi.startsWith("RInChI=")) {
    return Promise.reject(new Error("Invalid RInChI provided"));
  }

  if (!keyType || !["Long", "Short", "Web"].includes(keyType)) {
    return Promise.reject(
      new Error("Invalid key type provided. Must be 'Long', 'Short', or 'Web'")
    );
  }

  try {
    const module = await loadRinchiModule(version);

    const out_rinchi_keyPtr = module._malloc(4);
    const res = module.ccall(
      "rinchilib_rinchikey_from_rinchi",
      "number",
      ["string", "string", "number"],
      [rinchi, keyType, out_rinchi_keyPtr]
    );
    const rinchikey = module.UTF8ToString(module.getValue(out_rinchi_keyPtr, "i32"));
    module._free(out_rinchi_keyPtr);

    let error = "";
    if (res !== 0) {
      error = module.ccall("rinchilib_latest_err_msg", "string", [], []);
    }

    return {
      return_code: res,
      rinchikey: rinchikey,
      error: error,
    };
  } catch (err) {
    console.error(`Error generating ${keyType}-RInChIKey:`, err);
    throw new Error(`Failed to generate ${keyType}-RInChIKey: ${errorText(err)}`);
  }
};

/**
 * Converts a RInChI string to RXN/RD file
 * @param {string} rinchi - RInChI string
 * @param {string} rauxinfo - RAuxInfo string (optional)
 * @param {string} format - Output format ("RXN" or "RD")
 * @param {string} version - RInChI version to use
 * @returns {Promise<Object>} - File result
 */
export const convertRinchiToFileText = async (
  rinchi: string,
  rauxinfo: string,
  format: string,
  version: string
): Promise<RInChIFileResult> => {
  if (!rinchi || !rinchi.trim() || !rinchi.startsWith("RInChI=")) {
    return Promise.reject(new Error("Invalid RInChI provided"));
  }

  if (rauxinfo && !rauxinfo.startsWith("RAuxInfo=")) {
    return Promise.reject(new Error("Invalid RAuxInfo provided"));
  }

  if (!format || !["RXN", "RD"].includes(format)) {
    return Promise.reject(new Error("Invalid format provided. Must be 'RXN' or 'RD'"));
  }

  try {
    const module = await loadRinchiModule(version);

    const out_file_textPtr = module._malloc(4);
    const res = module.ccall(
      "rinchilib_file_text_from_rinchi",
      "number",
      ["string", "string", "string", "number"],
      [rinchi, rauxinfo || "", format, out_file_textPtr]
    );
    const fileText = module.UTF8ToString(module.getValue(out_file_textPtr, "i32"));
    module._free(out_file_textPtr);

    let error = "";
    if (res !== 0) {
      error = module.ccall("rinchilib_latest_err_msg", "string", [], []);
    }

    return {
      return_code: res,
      fileText: fileText,
      error: error,
    };
  } catch (err) {
    console.error("Error converting RInChI to file text:", err);
    throw new Error(`Failed to convert RInChI to ${format} file: ${errorText(err)}`);
  }
};
