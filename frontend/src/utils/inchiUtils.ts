/**
 * inchiUtils.js - Utility functions for InChI module loading and operations
 *
 * This file provides utilities for loading InChI WASM modules from remote URLs and performing conversions
 * like structure to InChI, InChI to structure, and generating InChIKeys.
 */

import { errorText } from "@/lib/error-messages";
import type { InChIKeyResult, InChIResult, InChIVersionConfig } from "@/types/api";

// Base URL for the GitHub Pages-hosted WASM files
const INCHI_BASE_URL = "https://iupac-inchi.github.io/InChI-Web-Demo";

// Configuration for available InChI modules with remote URLs
export const INCHI_VERSIONS: Record<string, InChIVersionConfig> = {
  1.06: {
    label: "InChI 1.06",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web106.js`,
    moduleName: "inchiModule106",
  },
  // File and module names follow pages/inchi_versions.json in IUPAC-InChI/InChI-Web-Demo.
  Latest: {
    label: "InChI 1.07.5",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-1075.js`,
    moduleName: "inchiModule1075",
    default: true,
  },
  "Latest-MoIn": {
    label: "InChI dev with Molecular Inorganics",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-dev-moin.js`,
    moduleName: "inchiModuleDevMoIn",
  },
  "Dev-EnhancedStereo": {
    label: "InChI dev with Enhanced Stereochemistry",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-dev-enhanced-stereo.js`,
    moduleName: "inchiModuleDevEnhancedStereo",
  },
  Polymer: {
    label: "InChI with Polymer Support",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-polymer-support.js`,
    moduleName: "inchiModulePolymerSupport",
  },
  Dev: {
    label: "InChI dev",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-dev.js`,
    moduleName: "inchiModuleDev",
  },
  "Dev-NoMetalH": {
    label: "InChI dev, no metal hydrogens",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-no-metal-h.js`,
    moduleName: "inchiModuleNoMetalH",
  },
  "Dev-ExplicitZeroValence": {
    label: "InChI dev, explicit zero valence",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-explicit-zero-valence.js`,
    moduleName: "inchiModuleExplicitZeroValence",
  },
};

/** [version, label] pairs for a version picker. */
export const VERSION_OPTIONS: [string, string][] = Object.entries(INCHI_VERSIONS).map(
  ([value, v]) => [value, v.label]
);

/** The version marked `default: true`. */
export const DEFAULT_VERSION =
  Object.keys(INCHI_VERSIONS).find((key) => INCHI_VERSIONS[key]?.default) || "Latest";

// Module references
const moduleInstances: Record<string, InchiWasmModule> = {};
const moduleLoadPromises: Record<string, Promise<InchiWasmModule>> = {};

/**
 * Loads an InChI module from remote URL
 * @param {string} version - Version identifier ("1.06", "Latest", or "Latest-MoIn")
 * @returns {Promise<Object>} - Promise that resolves to the module instance
 */
export const loadInchiModule = async (version: string): Promise<InchiWasmModule> => {
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
      const versionConfig = INCHI_VERSIONS[version];

      if (!versionConfig) {
        reject(new Error(`Unknown InChI version: ${version}`));
        return;
      }

      // Check if the script is already loaded
      const existingScript = document.getElementById(`inchi-script-${version}`);
      if (existingScript) {
        // Script exists, try to initialize the module
        const factory = window[versionConfig.moduleName];
        if (typeof factory === "function") {
          factory()
            .then((module) => {
              moduleInstances[version] = module;
              console.log(`InChI module ${version} initialized from existing script`);
              resolve(module);
            })
            .catch((err) => {
              console.error(`Error initializing existing InChI module ${version}:`, err);
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
      script.id = `inchi-script-${version}`;
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
              console.error(`Error initializing InChI module ${version}:`, err);
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
          `Failed to load InChI module ${version} from ${versionConfig.scriptSrc}. Check network connection or CORS configuration.`
        );
        console.error(error, e);
        reject(error);
      };

      document.body.appendChild(script);
    } catch (err) {
      console.error(`Unexpected error loading InChI module ${version}:`, err);
      reject(err);
    }
  });

  return moduleLoadPromises[version];
};

/**
 * Converts a molfile to InChI
 * @param {string} molfile - Molfile content
 * @param {string} options - InChI options string
 * @param {string} version - InChI version to use
 * @returns {Promise<Object>} - InChI result
 */
export const convertMolfileToInchi = async (
  molfile: string,
  options: string,
  version: string
): Promise<InChIResult> => {
  if (!molfile || !molfile.trim()) {
    return Promise.reject(new Error("No molfile provided"));
  }

  try {
    const module = await loadInchiModule(version);

    const ptr = module.ccall(
      "inchi_from_molfile",
      "number",
      ["string", "string"],
      [molfile, options]
    );

    const resultStr = module.UTF8ToString(ptr);
    module._free(ptr);

    return JSON.parse(resultStr);
  } catch (err) {
    console.error("Error converting molfile to InChI:", err);
    throw new Error(`Failed to convert structure to InChI: ${errorText(err)}`);
  }
};

/**
 * Generates an InChIKey from an InChI string
 * @param {string} inchi - InChI string
 * @param {string} version - InChI version to use
 * @returns {Promise<Object>} - InChIKey result
 */
export const generateInchiKey = async (inchi: string, version: string): Promise<InChIKeyResult> => {
  if (!inchi || !inchi.trim() || !inchi.startsWith("InChI=")) {
    return Promise.reject(new Error("Invalid InChI provided"));
  }

  try {
    const module = await loadInchiModule(version);

    const ptr = module.ccall("inchikey_from_inchi", "number", ["string"], [inchi]);

    const resultStr = module.UTF8ToString(ptr);
    module._free(ptr);

    return JSON.parse(resultStr);
  } catch (err) {
    console.error("Error generating InChIKey:", err);
    throw new Error(`Failed to generate InChIKey: ${errorText(err)}`);
  }
};

/**
 * Converts an InChI string to a molfile
 * @param {string} inchi - InChI string
 * @param {string} options - InChI options string (optional)
 * @param {string} version - InChI version to use
 * @returns {Promise<Object>} - Molfile result
 */
export const convertInchiToMolfile = async (
  inchi: string,
  options: string,
  version: string
): Promise<InChIResult> => {
  if (!inchi || !inchi.trim() || !inchi.startsWith("InChI=")) {
    return Promise.reject(new Error("Invalid InChI provided"));
  }

  try {
    const module = await loadInchiModule(version);

    const ptr = module.ccall(
      "molfile_from_inchi",
      "number",
      ["string", "string"],
      [inchi, options || ""]
    );

    const resultStr = module.UTF8ToString(ptr);
    module._free(ptr);

    return JSON.parse(resultStr);
  } catch (err) {
    console.error("Error converting InChI to molfile:", err);
    throw new Error(`Failed to convert InChI to structure: ${errorText(err)}`);
  }
};

/**
 * Converts AuxInfo to a molfile
 * @param {string} auxinfo - AuxInfo string
 * @param {number} doNotAddH - Option to control hydrogen addition (0 or 1)
 * @param {number} diffUnkUndfStereo - Option to differentiate unknown/undefined stereo (0 or 1)
 * @param {string} version - InChI version to use
 * @returns {Promise<Object>} - Molfile result
 */
export const convertAuxinfoToMolfile = async (
  auxinfo: string,
  doNotAddH = 0,
  diffUnkUndfStereo = 0,
  version: string
): Promise<InChIResult> => {
  if (!auxinfo || !auxinfo.trim() || !auxinfo.startsWith("AuxInfo=")) {
    return Promise.reject(new Error("Invalid AuxInfo provided"));
  }

  try {
    const module = await loadInchiModule(version);

    const ptr = module.ccall(
      "molfile_from_auxinfo",
      "number",
      ["string", "number", "number"],
      [auxinfo, doNotAddH, diffUnkUndfStereo]
    );

    const resultStr = module.UTF8ToString(ptr);
    module._free(ptr);

    return JSON.parse(resultStr);
  } catch (err) {
    console.error("Error converting AuxInfo to molfile:", err);
    throw new Error(`Failed to convert AuxInfo to structure: ${errorText(err)}`);
  }
};
