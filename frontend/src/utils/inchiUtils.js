/**
 * inchiUtils.js - Utility functions for InChI module loading and operations
 *
 * This file provides utilities for loading InChI WASM modules from remote URLs and performing conversions
 * like structure to InChI, InChI to structure, and generating InChIKeys.
 */

// Base URL for the GitHub Pages-hosted WASM files
const INCHI_BASE_URL = "https://iupac-inchi.github.io/InChI-Web-Demo";

// Configuration for available InChI modules with remote URLs
export const INCHI_VERSIONS = {
  1.06: {
    label: "InChI 1.06",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web106.js`,
    moduleName: "inchiModule106",
  },
  "Latest": {
    label: "InChI Latest",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-latest.js`,
    moduleName: "inchiModuleLatest",
    default: true,
  },
  "Latest-MoIn": {
    label: "InChI latest with Molecular inorganics",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web-latest-moin.js`,
    moduleName: "inchiModuleLatestMoIn",
  },
};

// Module references
const moduleInstances = {};
const moduleLoadPromises = {};

/**
 * Loads an InChI module from remote URL
 * @param {string} version - Version identifier ("1.06", "Latest", or "Latest-MoIn")
 * @returns {Promise<Object>} - Promise that resolves to the module instance
 */
export const loadInchiModule = async (version) => {
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
        if (typeof window[versionConfig.moduleName] === "function") {
          window[versionConfig.moduleName]()
            .then((module) => {
              moduleInstances[version] = module;
              console.log(
                `InChI module ${version} initialized from existing script`
              );
              resolve(module);
            })
            .catch((err) => {
              console.error(
                `Error initializing existing InChI module ${version}:`,
                err
              );
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
        if (typeof window[versionConfig.moduleName] === "function") {
          // Initialize the module
          window[versionConfig.moduleName]()
            .then((module) => {
              moduleInstances[version] = module;
              console.log(
                `InChI module ${version} loaded successfully from ${versionConfig.scriptSrc}`
              );
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
export const convertMolfileToInchi = async (molfile, options, version) => {
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
    throw new Error(`Failed to convert structure to InChI: ${err.message}`);
  }
};

/**
 * Generates an InChIKey from an InChI string
 * @param {string} inchi - InChI string
 * @param {string} version - InChI version to use
 * @returns {Promise<Object>} - InChIKey result
 */
export const generateInchiKey = async (inchi, version) => {
  if (!inchi || !inchi.trim() || !inchi.startsWith("InChI=")) {
    return Promise.reject(new Error("Invalid InChI provided"));
  }

  try {
    const module = await loadInchiModule(version);

    const ptr = module.ccall(
      "inchikey_from_inchi",
      "number",
      ["string"],
      [inchi]
    );

    const resultStr = module.UTF8ToString(ptr);
    module._free(ptr);

    return JSON.parse(resultStr);
  } catch (err) {
    console.error("Error generating InChIKey:", err);
    throw new Error(`Failed to generate InChIKey: ${err.message}`);
  }
};

/**
 * Converts an InChI string to a molfile
 * @param {string} inchi - InChI string
 * @param {string} options - InChI options string (optional)
 * @param {string} version - InChI version to use
 * @returns {Promise<Object>} - Molfile result
 */
export const convertInchiToMolfile = async (inchi, options, version) => {
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
    throw new Error(`Failed to convert InChI to structure: ${err.message}`);
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
  auxinfo,
  doNotAddH = 0,
  diffUnkUndfStereo = 0,
  version
) => {
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
    throw new Error(`Failed to convert AuxInfo to structure: ${err.message}`);
  }
};

/**
 * Converts a SMILES string to molfile format
 * This is a helper function to be implemented if needed
 * Not directly supported by InChI WASM modules
 */
export const convertSmilesToMolfile = async (smiles) => {
  // This would require integration with a SMILES parser or service
  // Not implemented in the original InChI WASM modules
  throw new Error(
    "SMILES to Molfile conversion not directly supported by InChI modules"
  );
};
