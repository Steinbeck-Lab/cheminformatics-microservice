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
  "1.07.3": {
    label: "InChI 1.07.3",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web107.js`,
    moduleName: "inchiModule107",
    default: true,
  },
  "1.07.3-orgmet": {
    label: "InChI 1.07.3 with Molecular inorganics",
    scriptSrc: `${INCHI_BASE_URL}/inchi/inchi-web107-orgmet.js`,
    moduleName: "inchiModule107OrgMet",
  },
};

// Module references
const moduleInstances = {};
const moduleLoadPromises = {};

/**
 * Loads an InChI module from remote URL
 * @param {string} version - Version identifier ("1.06", "1.07.3", or "1.07.3-orgmet")
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

// Set to true for development without WASM modules,
// set to false to use the remote WASM modules
export const useMockImplementation = false;

/**
 * Mock implementation of convertMolfileToInchi
 * @param {string} molfile - Molfile content
 * @param {string} options - InChI options string
 * @returns {Promise<Object>} - Mock InChI result
 */
export const mockConvertMolfileToInchi = async (molfile, options) => {
  await new Promise((resolve) => setTimeout(resolve, 500)); // Simulate processing time

  // Simple heuristic to categorize the structure based on molfile content
  const containsBenzene =
    molfile.includes("C6H6") ||
    (molfile.match(/C\s+0\s+0/g)?.length >= 6 &&
      molfile.match(/H\s+0\s+0/g)?.length >= 6);

  const containsEthanol =
    molfile.includes("C2H6O") ||
    (molfile.match(/C\s+0\s+0/g)?.length >= 2 &&
      molfile.match(/O\s+0\s+0/g)?.length >= 1 &&
      molfile.match(/H\s+0\s+0/g)?.length >= 6);

  if (containsBenzene) {
    return {
      return_code: 0,
      inchi: "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
      auxinfo:
        "AuxInfo=1/0/N:1,2,3,4,5,6/E:(1,2,3,4,5,6)/rA:6nCCCCCC/rB:d1;d2;d3;d4;d5;/rC:;;;;;",
      message: `Generated with options: ${options}`,
      log: "",
    };
  } else if (containsEthanol) {
    return {
      return_code: 0,
      inchi: "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
      auxinfo: "AuxInfo=1/0/N:1,2,3/E:(1,2)/rA:3nCCO/rB:s1;s2;/rC:;;;",
      message: `Generated with options: ${options}`,
      log: "",
    };
  } else {
    // Generic mock response
    return {
      return_code: 0,
      inchi: "InChI=1S/C10H16/c1-7-4-5-8(2)10(7)6-3-9(10)1/h7-9H,3-6H2,1-2H3",
      auxinfo:
        "AuxInfo=1/0/N:1,2,3,4,5,6,7,8,9,10/rA:10nCCCCCCCCCC/rB:s1;s2;s3;s4;s5;s6;s7;s8;s9;/rC:;;;;;;;;;;",
      message: `Generated with options: ${options}`,
      log: "",
    };
  }
};

/**
 * Mock implementation of generateInchiKey
 * @param {string} inchi - InChI string
 * @returns {Promise<Object>} - Mock InChIKey result
 */
export const mockGenerateInchiKey = async (inchi) => {
  await new Promise((resolve) => setTimeout(resolve, 300)); // Simulate processing time

  // Map of known InChI strings to their InChIKeys
  const knownKeys = {
    "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H": "UHOVQNZJYSORNB-UHFFFAOYSA-N",
    "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)":
      "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3":
      "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
    "InChI=1S/C10H16/c1-7-4-5-8(2)10(7)6-3-9(10)1/h7-9H,3-6H2,1-2H3":
      "YQIGLHPSVFAIAV-UHFFFAOYSA-N",
  };

  // Return the known key or a generated one based on the hash of the InChI string
  if (knownKeys[inchi]) {
    return {
      return_code: 0,
      inchikey: knownKeys[inchi],
      message: "",
    };
  } else {
    // Generate a simple hash value from the InChI string
    let hash = 0;
    for (let i = 0; i < inchi.length; i++) {
      hash = (hash << 5) - hash + inchi.charCodeAt(i);
      hash = hash & hash; // Convert to 32bit integer
    }

    // Generate a mock InChIKey
    const mockKey = `MOCK${Math.abs(hash)
      .toString(16)
      .padStart(10, "0")
      .toUpperCase()}-UHFFFAOYSA-N`;

    return {
      return_code: 0,
      inchikey: mockKey,
      message: "",
    };
  }
};

/**
 * Mock implementation of convertInchiToMolfile
 * @param {string} inchi - InChI string
 * @returns {Promise<Object>} - Mock molfile result
 */
export const mockConvertInchiToMolfile = async (inchi) => {
  await new Promise((resolve) => setTimeout(resolve, 500)); // Simulate processing time

  let mockMolfile = "";

  // Generate different mock molfiles based on the input InChI
  if (inchi === "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H") {
    // Benzene
    mockMolfile = `
Benzene
  
  6  6  0  0  0  0  0  0  0  0999 V2000
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
`;
  } else if (inchi === "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3") {
    // Ethanol
    mockMolfile = `
Ethanol
  
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
`;
  } else {
    // Generic structure for other InChI strings
    mockMolfile = `
Structure
  
  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
M  END
`;
  }

  return {
    return_code: 0,
    molfile: mockMolfile,
    message: "",
    log: "",
  };
};
