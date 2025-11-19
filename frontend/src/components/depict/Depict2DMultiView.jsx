// Description: This component allows users to generate 2D depictions of multiple molecules from SMILES strings.
import React, { useState } from "react"; // Removed useEffect
import {
  HiOutlinePhotograph,
  HiOutlineDownload,
  HiOutlineClipboard,
  HiOutlineCheck,
  HiOutlineSwitchHorizontal,
  HiOutlineExclamationCircle,
} from "react-icons/hi";
import { motion, AnimatePresence, LayoutGroup } from "framer-motion";
// Assuming these components are correctly implemented and styled for dark/light mode
import LoadingScreen from "../common/LoadingScreen";
// Assuming this service is configured correctly
import depictService from "../../services/depictService"; // Assuming this service exists
import api from "../../services/api";

// Animation variants
const resultsContainerVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { staggerChildren: 0.05 } }, // Stagger children animation
};

const depictionCardVariant = {
  hidden: { opacity: 0, y: 20, scale: 0.95 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { 
      duration: 0.5, 
      ease: [0.25, 0.46, 0.45, 0.94],
      scale: {
        type: "spring",
        stiffness: 300,
        damping: 25
      }
    },
  },
};

const toggleSpring = {
  type: "spring",
  stiffness: 600,
  damping: 30,
};

// Toggle Switch Component - Enhanced with stunning visuals
const ToggleSwitch = ({ id, checked, onChange, label, disabled = false }) => {
  return (
    <div className="flex items-center space-x-3 py-1 px-0.5">
      <div
        className={`relative flex items-center w-12 h-6 rounded-full p-0.5 cursor-pointer transition-all duration-300 ease-in-out ${
          checked
            ? "bg-gradient-to-r from-blue-600 via-blue-500 to-indigo-500 border-2 border-cyan-400 dark:border-cyan-300"
            : "bg-gradient-to-r from-gray-300 via-gray-350 to-gray-400 dark:from-gray-600 dark:via-gray-650 dark:to-gray-700"
        } ${disabled ? "opacity-50 cursor-not-allowed" : "hover:scale-105"}`}
        onClick={() => !disabled && onChange(!checked)}
        role="switch"
        aria-checked={checked}
        aria-label={label}
      >        
        <LayoutGroup>
          <motion.div
            className={`relative h-5 w-5 rounded-full ${
              checked
                ? "bg-gradient-to-br from-white via-blue-50 to-indigo-50"
                : "bg-gradient-to-br from-white via-gray-50 to-gray-100 dark:from-gray-200 dark:via-gray-300 dark:to-gray-400"
            }`}
            layout
            transition={toggleSpring}
            style={{
              left: checked ? "auto" : "2px",
              right: checked ? "2px" : "auto",
              position: "absolute",
            }}
          />
        </LayoutGroup>
      </div>
      <label
        htmlFor={id}
        className={`text-sm font-medium select-none transition-colors duration-200 ${
          disabled 
            ? "text-gray-400 dark:text-gray-500" 
            : checked
              ? "text-blue-700 dark:text-blue-400 cursor-pointer"
              : "text-gray-700 dark:text-gray-300 cursor-pointer hover:text-gray-900 dark:hover:text-gray-100"
        }`}
        onClick={() => !disabled && onChange(!checked)}
      >
        {label}
      </label>
    </div>
  );
};

const BatchDepictionView = () => {
  // Input state
  const [inputText, setInputText] = useState("");

  // Processing state
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [depictions, setDepictions] = useState([]); // Array of { smiles, title, imageUrl, id }

  // Toolkit selection
  const [toolkit, setToolkit] = useState("rdkit");

  // Depiction options - Basic (shared by both toolkits)
  const [width, setWidth] = useState(300);
  const [height, setHeight] = useState(200);
  const [showCIP, setShowCIP] = useState(false);
  const [useUnicolor, setUseUnicolor] = useState(false);
  const [highlight, setHighlight] = useState("");
  const [showAtomNumbers, setShowAtomNumbers] = useState(false);
  const [hydrogenDisplay, setHydrogenDisplay] = useState("Smart");

  // Depiction options - Enhanced features (CDK only)
  const [abbreviate, setAbbreviate] = useState("off");
  const [dative, setDative] = useState("metals");
  const [multicenter, setMulticenter] = useState("provided");
  const [style, setStyle] = useState("cow");
  const [annotate, setAnnotate] = useState("none");
  const [donuts, setDonuts] = useState(false);
  const [zoom, setZoom] = useState(1.0);
  const [ratio, setRatio] = useState(1.0);
  const [flip, setFlip] = useState(false);
  const [showtitle, setShowtitle] = useState(false);
  const [perceiveRadicals, setPerceiveRadicals] = useState(false);

  // UI state
  const [copiedSmiles, setCopiedSmiles] = useState(false); // For "Copy All SMILES" button
  const [copiedSingle, setCopiedSingle] = useState(null); // Tracks ID of molecule whose SMILES was copied
  const [downloadFormat, setDownloadFormat] = useState("svg"); // svg or png
  const [showToolsSection, setShowToolsSection] = useState(false); // Toggle visibility of options

  // Fixed radical results per depiction id -> { fixed_smiles, radicals_detected, radicals_fixed, fixedImageUrl }
  const [fixedResults, setFixedResults] = useState({});
  const [fetchingFixed, setFetchingFixed] = useState(false);

  // Per-molecule rotation state { [id]: rotationValue }
  const [rotations, setRotations] = useState({});

  // Parse SMILES (and optional titles) from input text
  const parseSmiles = (text) => {
    return text
      .split(/[\n\r]+/) // Split by new lines
      .map((line) => line.trim()) // Trim whitespace
      .filter((line) => line.length > 0 && !line.startsWith("#")); // Remove empty lines and comments
  };

  // Fetch fixed SMILES for each depiction by calling backend /chem/fixRadicals
  const fetchFixedForDepictions = async (depictionList) => {
    if (!Array.isArray(depictionList) || depictionList.length === 0) return;
    setFetchingFixed(true);
    const mapping = {};

    for (const dep of depictionList) {
      try {
        const resp = await api.get(`/chem/fixRadicals`, { params: { smiles: dep.smiles } });
        const data = resp.data || {};

        // Build a fixed image URL for the fixed smiles (use RDKit basic depiction)
        let fixedImageUrl = null;
        if (data.fixed_smiles) {
          try {
            fixedImageUrl = depictService.get2DDepictionUrl(data.fixed_smiles, {
              toolkit: "rdkit",
              width,
              height,
              rotate: 0,
              showAtomNumbers,
            });
          } catch (err) {
            console.warn("Failed to build fixed image url:", err);
          }
        }

        mapping[dep.id] = {
          fixed_smiles: data.fixed_smiles || null,
          radicals_detected: data.radicals_detected || 0,
          radicals_fixed: data.radicals_fixed || 0,
          fixedImageUrl,
        };
      } catch (err) {
        console.warn(`Error fetching fixed radicals for ${dep.smiles}:`, err.message || err);
        mapping[dep.id] = {
          fixed_smiles: null,
          radicals_detected: 0,
          radicals_fixed: 0,
          fixedImageUrl: null,
        };
      }
    }

    setFixedResults(mapping);
    setFetchingFixed(false);
  };

  // Regenerate depiction URLs when options change (or called manually)
  const regenerateDepictions = () => {
    if (depictions.length === 0) return;

    console.log(`Regenerating depictions with toolkit: ${toolkit}`);

    setDepictions((currentDepictions) =>
      currentDepictions.map((dep) => {
        const rotation = rotations[dep.id] || 0;
        
        // Use appropriate endpoint based on toolkit
        if (toolkit === "cdk") {
          // CDK uses enhanced endpoint with all features
          if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
            console.error("depictService.get2DDepictionUrlEnhanced is not available.");
            return dep;
          }
          
          const options = {
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            hydrogen_display: hydrogenDisplay,
            abbreviate,
            dative,
            multicenter,
            style,
            annotate,
            donuts,
            zoom,
            ratio,
            flip,
            showtitle,
            perceive_radicals: perceiveRadicals,
          };
          const updatedImageUrl = depictService.get2DDepictionUrlEnhanced(dep.smiles, options);
          return { ...dep, imageUrl: updatedImageUrl };
        } else {
          // RDKit uses basic endpoint
          if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
            console.error("depictService.get2DDepictionUrl is not available.");
            return dep;
          }
          
          const options = {
            toolkit: "rdkit",
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
          };
          const updatedImageUrl = depictService.get2DDepictionUrl(dep.smiles, options);
          return { ...dep, imageUrl: updatedImageUrl };
        }
      })
    );
  };

  // Handle toolkit change and regenerate
  const handleToolkitChange = (newToolkit) => {
    setToolkit(newToolkit);
    // Regenerate depictions if results already exist
    if (depictions.length > 0) {
      // We'll regenerate after state updates
      setTimeout(() => {
        regenerateDepictions();
      }, 0);
    }
  };

  // Switch between available toolkits
  const handleSwitchToolkit = () => {
    handleToolkitChange(toolkit === "rdkit" ? "cdk" : "rdkit");
  };

  // Handle text input change, toggle tools visibility
  const handleInputChange = (value) => {
    setInputText(value);
    // Show tools section automatically if there's input
    if (value.trim() && !showToolsSection) {
      setShowToolsSection(true);
    } else if (!value.trim() && showToolsSection) {
      // Keep tools section visible even if input is cleared for easier option changes
      // setShowToolsSection(false);
    }
  };

  // --- Removed useEffect hook for automatic regeneration on option changes ---

  // Handle rotation change for a specific molecule
  const handleRotationChange = (id, value) => {
    const rotation = Number(value);
    // Update rotation state for the specific molecule
    setRotations((prev) => ({ ...prev, [id]: rotation }));

    // Update the image URL immediately for this molecule
    setDepictions((currentDepictions) =>
      currentDepictions.map((dep) => {
        if (dep.id === id) {
          // Use appropriate endpoint based on toolkit
          if (toolkit === "cdk") {
            if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
              console.error("depictService.get2DDepictionUrlEnhanced is not available.");
              return dep;
            }
            
            const options = {
              width,
              height,
              rotate: rotation,
              CIP: showCIP,
              unicolor: useUnicolor,
              highlight: highlight || undefined,
              showAtomNumbers,
              hydrogen_display: hydrogenDisplay,
              abbreviate,
              dative,
              multicenter,
              style,
              annotate,
              donuts,
              zoom,
              ratio,
              flip,
              showtitle,
              perceive_radicals: perceiveRadicals,
            };
            const updatedImageUrl = depictService.get2DDepictionUrlEnhanced(dep.smiles, options);
            return { ...dep, imageUrl: updatedImageUrl };
          } else {
            if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
              console.error("depictService.get2DDepictionUrl is not available.");
              return dep;
            }
            
            const options = {
              toolkit: "rdkit",
              width,
              height,
              rotate: rotation,
              CIP: showCIP,
              unicolor: useUnicolor,
              highlight: highlight || undefined,
              showAtomNumbers,
            };
            const updatedImageUrl = depictService.get2DDepictionUrl(dep.smiles, options);
            return { ...dep, imageUrl: updatedImageUrl };
          }
        }
        return dep;
      })
    );
  };

  // Handle form submission to generate initial depictions
  const handleSubmit = async (e) => {
    e.preventDefault();

    const smilesList = parseSmiles(inputText);
    if (smilesList.length === 0) {
      setError("Please enter at least one valid SMILES string.");
      setDepictions([]);
      return;
    }

    setLoading(true);
    setError(null);
    setDepictions([]);
    setRotations({}); // Reset rotations

    try {
      const results = [];
      const maxToProcess = 50; // Limit processing
      const processingList = smilesList.slice(0, maxToProcess);

      if (smilesList.length > maxToProcess) {
        setError(`Processing limited to the first ${maxToProcess} molecules.`);
      }

      const initialRotations = {};
      for (let i = 0; i < processingList.length; i++) {
        const smilesLine = processingList[i];
        // Extract SMILES (first part before space/tab)
        const spaceIndex = smilesLine.search(/[\s\t]/);
        const smiles = spaceIndex > 0 ? smilesLine.substring(0, spaceIndex) : smilesLine;

        // Basic SMILES validation (can be improved)
        if (!smiles || smiles.length < 1) {
          console.warn(`Skipping invalid SMILES line: ${smilesLine}`);
          continue;
        }

        // Extract Title (rest of the line after SMILES)
        const title =
          spaceIndex > 0 ? smilesLine.substring(spaceIndex + 1).trim() : `Molecule ${i + 1}`;
        const id = `mol-${i}-${Date.now()}`; // Unique ID for key and state

        initialRotations[id] = 0; // Set initial rotation

        let imageUrl;
        
        // Use appropriate endpoint based on toolkit
        if (toolkit === "cdk") {
          if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
            console.error("depictService.get2DDepictionUrlEnhanced is not available.");
            setError("Depiction service is not configured correctly.");
            return;
          }
          
          const options = {
            width,
            height,
            rotate: 0,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            hydrogen_display: hydrogenDisplay,
            abbreviate,
            dative,
            multicenter,
            style,
            annotate,
            donuts,
            zoom,
            ratio,
            flip,
            showtitle,
            perceive_radicals: perceiveRadicals,
          };
          imageUrl = depictService.get2DDepictionUrlEnhanced(smiles, options);
        } else {
          if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
            console.error("depictService.get2DDepictionUrl is not available.");
            setError("Depiction service is not configured correctly.");
            return;
          }
          
          const options = {
            toolkit: "rdkit",
            width,
            height,
            rotate: 0,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
          };
          imageUrl = depictService.get2DDepictionUrl(smiles, options);
        }

        results.push({ smiles, title, imageUrl, id });
      }

      setRotations(initialRotations); // Set all initial rotations at once
      setDepictions(results); // Set all results at once

      // Fetch fixed SMILES for depictions (if radicals present)
      try {
        await fetchFixedForDepictions(results);
      } catch (e) {
        console.warn("Failed to fetch fixed radicals:", e);
      }

      if (results.length === 0 && !error) {
        setError("No valid SMILES found to process.");
      }
    } catch (err) {
      console.error("Error generating depictions:", err);
      setError(`Error generating depictions: ${err.message || "Unknown error"}`);
      setDepictions([]);
      setRotations({});
    } finally {
      setLoading(false);
    }
  };

  // Download all depictions as a ZIP file
  // NOTE: This requires the 'jszip' library to be installed (npm install jszip)
  const downloadAllDepictions = async () => {
    if (depictions.length === 0) return;

    setLoading(true); // Show loading indicator during zip creation
    setError(null);

    try {
      // Dynamically import JSZip only when needed
      const JSZip = (await import("jszip")).default;
      const zip = new JSZip();

      // Fetch all depiction blobs concurrently
      const fetchPromises = depictions.map(async (depiction) => {
        const safeTitle =
          depiction.title.replace(/[^a-z0-9]/gi, "_").toLowerCase() || `molecule_${depiction.id}`;
        const filename = `${safeTitle}.${downloadFormat}`; // Use selected format
        const rotation = rotations[depiction.id] || 0; // Get current rotation

        let downloadUrl;
        
        // Use appropriate endpoint based on toolkit
        if (toolkit === "cdk") {
          if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
            throw new Error("Depiction service (enhanced) is not configured correctly.");
          }
          
          const options = {
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            format: downloadFormat,
            hydrogen_display: hydrogenDisplay,
            abbreviate,
            dative,
            multicenter,
            style,
            annotate,
            donuts,
            zoom,
            ratio,
            flip,
            showtitle,
            perceive_radicals: perceiveRadicals,
          };
          downloadUrl = depictService.get2DDepictionUrlEnhanced(depiction.smiles, options);
        } else {
          if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
            throw new Error("Depiction service is not configured correctly.");
          }
          
          const options = {
            toolkit: "rdkit",
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            format: downloadFormat,
          };
          downloadUrl = depictService.get2DDepictionUrl(depiction.smiles, options);
        }

        try {
          const response = await fetch(downloadUrl);
          if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);
          const blob = await response.blob();
          zip.file(filename, blob); // Add file to zip
        } catch (fetchErr) {
          console.error(`Failed to fetch or add ${filename}:`, fetchErr);
          // Add an error file to the zip instead
          zip.file(
            `${filename}.error.txt`,
            `Failed to download image for ${depiction.smiles}. Error: ${fetchErr.message}`
          );
        }
      });

      await Promise.all(fetchPromises); // Wait for all fetches

      // Generate zip file content
      const content = await zip.generateAsync({ type: "blob" });
      const downloadUrl = URL.createObjectURL(content);

      // Trigger download
      const a = document.createElement("a");
      a.href = downloadUrl;
      a.download = `depictions_${toolkit}_${downloadFormat}.zip`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(downloadUrl); // Clean up
    } catch (err) {
      console.error("Error creating zip file:", err);
      setError(`Error downloading depictions: ${err.message || "Failed to create ZIP."}`);
    } finally {
      setLoading(false); // Hide loading indicator
    }
  };

  // Copy all SMILES strings to clipboard
  const copyAllSmiles = () => {
    if (depictions.length === 0 || !navigator.clipboard) return;
    const smilesText = depictions.map((d) => d.smiles).join("\n");
    navigator.clipboard
      .writeText(smilesText)
      .then(() => {
        setCopiedSmiles(true);
        setTimeout(() => setCopiedSmiles(false), 2000);
      })
      .catch((err) => {
        console.error("Failed to copy SMILES:", err);
        setError("Failed to copy SMILES to clipboard.");
      });
  };

  // Copy a single SMILES string
  const copySingleSmiles = (smiles, id) => {
    if (!navigator.clipboard) return;
    navigator.clipboard
      .writeText(smiles)
      .then(() => {
        setCopiedSingle(id);
        setTimeout(() => setCopiedSingle(null), 1500);
      })
      .catch((err) => {
        console.error("Failed to copy SMILES:", err);
        setError("Failed to copy SMILES.");
      });
  };

  // Download a single depiction image
  const downloadSingleDepiction = async (depiction) => {
    try {
      setLoading(true); // Show loading indicator while downloading
      const rotation = rotations[depiction.id] || 0;

      let url;
      
      // Use appropriate endpoint based on toolkit
      if (toolkit === "cdk") {
        if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
          console.error("depictService.get2DDepictionUrlEnhanced is not available.");
          setError("Depiction service is not configured correctly.");
          return;
        }
        
        const options = {
          width,
          height,
          rotate: rotation,
          CIP: showCIP,
          unicolor: useUnicolor,
          highlight: highlight || undefined,
          showAtomNumbers,
          format: downloadFormat,
          hydrogen_display: hydrogenDisplay,
          abbreviate,
          dative,
          multicenter,
          style,
          annotate,
          donuts,
          zoom,
          ratio,
          flip,
          showtitle,
          perceive_radicals: perceiveRadicals,
        };
        url = depictService.get2DDepictionUrlEnhanced(depiction.smiles, options);
      } else {
        if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
          console.error("depictService.get2DDepictionUrl is not available.");
          setError("Depiction service is not configured correctly.");
          return;
        }
        
        const options = {
          toolkit: "rdkit",
          width,
          height,
          rotate: rotation,
          CIP: showCIP,
          unicolor: useUnicolor,
          highlight: highlight || undefined,
          showAtomNumbers,
          format: downloadFormat,
        };
        url = depictService.get2DDepictionUrl(depiction.smiles, options);
      }

      // Fetch the image data as a blob
      const response = await fetch(url);
      if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);

      const blob = await response.blob();
      const downloadUrl = URL.createObjectURL(blob);

      // Trigger download using an anchor tag
      const a = document.createElement("a");
      a.href = downloadUrl;
      a.download = `${depiction.title.replace(/[^a-z0-9]/gi, "_").toLowerCase() || "molecule"}.${downloadFormat}`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);

      // Clean up the blob URL
      URL.revokeObjectURL(downloadUrl);
    } catch (err) {
      console.error("Error downloading single depiction:", err);
      setError(`Error downloading depiction: ${err.message || "Unknown error"}`);
    } finally {
      setLoading(false); // Hide loading indicator
    }
  };

  // Generate example molecules and populate input
  const generateExamples = () => {
    const examples = [
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C Caffeine",
      "CC(=O)OC1=CC=CC=C1C(=O)O Aspirin",
      "CC(=O)NC1=CC=C(C=C1)O Paracetamol",
      "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O Ibuprofen",
      "COc1cc2c(cc1OC)NC(=O)C3=C2C=CN3C Emetine",
      "CC(C)(C)NCC(O)c1ccc(O)c(O)c1 Salbutamol",
      "C1=CC=C(C=C1)C(=O)C(=O)O Phenylglyoxylic acid",
      "C1=CC=C(C=C1)C2=CC=C(C=C2)C(=O)O Biphenyl-4-carboxylic acid",
      "C1=CC=C2C(=C1)C=CC=C2 Naphthalene",
      "CC1=CC=C(C=C1)C Toluene",
    ];
    const exampleText = examples.join("\n");
    setInputText(exampleText);
    setDepictions([]); // Clear previous depictions
    setError(null);
    setShowToolsSection(true); // Show tools when examples are loaded
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Section Card */}
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5, ease: [0.25, 0.46, 0.45, 0.94] }}
        className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-xl dark:shadow-2xl border border-gray-200 dark:border-gray-700"
      >
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Batch 2D Depiction
        </h2>
        <p className="text-sm text-gray-600 dark:text-gray-400 mb-4">
          Generate 2D molecular depictions with RDKit or CDK. CDK offers enhanced features like abbreviations and advanced styling.
        </p>

        <form onSubmit={handleSubmit} className="space-y-4">
          {/* Input Text Area */}
          <div>
            <div className="flex justify-between items-center mb-1">
              <label
                htmlFor="smiles-input-batch"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300"
              >
                Enter SMILES (one per line, optional title after space/tab)
              </label>
              {/* Load Examples Button */}
              <button
                type="button"
                onClick={generateExamples}
                className="text-sm font-medium text-blue-600 hover:text-blue-700 dark:text-blue-400 dark:hover:text-blue-300 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 rounded"
              >
                Load Examples
              </button>
            </div>
            {/* Textarea Styling */}
            <textarea
              id="smiles-input-batch"
              value={inputText}
              onChange={(e) => handleInputChange(e.target.value)}
              className="w-full h-40 font-mono text-xs sm:text-sm px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 dark:focus:border-blue-500 dark:focus:ring-blue-500 resize-y"
              placeholder="CN1C=NC2=C1C(=O)N(C(=O)N2C)C Caffeine&#10;CCO Ethanol&#10;# Lines starting with # are ignored"
              required
              aria-required="true"
            />
            <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
              Format: SMILES [space/tab] Title (optional). Max 50 molecules processed.
            </p>
          </div>

          {/* Tools Section - Animated visibility */}
          <AnimatePresence>
            {showToolsSection && (
              <motion.div
                key="tools-section"
                initial={{ opacity: 0, height: 0, scale: 0.98 }}
                animate={{ 
                  opacity: 1, 
                  height: "auto", 
                  scale: 1,
                  transition: {
                    duration: 0.4,
                    ease: [0.25, 0.46, 0.45, 0.94],
                    opacity: { duration: 0.3 },
                    height: { duration: 0.4 },
                    scale: { duration: 0.3, delay: 0.1 }
                  }
                }}
                exit={{ 
                  opacity: 0, 
                  height: 0,
                  scale: 0.98,
                  transition: {
                    duration: 0.3,
                    ease: "easeIn"
                  }
                }}
                className="space-y-4 overflow-hidden border-t border-gray-200 dark:border-gray-700 pt-4"
              >
                <h3 className="text-md font-medium text-gray-700 dark:text-gray-300">
                  Depiction Options
                </h3>
                {/* Tools Grid */}
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4 px-0.5">
                  {/* Toolkit Select */}
                  <div>
                    <label
                      htmlFor="toolkit-select-batch"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Toolkit
                    </label>
                    <div className="flex items-center">
                      <select
                        id="toolkit-select-batch"
                        value={toolkit}
                        onChange={(e) => handleToolkitChange(e.target.value)}
                        className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                      >
                        <option value="rdkit">RDKit</option>
                        <option value="cdk">CDK (Enhanced)</option>
                      </select>
                      <button
                        type="button"
                        onClick={handleSwitchToolkit}
                        className="ml-2 p-2 bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-600 dark:text-gray-300 rounded-md border border-gray-300 dark:border-gray-600 shadow-sm"
                        title="Switch toolkit"
                      >
                        <HiOutlineSwitchHorizontal className="h-5 w-5" />
                      </button>
                    </div>
                  </div>

                  {/* Width Input */}
                  <div>
                    <label
                      htmlFor="width-input"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Width (px)
                    </label>
                    <div className="flex items-center space-x-2">
                      <motion.button
                        type="button"
                        onClick={() => setWidth(Math.max(50, width - 10))}
                        whileHover={{ scale: 1.05 }}
                        whileTap={{ scale: 0.95 }}
                        className="px-3 py-2 bg-gradient-to-r from-gray-200 to-gray-300 dark:from-gray-600 dark:to-gray-700 hover:from-gray-300 hover:to-gray-400 dark:hover:from-gray-500 dark:hover:to-gray-600 text-gray-700 dark:text-gray-200 font-bold rounded-md border border-gray-400 dark:border-gray-500 shadow-md transition-all duration-200"
                        aria-label="Decrease width"
                      >
                        −
                      </motion.button>
                      <input
                        id="width-input"
                        type="number"
                        value={width}
                        onChange={(e) =>
                          setWidth(Math.max(50, Math.min(2000, Number(e.target.value))))
                        }
                        min="50"
                        max="2000"
                        step="10"
                        className="flex-1 px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm text-center [appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none"
                      />
                      <motion.button
                        type="button"
                        onClick={() => setWidth(Math.min(2000, width + 10))}
                        whileHover={{ scale: 1.05 }}
                        whileTap={{ scale: 0.95 }}
                        className="px-3 py-2 bg-gradient-to-r from-gray-200 to-gray-300 dark:from-gray-600 dark:to-gray-700 hover:from-gray-300 hover:to-gray-400 dark:hover:from-gray-500 dark:hover:to-gray-600 text-gray-700 dark:text-gray-200 font-bold rounded-md border border-gray-400 dark:border-gray-500 shadow-md transition-all duration-200"
                        aria-label="Increase width"
                      >
                        +
                      </motion.button>
                    </div>
                  </div>

                  {/* Height Input */}
                  <div>
                    <label
                      htmlFor="height-input"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Height (px)
                    </label>
                    <div className="flex items-center space-x-2">
                      <motion.button
                        type="button"
                        onClick={() => setHeight(Math.max(50, height - 10))}
                        whileHover={{ scale: 1.05 }}
                        whileTap={{ scale: 0.95 }}
                        className="px-3 py-2 bg-gradient-to-r from-gray-200 to-gray-300 dark:from-gray-600 dark:to-gray-700 hover:from-gray-300 hover:to-gray-400 dark:hover:from-gray-500 dark:hover:to-gray-600 text-gray-700 dark:text-gray-200 font-bold rounded-md border border-gray-400 dark:border-gray-500 shadow-md transition-all duration-200"
                        aria-label="Decrease height"
                      >
                        −
                      </motion.button>
                      <input
                        id="height-input"
                        type="number"
                        value={height}
                        onChange={(e) =>
                          setHeight(Math.max(50, Math.min(2000, Number(e.target.value))))
                        }
                        min="50"
                        max="2000"
                        step="10"
                        className="flex-1 px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm text-center [appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none"
                      />
                      <motion.button
                        type="button"
                        onClick={() => setHeight(Math.min(2000, height + 10))}
                        whileHover={{ scale: 1.05 }}
                        whileTap={{ scale: 0.95 }}
                        className="px-3 py-2 bg-gradient-to-r from-gray-200 to-gray-300 dark:from-gray-600 dark:to-gray-700 hover:from-gray-300 hover:to-gray-400 dark:hover:from-gray-500 dark:hover:to-gray-600 text-gray-700 dark:text-gray-200 font-bold rounded-md border border-gray-400 dark:border-gray-500 shadow-md transition-all duration-200"
                        aria-label="Increase height"
                      >
                        +
                      </motion.button>
                    </div>
                  </div>

                  {/* Hydrogen Display (CDK only) */}
                  {toolkit === "cdk" && (
                    <div>
                      <label
                        htmlFor="hydrogen-display-select"
                        className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                      >
                        Hydrogen Display
                      </label>
                      <select
                        id="hydrogen-display-select"
                        value={hydrogenDisplay}
                        onChange={(e) => setHydrogenDisplay(e.target.value)}
                        className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                      >
                        <option value="Smart">Smart</option>
                        <option value="Provided">Provided</option>
                        <option value="Minimal">Minimal</option>
                        <option value="Explicit">Explicit</option>
                        <option value="Stereo">Stereo</option>
                      </select>
                    </div>
                  )}
                </div>

                {/* CDK Enhanced Options - Animated */}
                <AnimatePresence>
                  {toolkit === "cdk" && (
                    <motion.div
                      initial={{ opacity: 0, height: 0, y: -10 }}
                      animate={{ 
                        opacity: 1, 
                        height: "auto", 
                        y: 0,
                        transition: {
                          duration: 0.5,
                          ease: [0.25, 0.46, 0.45, 0.94],
                          opacity: { duration: 0.4, delay: 0.1 },
                          height: { duration: 0.5 },
                          y: { duration: 0.4, delay: 0.1 }
                        }
                      }}
                      exit={{ 
                        opacity: 0, 
                        height: 0, 
                        y: -10,
                        transition: {
                          duration: 0.3,
                          ease: "easeIn"
                        }
                      }}
                      className="space-y-4 overflow-hidden"
                    >
                      <div className="border-t border-blue-200 dark:border-blue-800 pt-4 mt-2">
                        <h4 className="text-sm font-medium text-blue-700 dark:text-blue-400 mb-3 flex items-center">
                          <span className="mr-2">✨</span>
                          Enhanced CDK Features
                        </h4>
                        
                        <div className="grid grid-cols-1 md:grid-cols-3 gap-4 px-0.5">
                          {/* Style Preset */}
                          <div>
                            <label
                              htmlFor="style-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Color Scheme
                            </label>
                            <select
                              id="style-select"
                              value={style}
                              onChange={(e) => setStyle(e.target.value)}
                              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                            >
                              <option value="cow">Color on White</option>
                              <option value="cob">Color on Black</option>
                              <option value="cot">Color on Transparent</option>
                              <option value="bow">Black on White</option>
                              <option value="bot">Black on Transparent</option>
                              <option value="wob">White on Black</option>
                              <option value="nob">Neon on Black</option>
                            </select>
                          </div>

                          {/* Abbreviations */}
                          <div>
                            <label
                              htmlFor="abbreviate-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Abbreviations
                            </label>
                            <select
                              id="abbreviate-select"
                              value={abbreviate}
                              onChange={(e) => setAbbreviate(e.target.value)}
                              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                            >
                              <option value="off">None</option>
                              <option value="groups">Functional Groups (Ph, Me, Et)</option>
                              <option value="reagents">Reagents (THF, DMF)</option>
                              <option value="on">Both Groups & Reagents</option>
                            </select>
                          </div>

                          {/* Annotations */}
                          <div>
                            <label
                              htmlFor="annotate-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Annotations
                            </label>
                            <select
                              id="annotate-select"
                              value={annotate}
                              onChange={(e) => setAnnotate(e.target.value)}
                              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                            >
                              <option value="none">None</option>
                              <option value="number">Atom Numbers</option>
                              <option value="cip">CIP Labels (R/S, E/Z)</option>
                              <option value="mapidx">Atom Mapping</option>
                              <option value="colmap">Color-Coded Mapping</option>
                              <option value="bondnumber">Bond Numbers</option>
                              <option value="atomvalue">Atom Values</option>
                            </select>
                          </div>

                          {/* Dative Bonds */}
                          <div>
                            <label
                              htmlFor="dative-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Dative Bonds
                            </label>
                            <select
                              id="dative-select"
                              value={dative}
                              onChange={(e) => setDative(e.target.value)}
                              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                            >
                              <option value="never">Never</option>
                              <option value="metals">Metal Complexes</option>
                              <option value="always">Always</option>
                            </select>
                          </div>

                          {/* Multicenter Bonds */}
                          <div>
                            <label
                              htmlFor="multicenter-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Multicenter Bonds
                            </label>
                            <select
                              id="multicenter-select"
                              value={multicenter}
                              onChange={(e) => setMulticenter(e.target.value)}
                              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                            >
                              <option value="provided">As Provided</option>
                              <option value="dative">Dative Arrows</option>
                              <option value="dashed">Dashed Lines</option>
                              <option value="dashed_neutral">Dashed (Neutral)</option>
                              <option value="hidden">Hidden</option>
                              <option value="hidden_neutral">Hidden (Neutral)</option>
                            </select>
                          </div>

                          {/* Zoom */}
                          <div>
                            <label
                              htmlFor="zoom-input"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Zoom: {zoom.toFixed(1)}x
                            </label>
                            <input
                              id="zoom-input"
                              type="range"
                              value={zoom}
                              onChange={(e) => setZoom(Number(e.target.value))}
                              min="0.5"
                              max="3.0"
                              step="0.1"
                              className="w-full h-2 bg-gray-200 dark:bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500 dark:accent-blue-400"
                            />
                          </div>

                          {/* Ratio */}
                          <div>
                            <label
                              htmlFor="ratio-input"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Bond Thickness: {ratio.toFixed(1)}x
                            </label>
                            <input
                              id="ratio-input"
                              type="range"
                              value={ratio}
                              onChange={(e) => setRatio(Number(e.target.value))}
                              min="0.5"
                              max="2.0"
                              step="0.1"
                              className="w-full h-2 bg-gray-200 dark:bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500 dark:accent-blue-400"
                            />
                          </div>
                        </div>

                        {/* CDK Advanced Toggles */}
                        <div className="flex flex-wrap gap-6 pt-4 border-t border-gray-200 dark:border-gray-700 mt-4">
                          <ToggleSwitch
                            id="flip"
                            checked={flip}
                            onChange={setFlip}
                            label="Flip structure horizontally"
                          />
                          
                          <ToggleSwitch
                            id="showtitle"
                            checked={showtitle}
                            onChange={setShowtitle}
                            label="Show titles in depiction"
                          />
                          
                          <ToggleSwitch
                            id="perceive-radicals"
                            checked={perceiveRadicals}
                            onChange={setPerceiveRadicals}
                            label="Perceive radicals"
                          />
                        </div>
                      </div>
                    </motion.div>
                  )}
                </AnimatePresence>

                {/* Highlight Input */}
                <div>
                  <label
                    htmlFor="highlight-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                  >
                    Highlight Substructure (SMARTS - optional)
                  </label>
                  <input
                    id="highlight-input"
                    type="text"
                    value={highlight}
                    onChange={(e) => setHighlight(e.target.value)}
                    placeholder="e.g., c1ccccc1"
                    className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                  />
                </div>

                {/* Checkboxes */}
                <div className="flex flex-wrap gap-6 pt-3">
                  <ToggleSwitch
                    id="cip"
                    checked={showCIP}
                    onChange={setShowCIP}
                    label="Show CIP (R/S, E/Z)"
                  />
                  <ToggleSwitch
                    id="unicolor"
                    checked={useUnicolor}
                    onChange={setUseUnicolor}
                    label="Black & white"
                  />
                  <ToggleSwitch
                    id="atomNumbers"
                    checked={showAtomNumbers}
                    onChange={setShowAtomNumbers}
                    label="Show atom numbers"
                  />
                  {/* Aromatic Donuts Toggle (CDK only) */}
                  {toolkit === "cdk" && (
                    <motion.div
                      initial={{ opacity: 0, scale: 0.9 }}
                      animate={{ opacity: 1, scale: 1 }}
                      exit={{ opacity: 0, scale: 0.9 }}
                      transition={{ duration: 0.2 }}
                    >
                      <ToggleSwitch
                        id="donuts"
                        checked={donuts}
                        onChange={setDonuts}
                        label="Aromatic rings (circles)"
                      />
                    </motion.div>
                  )}
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Action Buttons Row */}
          <div className="flex flex-wrap items-center gap-3 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Generate Button */}
            <motion.button
              type="submit"
              disabled={!inputText.trim() || loading}
              whileHover={!inputText.trim() || loading ? {} : { scale: 1.02, y: -1 }}
              whileTap={!inputText.trim() || loading ? {} : { scale: 0.98 }}
              transition={{ type: "spring", stiffness: 400, damping: 17 }}
              className={`px-5 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !inputText.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed" // Disabled
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-lg hover:shadow-xl" // Enabled
              }`}
            >
              <HiOutlinePhotograph className="mr-2 h-5 w-5" />
              {loading ? "Generating..." : "Generate Depictions"}
            </motion.button>

            {/* Regenerate Button - only shown when depictions exist and not loading */}
            {depictions.length > 0 && !loading && (
              <motion.button
                onClick={() => regenerateDepictions()}
                initial={{ opacity: 0, scale: 0.9 }}
                animate={{ opacity: 1, scale: 1 }}
                exit={{ opacity: 0, scale: 0.9 }}
                whileHover={{ scale: 1.02, y: -1 }}
                whileTap={{ scale: 0.98 }}
                transition={{ type: "spring", stiffness: 400, damping: 17 }}
                className="px-4 py-2 rounded-lg text-gray-700 dark:text-gray-300 bg-gray-100 hover:bg-gray-200 dark:bg-gray-700 dark:hover:bg-gray-600 border border-gray-300 dark:border-gray-600 font-medium flex items-center justify-center transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-gray-500 shadow-md hover:shadow-lg"
                title="Apply current options to all depictions"
              >
                <HiOutlinePhotograph className="mr-2 h-4 w-4" />
                Update All
              </motion.button>
            )}

            {/* Download and Copy buttons - only shown when depictions exist */}
            {depictions.length > 0 && !loading && (
              <>
                {/* Download Format Selector and Button */}
                <div className="flex items-center shadow-md rounded-lg overflow-hidden px-0.5">
                  {/* Select Styling */}
                  <select
                    value={downloadFormat}
                    onChange={(e) => setDownloadFormat(e.target.value)}
                    className="h-full px-3 py-2 bg-white dark:bg-gray-700 border-r border-gray-300 dark:border-gray-600 text-gray-900 dark:text-white text-sm focus:ring-green-500 focus:border-green-500 transition-all duration-200"
                    aria-label="Download format"
                  >
                    <option value="svg">SVG</option>
                    <option value="png">PNG</option>
                  </select>
                  {/* Download Button with Animation */}
                  <motion.button
                    type="button"
                    onClick={downloadAllDepictions}
                    whileHover={{ scale: 1.02 }}
                    whileTap={{ scale: 0.98 }}
                    transition={{ type: "spring", stiffness: 400, damping: 17 }}
                    className="px-4 py-2 bg-green-600 hover:bg-green-700 dark:bg-green-500 dark:hover:bg-green-600 text-white font-medium flex items-center text-sm transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-gray-800 focus:ring-green-500"
                    title="Download all depictions as a ZIP file"
                  >
                    <HiOutlineDownload className="mr-1.5 h-4 w-4" />
                    Download All (.zip)
                  </motion.button>
                </div>

                {/* Copy All SMILES Button */}
                <motion.button
                  type="button"
                  onClick={copyAllSmiles}
                  whileHover={{ scale: 1.02 }}
                  whileTap={{ scale: 0.98 }}
                  transition={{ type: "spring", stiffness: 400, damping: 17 }}
                  className={`px-4 py-2 font-medium rounded-lg flex items-center text-sm transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 shadow-md hover:shadow-lg ${
                    copiedSmiles
                      ? "bg-green-100 dark:bg-green-700 text-green-700 dark:text-green-200"
                      : "bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600"
                  }`}
                  title="Copy all input SMILES to clipboard"
                >
                  <motion.div
                    animate={{ rotate: copiedSmiles ? [0, -10, 10, -10, 0] : 0 }}
                    transition={{ duration: 0.5 }}
                  >
                    {copiedSmiles ? (
                      <HiOutlineCheck className="mr-1.5 h-5 w-5" />
                    ) : (
                      <HiOutlineClipboard className="mr-1.5 h-5 w-5" />
                    )}
                  </motion.div>
                  {copiedSmiles ? "Copied!" : "Copy All SMILES"}
                </motion.button>
              </>
            )}
          </div>
        </form>
      </motion.div>

      {/* Loading Screen */}
      {loading && <LoadingScreen text="Generating depictions..." />}

      {/* Error Display */}
      {error && !loading && (
        // Error message styling
        <div
          className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow"
          role="alert"
        >
          <HiOutlineExclamationCircle
            className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
            aria-hidden="true"
          />
          <span>{error}</span>
        </div>
      )}

      {/* Results Grid */}
      {depictions.length > 0 && !loading && (
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, ease: [0.25, 0.46, 0.45, 0.94] }}
          className="space-y-4"
        >
          {/* Results Header */}
          <motion.h3
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.4, delay: 0.1 }}
            className="text-lg font-semibold text-gray-800 dark:text-gray-200 flex items-center"
          >
            <span className="mr-2">🧪</span>
            Generated Depictions ({depictions.length})
          </motion.h3>

          {/* Grid Container with Animation */}
          <motion.div
            className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4" // Responsive grid
            variants={resultsContainerVariants}
            initial="hidden"
            animate="visible"
          >
            {depictions.map((depiction) => (
              // Individual Card with Animation
              <motion.div
                key={depiction.id}
                className="bg-white dark:bg-gray-800 rounded-xl overflow-hidden shadow-lg dark:shadow-2xl border border-gray-200 dark:border-gray-700 flex flex-col hover:shadow-2xl dark:hover:shadow-3xl transition-shadow duration-300"
                variants={depictionCardVariant}
                layout
                whileHover={{ y: -4, transition: { duration: 0.2 } }}
              >
                {/* Image container - Kept white background for consistency */}
                <div className="p-2 flex-grow">
                  {/* Ensure depiction image background is consistent */}
                  <div className="bg-white rounded-md p-2 flex items-center justify-center h-48">
                    <img
                      src={depiction.imageUrl}
                      alt={depiction.title}
                      className="max-w-full max-h-full object-contain"
                      loading="lazy" // Lazy load images
                      onError={(e) => {
                        console.error(
                          `Error loading image for ${depiction.title}: ${depiction.imageUrl}`
                        );
                        e.target.onerror = null;
                        // Use neutral fallback SVG
                        e.target.src =
                          "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48dGV4dCB4PSI1MCIgeT0iNTAiIGRvbWluYW50LWJhc2VsaW5lPSJtaWRkbGUiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTAiIGZpbGw9IiM4ODg4ODg4Ij5FcnJvciBsb2FkaW5nPC90ZXh0Pjwvc3ZnPg==";
                        e.target.classList.add("p-4"); // Add padding to text
                      }}
                    />
                  </div>
                </div>

                {/* Info Section */}
                <div className="px-3 py-2 border-t border-gray-200 dark:border-gray-700">
                  <h4
                    className="font-medium text-gray-800 dark:text-gray-100 text-sm truncate"
                    title={depiction.title}
                  >
                    {depiction.title}
                  </h4>
                  <p
                    className="text-gray-500 dark:text-gray-400 text-xs truncate font-mono"
                    title={depiction.smiles}
                  >
                    {depiction.smiles}
                  </p>
                </div>

                {/* Rotation Slider */}
                <div className="px-3 py-2 bg-gray-50 dark:bg-gray-900 border-t border-gray-200 dark:border-gray-700">
                  <label
                    htmlFor={`rotate-${depiction.id}`}
                    className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1"
                  >
                    Rotate: {rotations[depiction.id] || 0}°
                  </label>
                  {/* Slider Styling */}
                  <input
                    id={`rotate-${depiction.id}`}
                    type="range"
                    value={rotations[depiction.id] || 0}
                    onChange={(e) => handleRotationChange(depiction.id, e.target.value)}
                    min="0"
                    max="359"
                    step="1"
                    // Enhanced slider styling
                    className="w-full h-2 bg-gray-200 dark:bg-gray-700 rounded-lg appearance-none cursor-pointer range-sm accent-blue-500 dark:accent-blue-400 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-900 focus:ring-blue-500"
                    aria-label={`Rotation for ${depiction.title}`}
                  />
                  {/* Slider Labels */}
                  <div className="flex justify-between text-xs text-gray-400 dark:text-gray-500 mt-1">
                    <span>0°</span>
                    <span>359°</span>
                  </div>
                </div>

                {/* Action Buttons */}
                <div className="px-3 py-2 bg-gradient-to-r from-gray-50 to-gray-100 dark:from-gray-700 dark:to-gray-800 border-t border-gray-200 dark:border-gray-600 flex justify-end space-x-2">
                  {/* Download Button */}
                  <motion.button
                    onClick={() => downloadSingleDepiction(depiction)}
                    whileHover={{ scale: 1.1, y: -2 }}
                    whileTap={{ scale: 0.95 }}
                    transition={{ type: "spring", stiffness: 400, damping: 17 }}
                    className="p-2 rounded-lg text-gray-600 dark:text-gray-400 hover:text-blue-600 dark:hover:text-blue-400 hover:bg-blue-50 dark:hover:bg-blue-900/30 transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-blue-500 shadow-sm hover:shadow-md"
                    title={`Download ${downloadFormat.toUpperCase()}`}
                    aria-label={`Download ${depiction.title} as ${downloadFormat.toUpperCase()}`}
                    disabled={loading}
                  >
                    <HiOutlineDownload className="h-5 w-5" />
                  </motion.button>

                  {/* Copy SMILES Button */}
                  <motion.button
                    onClick={() => copySingleSmiles(depiction.smiles, depiction.id)}
                    whileHover={{ scale: 1.1, y: -2 }}
                    whileTap={{ scale: 0.95 }}
                    transition={{ type: "spring", stiffness: 400, damping: 17 }}
                    className="p-2 rounded-lg text-gray-600 dark:text-gray-400 hover:text-green-600 dark:hover:text-green-400 hover:bg-green-50 dark:hover:bg-green-900/30 relative transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-green-500 shadow-sm hover:shadow-md"
                    title="Copy SMILES"
                    aria-label={`Copy SMILES for ${depiction.title}`}
                  >
                    {/* Animated Check/Clipboard Icon */}
                    <AnimatePresence mode="wait">
                      {copiedSingle === depiction.id ? (
                        <motion.div
                          key="check"
                          initial={{ scale: 0, rotate: -180 }}
                          animate={{ scale: 1, rotate: 0 }}
                          exit={{ scale: 0, rotate: 180 }}
                          transition={{ 
                            type: "spring", 
                            stiffness: 500, 
                            damping: 25,
                            duration: 0.3
                          }}
                          className="flex items-center justify-center"
                        >
                          <HiOutlineCheck className="h-5 w-5 text-green-500 dark:text-green-400" />
                        </motion.div>
                      ) : (
                        <motion.div
                          key="clipboard"
                          initial={{ scale: 1, rotate: 0 }}
                          exit={{ scale: 0, rotate: -180 }}
                          transition={{ duration: 0.2 }}
                          className="flex items-center justify-center"
                        >
                          <HiOutlineClipboard className="h-5 w-5" />
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </motion.button>
                </div>
              </motion.div>
            ))}
          </motion.div>
        </motion.div>
      )}

      {/* Initial Placeholder */}
      {!depictions.length && !loading && !error && (
        // Initial state card styling
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg flex flex-col items-center justify-center text-center min-h-[300px] border border-gray-200 dark:border-gray-700">
          <HiOutlinePhotograph className="h-16 w-16 text-gray-400 dark:text-gray-600 mb-4" />
          <p className="text-gray-500 dark:text-gray-400">
            Enter SMILES strings (one per line, optionally with titles) and click "Generate
            Depictions".
          </p>
          <button
            type="button"
            onClick={generateExamples}
            className="mt-4 text-sm font-medium text-blue-600 hover:text-blue-700 dark:text-blue-400 dark:hover:text-blue-300 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 rounded"
          >
            Or load example molecules
          </button>
        </div>
      )}
    </div>
  );
};

export default BatchDepictionView;
