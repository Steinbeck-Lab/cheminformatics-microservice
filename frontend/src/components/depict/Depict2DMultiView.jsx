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
import { motion, AnimatePresence } from "framer-motion";
// Assuming these components are correctly implemented and styled for dark/light mode
import LoadingScreen from "../common/LoadingScreen";
// Assuming this service is configured correctly
import depictService from "../../services/depictService"; // Assuming this service exists

// Animation variants
const resultsContainerVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { staggerChildren: 0.05 } }, // Stagger children animation
};

const depictionCardVariant = {
  hidden: { opacity: 0, y: 20, scale: 0.98 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { duration: 0.4, ease: "easeOut" },
  },
};

const BatchDepictionView = () => {
  // Input state
  const [inputText, setInputText] = useState("");

  // Processing state
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [depictions, setDepictions] = useState([]); // Array of { smiles, title, imageUrl, id }

  // Depiction options
  const [toolkit, setToolkit] = useState("rdkit");
  const [width, setWidth] = useState(300);
  const [height, setHeight] = useState(200);
  const [showCIP, setShowCIP] = useState(false);
  const [useUnicolor, setUseUnicolor] = useState(false);
  const [highlight, setHighlight] = useState("");
  const [showAtomNumbers, setShowAtomNumbers] = useState(false);

  // UI state
  const [copiedSmiles, setCopiedSmiles] = useState(false); // For "Copy All SMILES" button
  const [copiedSingle, setCopiedSingle] = useState(null); // Tracks ID of molecule whose SMILES was copied
  const [downloadFormat, setDownloadFormat] = useState("svg"); // svg or png
  const [showToolsSection, setShowToolsSection] = useState(false); // Toggle visibility of options

  // Per-molecule rotation state { [id]: rotationValue }
  const [rotations, setRotations] = useState({});

  // Parse SMILES (and optional titles) from input text
  const parseSmiles = (text) => {
    return text
      .split(/[\n\r]+/) // Split by new lines
      .map((line) => line.trim()) // Trim whitespace
      .filter((line) => line.length > 0 && !line.startsWith("#")); // Remove empty lines and comments
  };

  // Regenerate depiction URLs when toolkit changes (or called manually)
  const regenerateDepictions = (currentToolkit = toolkit) => {
    // Check if depictService and the method exist before proceeding
    if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
      console.error("depictService.get2DDepictionUrl is not available.");
      setError("Depiction service is not configured correctly.");
      return;
    }
    if (depictions.length === 0) return;

    console.log(`Regenerating depictions with toolkit: ${currentToolkit}`); // Debug log

    setDepictions((currentDepictions) =>
      currentDepictions.map((dep) => {
        const rotation = rotations[dep.id] || 0;
        const options = {
          toolkit: currentToolkit, // Use potentially updated toolkit
          width,
          height,
          rotate: rotation,
          CIP: currentToolkit === "cdk" ? showCIP : undefined, // Check updated toolkit
          unicolor: useUnicolor,
          highlight: highlight || undefined,
          showAtomNumbers,
        };
        // DEBUG: Log options and URL generation
        // console.log(`Regen options for ${dep.smiles}:`, options);
        const updatedImageUrl = depictService.get2DDepictionUrl(dep.smiles, options);
        // console.log(`Regen URL for ${dep.smiles}:`, updatedImageUrl);
        return { ...dep, imageUrl: updatedImageUrl };
      })
    );
  };

  // Handle toolkit change, reset incompatible options, and regenerate
  const handleToolkitChange = (newToolkit) => {
    setToolkit(newToolkit);
    // RDKit doesn't support CIP labels in this context
    if (newToolkit === "rdkit") {
      setShowCIP(false);
    }
    // Regenerate depictions if results already exist
    if (depictions.length > 0) {
      regenerateDepictions(newToolkit);
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
    // Check if depictService and the method exist before proceeding
    if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
      console.error("depictService.get2DDepictionUrl is not available for rotation.");
      setError("Depiction service is not configured correctly.");
      return;
    }
    setDepictions((currentDepictions) =>
      currentDepictions.map((dep) => {
        if (dep.id === id) {
          const options = {
            toolkit,
            width,
            height,
            rotate: rotation,
            CIP: toolkit === "cdk" ? showCIP : undefined,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
          };
          const updatedImageUrl = depictService.get2DDepictionUrl(dep.smiles, options);
          return { ...dep, imageUrl: updatedImageUrl };
        }
        return dep;
      })
    );
  };

  // Handle form submission to generate initial depictions
  const handleSubmit = async (e) => {
    e.preventDefault();

    // Check if depictService and the method exist before proceeding
    if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
      console.error("depictService.get2DDepictionUrl is not available.");
      setError("Depiction service is not configured correctly.");
      return;
    }

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

        // Prepare options for initial depiction using CURRENT state values
        const options = {
          toolkit,
          width,
          height,
          rotate: 0, // Initial rotation
          CIP: toolkit === "cdk" ? showCIP : undefined,
          unicolor: useUnicolor,
          highlight: highlight || undefined,
          showAtomNumbers,
        };

        // Get the URL (assuming service returns URL directly)
        const imageUrl = depictService.get2DDepictionUrl(smiles, options);

        results.push({ smiles, title, imageUrl, id });
      }

      setRotations(initialRotations); // Set all initial rotations at once
      setDepictions(results); // Set all results at once

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

      // Check if depictService and the method exist before proceeding
      if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
        throw new Error("Depiction service is not configured correctly for download.");
      }

      // Fetch all depiction blobs concurrently
      const fetchPromises = depictions.map(async (depiction) => {
        const safeTitle =
          depiction.title.replace(/[^a-z0-9]/gi, "_").toLowerCase() || `molecule_${depiction.id}`;
        const filename = `${safeTitle}.${downloadFormat}`; // Use selected format
        const rotation = rotations[depiction.id] || 0; // Get current rotation

        // Generate URL with current options and selected download format
        const options = {
          toolkit,
          width,
          height,
          rotate: rotation,
          CIP: toolkit === "cdk" ? showCIP : undefined,
          unicolor: useUnicolor,
          highlight: highlight || undefined,
          showAtomNumbers,
          format: downloadFormat, // Pass format to service if needed
        };

        // Use a potentially different URL for download if format differs or specific endpoint exists
        const downloadUrl = depictService.get2DDepictionUrl(depiction.smiles, options);

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
      // Check if depictService and the method exist before proceeding
      if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
        console.error("depictService.get2DDepictionUrl is not available for single download.");
        setError("Depiction service is not configured correctly.");
        return;
      }

      setLoading(true); // Show loading indicator while downloading
      const rotation = rotations[depiction.id] || 0;

      // Generate URL with current options and selected download format
      const options = {
        toolkit,
        width,
        height,
        rotate: rotation,
        CIP: toolkit === "cdk" ? showCIP : undefined,
        unicolor: useUnicolor,
        highlight: highlight || undefined,
        showAtomNumbers,
        format: downloadFormat, // Use selected format for single download too
      };

      // Get the URL for fetching the image
      const url = depictService.get2DDepictionUrl(depiction.smiles, options);

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
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Batch 2D Depiction
        </h2>

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
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: "auto" }}
                exit={{ opacity: 0, height: 0 }}
                transition={{ duration: 0.3, ease: "easeInOut" }}
                className="space-y-4 overflow-hidden border-t border-gray-200 dark:border-gray-700 pt-4"
              >
                <h3 className="text-md font-medium text-gray-700 dark:text-gray-300">
                  Depiction Options
                </h3>
                {/* Tools Grid */}
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  {/* Toolkit Select */}
                  <div>
                    <label
                      htmlFor="toolkit-select-batch"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Toolkit
                    </label>
                    <div className="flex items-center">
                      {/* Select Styling */}
                      <select
                        id="toolkit-select-batch"
                        value={toolkit}
                        onChange={(e) => handleToolkitChange(e.target.value)}
                        className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                      >
                        <option value="rdkit">RDKit</option>
                        <option value="cdk">CDK</option>
                      </select>
                      {/* Switch Button Styling */}
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
                    {/* Input Styling */}
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
                      className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                    />
                  </div>

                  {/* Height Input */}
                  <div>
                    <label
                      htmlFor="height-input"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Height (px)
                    </label>
                    {/* Input Styling */}
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
                      className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                    />
                  </div>
                </div>

                {/* Highlight Input */}
                <div>
                  <label
                    htmlFor="highlight-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                  >
                    Highlight Substructure (SMARTS - optional)
                  </label>
                  {/* Input Styling */}
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
                <div className="flex flex-wrap gap-x-6 gap-y-2 pt-2">
                  {/* CIP Checkbox */}
                  <div className="flex items-center">
                    <input
                      id="cip"
                      type="checkbox"
                      checked={showCIP}
                      onChange={(e) => setShowCIP(e.target.checked)}
                      disabled={toolkit === "rdkit"}
                      // Checkbox Styling (including disabled state)
                      className={`h-4 w-4 rounded border-gray-300 dark:border-gray-600 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 ${
                        toolkit === "rdkit"
                          ? "bg-gray-100 dark:bg-gray-700 opacity-50 cursor-not-allowed"
                          : "bg-white dark:bg-gray-700 text-blue-600 dark:text-blue-500"
                      }`}
                    />
                    <label
                      htmlFor="cip"
                      // Label Styling (including disabled state)
                      className={`ml-2 text-sm ${toolkit === "rdkit" ? "text-gray-400 dark:text-gray-500" : "text-gray-700 dark:text-gray-300"}`}
                    >
                      Show CIP (R/S, E/Z)
                      {toolkit === "rdkit" && <span className="ml-2 text-xs">(CDK only)</span>}
                    </label>
                  </div>
                  {/* Unicolor Checkbox */}
                  <div className="flex items-center">
                    <input
                      id="unicolor"
                      type="checkbox"
                      checked={useUnicolor}
                      onChange={(e) => setUseUnicolor(e.target.checked)}
                      // Checkbox Styling
                      className="h-4 w-4 rounded border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
                    />
                    <label
                      htmlFor="unicolor"
                      className="ml-2 text-sm text-gray-700 dark:text-gray-300"
                    >
                      Use black & white color scheme
                    </label>
                  </div>
                  {/* Atom Numbers Checkbox */}
                  <div className="flex items-center">
                    <input
                      id="atomNumbers"
                      type="checkbox"
                      checked={showAtomNumbers}
                      onChange={(e) => setShowAtomNumbers(e.target.checked)}
                      // Checkbox Styling
                      className="h-4 w-4 rounded border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
                    />
                    <label
                      htmlFor="atomNumbers"
                      className="ml-2 text-sm text-gray-700 dark:text-gray-300"
                    >
                      Show atom numbers
                    </label>
                  </div>
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Action Buttons Row */}
          <div className="flex flex-wrap items-center gap-3 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Generate Button */}
            <button
              type="submit"
              disabled={!inputText.trim() || loading}
              // Button Styling
              className={`px-5 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !inputText.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed" // Disabled
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm" // Enabled
              }`}
            >
              <HiOutlinePhotograph className="mr-2 h-5 w-5" />
              {loading ? "Generating..." : "Generate Depictions"}
            </button>

            {/* Regenerate Button - only shown when depictions exist and not loading */}
            {depictions.length > 0 && !loading && (
              <button
                onClick={() => regenerateDepictions()}
                // Button Styling
                className="px-4 py-2 rounded-lg text-gray-700 dark:text-gray-300 bg-gray-100 hover:bg-gray-200 dark:bg-gray-700 dark:hover:bg-gray-600 border border-gray-300 dark:border-gray-600 font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-gray-500 shadow-sm"
                title="Apply current options to all depictions"
              >
                <HiOutlineSwitchHorizontal className="mr-2 h-4 w-4" />
                Update All
              </button>
            )}

            {/* Download and Copy buttons - only shown when depictions exist */}
            {depictions.length > 0 && !loading && (
              <>
                {/* Download Format Selector and Button */}
                <div className="flex items-center shadow-sm rounded-md">
                  {/* Select Styling */}
                  <select
                    value={downloadFormat}
                    onChange={(e) => setDownloadFormat(e.target.value)}
                    className="h-full px-3 py-2 bg-white dark:bg-gray-700 border border-r-0 border-gray-300 dark:border-gray-600 rounded-l-md text-gray-900 dark:text-white text-sm focus:ring-indigo-500 focus:border-indigo-500"
                    aria-label="Download format"
                  >
                    <option value="svg">SVG</option>
                    <option value="png">PNG</option>
                  </select>
                  {/* Download Button Styling */}
                  <button
                    type="button"
                    onClick={downloadAllDepictions}
                    className="px-4 py-2 bg-green-600 hover:bg-green-700 dark:bg-green-500 dark:hover:bg-green-600 text-white font-medium rounded-r-md flex items-center text-sm transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-gray-800 focus:ring-green-500"
                    title="Download all depictions as a ZIP file"
                  >
                    <HiOutlineDownload className="mr-1.5 h-4 w-4" />
                    Download All (.zip)
                  </button>
                </div>

                {/* Copy All SMILES Button */}
                <button
                  type="button"
                  onClick={copyAllSmiles}
                  // Button Styling with copied state
                  className={`px-4 py-2 font-medium rounded-md flex items-center text-sm transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 ${
                    copiedSmiles
                      ? "bg-green-100 dark:bg-green-700 text-green-700 dark:text-green-200"
                      : "bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600 shadow-sm"
                  }`}
                  title="Copy all input SMILES to clipboard"
                >
                  {copiedSmiles ? (
                    <HiOutlineCheck className="mr-1.5 h-5 w-5" />
                  ) : (
                    <HiOutlineClipboard className="mr-1.5 h-5 w-5" />
                  )}
                  {copiedSmiles ? "Copied!" : "Copy All SMILES"}
                </button>
              </>
            )}
          </div>
        </form>
      </div>

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
        <div className="space-y-4">
          {/* Results Header */}
          <h3 className="text-lg font-medium text-gray-800 dark:text-gray-200">
            Generated Depictions ({depictions.length})
          </h3>

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
                // Card Styling
                className="bg-white dark:bg-gray-800 rounded-lg overflow-hidden shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700 flex flex-col"
                variants={depictionCardVariant}
                layout // Animate layout changes smoothly
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
                <div className="px-3 py-1.5 bg-gray-100 dark:bg-gray-700 border-t border-gray-200 dark:border-gray-600 flex justify-end space-x-1">
                  {/* Download Button */}
                  <button
                    onClick={() => downloadSingleDepiction(depiction)}
                    className="p-1.5 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500"
                    title={`Download ${downloadFormat.toUpperCase()}`}
                    aria-label={`Download ${depiction.title} as ${downloadFormat.toUpperCase()}`}
                    disabled={loading}
                  >
                    <HiOutlineDownload className="h-5 w-5" />
                  </button>

                  {/* Copy SMILES Button */}
                  <button
                    onClick={() => copySingleSmiles(depiction.smiles, depiction.id)}
                    className="p-1.5 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-200 dark:hover:bg-gray-600 relative transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500"
                    title="Copy SMILES"
                    aria-label={`Copy SMILES for ${depiction.title}`}
                  >
                    {/* Animated Check/Clipboard Icon */}
                    <AnimatePresence mode="wait">
                      {copiedSingle === depiction.id ? (
                        <motion.div
                          key="check"
                          initial={{ scale: 0.5, opacity: 0 }}
                          animate={{ scale: 1, opacity: 1 }}
                          exit={{ scale: 0.5, opacity: 0 }}
                          transition={{ duration: 0.15 }}
                          className="flex items-center justify-center" // Ensure icon is centered
                        >
                          <HiOutlineCheck className="h-5 w-5 text-green-500 dark:text-green-400" />
                        </motion.div>
                      ) : (
                        <motion.div
                          key="clipboard"
                          initial={{ scale: 1, opacity: 1 }}
                          exit={{ scale: 0.5, opacity: 0 }}
                          transition={{ duration: 0.1 }}
                          className="flex items-center justify-center"
                        >
                          <HiOutlineClipboard className="h-5 w-5" />
                        </motion.div>
                      )}
                    </AnimatePresence>
                  </button>
                </div>
              </motion.div>
            ))}
          </motion.div>
        </div>
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
