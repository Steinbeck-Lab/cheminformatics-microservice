// Description: StructureGenView component for generating chemical structures from molecular formulas.
import React, { useState } from "react";
import {
  HiOutlineBeaker,
  HiOutlineDocumentDuplicate,
  HiOutlineDownload,
  HiOutlineInformationCircle,
  HiOutlineCheck,
  HiOutlineXCircle,
  HiOutlineExclamationCircle, // FIX: Added missing import
} from "react-icons/hi";
import { FaAtom } from "react-icons/fa";
import { motion, AnimatePresence } from "framer-motion";
import LoadingScreen from "../common/LoadingScreen";
import MoleculeCard from "../common/MoleculeCard";
import toolsService from "../../services/toolsService";

// --- Animation Variants ---
const containerVariant = {
  // For overall component fade-in
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.5, ease: "easeOut" } },
};
const contentSectionVariant = {
  // For main sections like form, results
  hidden: { opacity: 0, y: 20 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.6, ease: [0.25, 1, 0.5, 1] },
  },
};
const resultsContainerVariants = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.05, delayChildren: 0.1 },
  },
};
const resultCardVariant = {
  // Renamed from depictionCardVariant
  hidden: { opacity: 0, y: 15, scale: 0.97 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: { duration: 0.4, ease: "easeOut" },
  },
};

const StructureGenView = () => {
  // --- State ---
  const [formula, setFormula] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [structures, setStructures] = useState([]);
  const [selectedStructures, setSelectedStructures] = useState([]);
  const [generationResult, setGenerationResult] = useState(null); // Store full result with metadata

  // --- Functions ---
  const handleGenerateStructures = async (e) => {
    e.preventDefault();
    if (!formula) {
      setError("Please enter a molecular formula");
      return;
    }
    setIsLoading(true);
    setError(null);
    setStructures([]);
    setSelectedStructures([]);
    setGenerationResult(null);

    try {
      const result = await toolsService.generateStructures(formula);
      // Handle different API response formats
      if (typeof result === "string") {
        setError(result);
        setStructures([]);
        setGenerationResult(null);
      } else if (result?.output?.structures && Array.isArray(result.output.structures)) {
        setStructures(result.output.structures);
        setGenerationResult(result.output);
      } else if (result?.output && Array.isArray(result.output)) {
        // Legacy format - just structures array
        setStructures(result.output);
        setGenerationResult(null);
      } else if (Array.isArray(result)) {
        // Direct array format
        setStructures(result);
        setGenerationResult(null);
      } else {
        console.error("Unexpected API response format:", result);
        setError("Unexpected response format from API");
        setStructures([]);
        setGenerationResult(null);
      }
    } catch (err) {
      console.error("Error generating structures:", err);
      setError(`Error generating structures: ${err.message || "Unknown error"}`);
      setStructures([]);
      setGenerationResult(null);
    } finally {
      setIsLoading(false);
    }
  };

  const toggleSelectStructure = (smiles) => {
    setSelectedStructures((prev) =>
      prev.includes(smiles) ? prev.filter((s) => s !== smiles) : [...prev, smiles]
    );
  };

  const downloadSelectedStructures = () => {
    if (selectedStructures.length === 0) return;
    const content = selectedStructures.join("\n");
    const blob = new Blob([content], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${formula || "generated"}-structures.smi`; // Changed extension to smi
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };
  // --- End Functions ---

  return (
    <motion.div
      className="space-y-6 lg:space-y-8"
      variants={containerVariant} // Animate overall container
      initial="hidden"
      animate="visible"
    >
      {/* Form Container - Apply glass, adaptive styles */}
      <motion.div
        className="glass p-5 sm:p-6 rounded-xl shadow-lg border border-[var(--border-secondary)] dark:border-[var(--border-primary)]"
        variants={contentSectionVariant} // Animate form section
      >
        {/* Use adaptive text color */}
        <h2 className="text-xl font-semibold text-[var(--text-accent)] mb-4">
          Structure Generation (from Formula)
        </h2>

        <form onSubmit={handleGenerateStructures} className="space-y-4">
          {/* Input */}
          <div>
            {/* Use adaptive label color */}
            <label
              htmlFor="formula-input"
              className="block text-sm font-medium text-[var(--text-secondary)] mb-1"
            >
              Molecular Formula
            </label>
            {/* Form input with proper dark mode styling */}
            <input
              id="formula-input"
              type="text"
              value={formula}
              onChange={(e) => setFormula(e.target.value)}
              placeholder="e.g., C6H6"
              className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md bg-white dark:bg-gray-700 text-gray-900 dark:text-white placeholder-gray-500 dark:placeholder-gray-400 focus:ring-2 focus:ring-blue-500 focus:border-blue-500 transition-colors duration-200"
              required
            />
            {/* Use adaptive text color */}
            <p className="mt-1 text-xs text-slate-500 dark:text-slate-400">
              Enter a molecular formula (e.g., C6H6, C8H10). Max 10 heavy atoms. Results limited to
              first 1000 structures when total exceeds this limit.
            </p>
          </div>

          {/* Buttons */}
          <div className="flex flex-wrap gap-3 pt-2">
            {/* Use adaptive button classes */}
            <button
              type="submit"
              disabled={!formula || isLoading}
              className="btn btn-primary disabled:opacity-60 disabled:cursor-not-allowed inline-flex items-center"
            >
              <HiOutlineBeaker className="mr-2 h-5 w-5" />
              {isLoading ? "Generating..." : "Generate Structures"}
            </button>

            {/* Show selection buttons only when structures are present */}
            <AnimatePresence>
              {Array.isArray(structures) && structures.length > 0 && (
                <motion.div
                  className="flex flex-wrap gap-3"
                  initial={{ opacity: 0, y: 10 }}
                  animate={{ opacity: 1, y: 0 }}
                  exit={{ opacity: 0 }}
                  transition={{ delay: 0.1 }}
                >
                  <button
                    type="button"
                    onClick={downloadSelectedStructures}
                    disabled={selectedStructures.length === 0}
                    className="btn btn-secondary btn-sm disabled:opacity-60 disabled:cursor-not-allowed inline-flex items-center" // Use btn-secondary
                  >
                    <HiOutlineDownload className="mr-1.5 h-4 w-4" />
                    Download Selected ({selectedStructures.length})
                  </button>
                  <button
                    type="button"
                    onClick={() => setSelectedStructures([...structures])}
                    className="btn btn-secondary btn-sm inline-flex items-center" // Use btn-secondary
                  >
                    <HiOutlineDocumentDuplicate className="mr-1.5 h-4 w-4" />
                    Select All
                  </button>
                  <button
                    type="button"
                    onClick={() => setSelectedStructures([])}
                    disabled={selectedStructures.length === 0}
                    className="btn btn-secondary btn-sm disabled:opacity-60 disabled:cursor-not-allowed inline-flex items-center" // Use btn-secondary
                  >
                    <HiOutlineXCircle className="mr-1.5 h-4 w-4" /> {/* Changed Icon */}
                    Deselect All
                  </button>
                </motion.div>
              )}
            </AnimatePresence>
          </div>
        </form>
      </motion.div>

      {/* Loading Screen */}
      <AnimatePresence>
        {isLoading && <LoadingScreen text="Generating structures..." />}
      </AnimatePresence>

      {/* Error Display - Adaptive */}
      {error && (
        <motion.div
          className="p-4 bg-red-100 dark:bg-red-900/60 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700/60 rounded-lg text-sm flex items-start gap-2"
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
        >
          {/* FIX: Use imported HiOutlineExclamationCircle */}
          <HiOutlineExclamationCircle className="h-5 w-5 flex-shrink-0 mt-0.5" />
          <span>
            <span className="font-medium">Error:</span> {error}
          </span>
        </motion.div>
      )}

      {/* Results Grid Section - Animate entrance */}
      {Array.isArray(structures) && structures.length > 0 && !isLoading && (
        <motion.div
          className="space-y-4"
          variants={contentSectionVariant} // Use section variant
          initial="hidden"
          animate="visible"
        >
          {/* Generation Summary */}
          {generationResult && (
            <motion.div
              className="glass p-4 rounded-lg border border-[var(--border-secondary)] dark:border-[var(--border-primary)]"
              initial={{ opacity: 0, y: 10 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ delay: 0.1 }}
            >
              <h3 className="text-lg font-medium text-[var(--text-primary)] mb-3">
                Generation Summary
              </h3>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4 mb-4">
                <div className="text-center">
                  <div className="text-2xl font-bold text-[var(--text-accent)]">
                    {generationResult.total_count?.toLocaleString()}
                  </div>
                  <div className="text-sm text-[var(--text-secondary)]">Total Possible</div>
                </div>
                <div className="text-center">
                  <div className="text-2xl font-bold text-[var(--text-accent)]">
                    {generationResult.generated_count?.toLocaleString()}
                  </div>
                  <div className="text-sm text-[var(--text-secondary)]">Generated</div>
                </div>
                <div className="text-center">
                  <div className="text-2xl font-bold text-[var(--text-accent)]">
                    {selectedStructures.length}
                  </div>
                  <div className="text-sm text-[var(--text-secondary)]">Selected</div>
                </div>
                <div className="text-center">
                  <div
                    className={`text-sm px-2 py-1 rounded ${generationResult.limit_applied ? "bg-orange-100 dark:bg-orange-900/30 text-orange-700 dark:text-orange-300" : "bg-green-100 dark:bg-green-900/30 text-green-700 dark:text-green-300"}`}
                  >
                    {generationResult.limit_applied ? "Limited" : "Complete"}
                  </div>
                  <div className="text-sm text-[var(--text-secondary)]">Results</div>
                </div>
              </div>

              {/* Generation Settings */}
              {generationResult.settings && (
                <div>
                  <h4 className="text-sm font-medium text-[var(--text-primary)] mb-2">
                    Generation Settings
                  </h4>
                  <div className="space-y-1">
                    {Object.entries(generationResult.settings).map(([flag, description]) => (
                      <div key={flag} className="text-xs flex">
                        <code className="bg-slate-200 dark:bg-slate-700 px-1.5 py-0.5 rounded mr-2 font-mono text-[var(--text-accent)]">
                          {flag}
                        </code>
                        <span className="text-[var(--text-secondary)]">{description}</span>
                      </div>
                    ))}
                  </div>
                </div>
              )}

              {generationResult.limit_applied && (
                <div className="mt-3 p-3 bg-orange-50 dark:bg-orange-900/20 border border-orange-200 dark:border-orange-800/50 rounded text-sm">
                  <HiOutlineInformationCircle className="inline h-4 w-4 mr-1.5 text-orange-600 dark:text-orange-400" />
                  <span className="text-orange-700 dark:text-orange-300">
                    Results limited to first {generationResult.generated_count} structures out of{" "}
                    {generationResult.total_count?.toLocaleString()} total possible isomers.
                  </span>
                </div>
              )}
            </motion.div>
          )}

          <div className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-2">
            {/* Use adaptive text color */}
            <h3 className="text-lg font-medium text-[var(--text-primary)]">
              Generated Structures ({structures.length})
            </h3>
            {/* Use adaptive text color */}
            <div className="text-[var(--text-secondary)] flex items-center text-sm flex-shrink-0">
              <FaAtom className="mr-1.5 h-4 w-4" />
              Formula:{" "}
              <code className="ml-1.5 text-xs bg-slate-200 dark:bg-slate-700 px-1.5 py-0.5 rounded">
                {formula}
              </code>
            </div>
          </div>

          {/* Apply stagger container to grid */}
          <motion.div
            className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-4" // Responsive grid
            variants={resultsContainerVariants}
            initial="hidden"
            animate="visible"
          >
            {structures.map((smiles, index) => (
              // Animate each card
              <motion.div
                key={smiles + index} // Use SMILES in key for better stability if list order changes
                variants={resultCardVariant}
                layout // Animate layout changes if selection modifies appearance drastically
                // Apply adaptive border for selection state
                className={`rounded-lg transition-all duration-200 cursor-pointer relative group ${
                  selectedStructures.includes(smiles)
                    ? "ring-2 ring-offset-2 ring-offset-[var(--bg-primary)] ring-[var(--text-accent)]" // Use ring for selection
                    : "ring-1 ring-transparent hover:ring-slate-300 dark:hover:ring-slate-600" // Subtle hover ring
                }`}
                onClick={() => toggleSelectStructure(smiles)}
              >
                {/* Ensure MoleculeCard is theme-aware and handles clicks appropriately */}
                <MoleculeCard
                  smiles={smiles}
                  title={`Isomer ${index + 1}`}
                  showActions={false} // Disable internal actions if selection is handled here
                  className="pointer-events-none" // Prevent card actions interfering with selection click
                />
                {/* Selection Checkmark Overlay */}
                <AnimatePresence>
                  {selectedStructures.includes(smiles) && (
                    <motion.div
                      initial={{ opacity: 0, scale: 0.5 }}
                      animate={{ opacity: 1, scale: 1 }}
                      exit={{ opacity: 0, scale: 0.5 }}
                      className="absolute top-2 right-2 z-10 p-1 bg-sky-600 dark:bg-blue-600 rounded-full shadow-md"
                    >
                      <HiOutlineCheck className="h-4 w-4 text-white" />
                    </motion.div>
                  )}
                </AnimatePresence>
              </motion.div>
            ))}
          </motion.div>
        </motion.div>
      )}

      {/* Initial Placeholder / Info - Adaptive */}
      {/* This section is preserved */}
      {Array.isArray(structures) && structures.length === 0 && !isLoading && !error && (
        // Use adaptive styles
        <motion.div
          className="bg-sky-50 dark:bg-blue-900/20 border border-sky-200 dark:border-blue-800/50 rounded-lg p-6 text-slate-600 dark:text-gray-300 text-sm flex items-start space-x-4"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ delay: 0.3 }} // Simple fade-in
        >
          <HiOutlineInformationCircle className="h-6 w-6 text-sky-700 dark:text-blue-400 flex-shrink-0 mt-0.5" />
          <div>
            <h3 className="text-base font-medium text-sky-800 dark:text-blue-300 mb-1.5">
              About Structure Generation
            </h3>
            <p className="mb-2">
              This tool uses the SURGE algorithm to generate possible chemical structures matching
              the given molecular formula.
            </p>
            <p className="mb-2 text-xs">
              Due to computational complexity, generation is limited (e.g., max 10 heavy atoms).
            </p>
            <p className="text-xs text-slate-500 dark:text-gray-400">
              Ref: McKay, B.D., Yirik, M.A. &amp; Steinbeck, C. Surge: a fast open-source chemical
              graph generator. J Cheminform 14, 24 (2022).
            </p>
          </div>
        </motion.div>
      )}
    </motion.div> // End main container div
  );
};

export default StructureGenView;
