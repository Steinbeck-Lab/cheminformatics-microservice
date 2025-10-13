// Description: Production-grade Sugar Removal component with comprehensive detection, removal, and extraction capabilities
import React, { useState } from "react";
import {
  HiOutlineSearch,
  HiOutlineTrash,
  HiOutlineBeaker,
  HiOutlineInformationCircle,
  HiOutlineExclamationCircle,
  HiChevronDown,
  HiChevronUp,
  HiOutlineRefresh,
  HiOutlineBookOpen,
  HiOutlineLightBulb,
} from "react-icons/hi";
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import LoadingScreen from "../common/LoadingScreen";
import SMILESDisplay from "../common/SMILESDisplay";
import {
  getSugarInfo,
  removeLinearSugars,
  removeCircularSugars,
  removeAllSugars,
  extractAglyconeAndSugars,
} from "../../services/toolsService";

// Default options configuration
const DEFAULT_OPTIONS = {
  // Basic options (apply to both sugar types)
  only_terminal: true,
  preservation_mode: 2,
  preservation_threshold: 5,
  mark_attach_points: false,

  // Circular sugar options
  gly_bond: false,
  oxygen_atoms: true,
  oxygen_atoms_threshold: 0.5,
  spiro_sugars: false,
  keto_sugars: false,

  // Linear sugar options
  linear_sugars_in_rings: false,
  linear_sugars_min_size: 4,
  linear_sugars_max_size: 7,
  linear_acidic_sugars: false,

  // Extract-specific options
  extract_circular_sugars: true,
  extract_linear_sugars: false,
  post_process_sugars: false,
  limit_post_process_by_size: false,
};

const SugarRemovalView = () => {
  // State management
  const [smiles, setSmiles] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);

  // Operation mode: 'detect', 'remove', 'extract'
  const [operationMode, setOperationMode] = useState("detect");

  // Remove sub-options: which sugar types to remove
  const [removeCircular, setRemoveCircular] = useState(true);
  const [removeLinear, setRemoveLinear] = useState(false);

  // Options state
  const [options, setOptions] = useState(DEFAULT_OPTIONS);

  // Collapsible sections state
  const [expandedSections, setExpandedSections] = useState({
    basic: true,
    circular: false,
    linear: false,
    extract: false,
  });

  // Results state
  const [results, setResults] = useState(null);

  // Sugar highlighting toggle for detect mode (circular/linear)
  const [highlightMode, setHighlightMode] = useState("circular");

  // Toggle section expansion
  const toggleSection = (section) => {
    setExpandedSections((prev) => ({
      ...prev,
      [section]: !prev[section],
    }));
  };

  // Reset all options to defaults
  const resetAllOptions = () => {
    setOptions(DEFAULT_OPTIONS);
  };

  // Reset individual option
  const resetOption = (optionKey) => {
    setOptions((prev) => ({
      ...prev,
      [optionKey]: DEFAULT_OPTIONS[optionKey],
    }));
  };

  // Update option value
  const updateOption = (key, value) => {
    setOptions((prev) => ({
      ...prev,
      [key]: value,
    }));
  };

  // Execute the selected operation
  const handleExecute = async () => {
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setResults(null);
      return;
    }

    setIsLoading(true);
    setError(null);
    setResults(null);

    try {
      let result;

      switch (operationMode) {
        case "detect":
          // First get detection info
          result = await getSugarInfo(trimmedSmiles, options);

          // Determine which sugars are present from the message
          const hasCircular = result.toLowerCase().includes("circular");
          const hasLinear = result.toLowerCase().includes("linear");

          let circularSugars = [];
          let linearSugars = [];

          // Extract sugars for highlighting if any are present
          if (hasCircular || hasLinear) {
            try {
              // Separate circular and linear sugars by extracting them individually
              if (hasCircular) {
                const circularOptions = {
                  ...options,
                  extract_circular_sugars: true,
                  extract_linear_sugars: false,
                };
                const circularResult = await extractAglyconeAndSugars(
                  trimmedSmiles,
                  circularOptions
                );
                if (
                  Array.isArray(circularResult) &&
                  circularResult.length > 1
                ) {
                  circularSugars = circularResult.slice(1);
                }
              }

              if (hasLinear) {
                const linearOptions = {
                  ...options,
                  extract_circular_sugars: false,
                  extract_linear_sugars: true,
                };
                const linearResult = await extractAglyconeAndSugars(
                  trimmedSmiles,
                  linearOptions
                );
                if (Array.isArray(linearResult) && linearResult.length > 1) {
                  linearSugars = linearResult.slice(1);
                }
              }
            } catch (extractError) {
              console.warn(
                "Could not extract sugars for highlighting:",
                extractError
              );
            }
          }

          setResults({
            type: "detect",
            message: result,
            originalSmiles: trimmedSmiles,
            hasCircular,
            hasLinear,
            circularSugars,
            linearSugars,
          });
          break;

        case "remove":
          // Determine which API to call based on checkbox selections
          if (removeCircular && removeLinear) {
            // Both selected - remove all sugars
            result = await removeAllSugars(trimmedSmiles, options);
          } else if (removeLinear) {
            // Only linear sugars
            result = await removeLinearSugars(trimmedSmiles, options);
          } else if (removeCircular) {
            // Only circular sugars
            result = await removeCircularSugars(trimmedSmiles, options);
          } else {
            // None selected
            setError("Please select at least one sugar type to remove");
            setIsLoading(false);
            return;
          }

          setResults({
            type: "remove",
            removeCircular: removeCircular,
            removeLinear: removeLinear,
            originalSmiles: trimmedSmiles,
            resultSmiles: result,
          });
          break;

        case "extract":
          result = await extractAglyconeAndSugars(trimmedSmiles, options);
          if (Array.isArray(result) && result.length > 0) {
            setResults({
              type: "extract",
              originalSmiles: trimmedSmiles,
              aglycone: result[0],
              sugars: result.slice(1),
            });
          } else {
            setError("No sugars were extracted from the molecule.");
          }
          break;

        default:
          setError("Invalid operation mode selected.");
      }
    } catch (err) {
      console.error(`Error during ${operationMode} operation:`, err);
      setError(`Error: ${err.message || "Unknown error occurred"}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Render individual option control
  const renderOption = (key, label, type, min, max, step, description) => {
    const value = options[key];
    const defaultValue = DEFAULT_OPTIONS[key];
    const isModified = value !== defaultValue;

    return (
      <div className="space-y-2">
        <div className="flex items-center justify-between">
          <label className="block text-sm font-medium text-gray-700 dark:text-gray-300">
            {label}
            {isModified && (
              <span className="ml-2 text-xs text-blue-600 dark:text-blue-400">
                (Modified)
              </span>
            )}
          </label>
          <button
            type="button"
            onClick={() => resetOption(key)}
            disabled={!isModified}
            className={`p-1 rounded hover:bg-gray-100 dark:hover:bg-gray-700 transition-colors ${
              isModified
                ? "text-blue-600 dark:text-blue-400"
                : "text-gray-400 dark:text-gray-600 cursor-not-allowed"
            }`}
            title="Reset to default"
          >
            <HiOutlineRefresh className="h-4 w-4" />
          </button>
        </div>

        {type === "boolean" ? (
          <div className="flex items-center">
            <input
              type="checkbox"
              checked={value}
              onChange={(e) => updateOption(key, e.target.checked)}
              className="h-4 w-4 text-blue-600 border-gray-300 dark:border-gray-600 rounded focus:ring-blue-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
            />
            <span className="ml-2 text-sm text-gray-600 dark:text-gray-400">
              {description}
            </span>
          </div>
        ) : type === "number" ? (
          <div className="space-y-1">
            <input
              type="number"
              value={value}
              onChange={(e) => updateOption(key, parseFloat(e.target.value))}
              min={min}
              max={max}
              step={step || 1}
              className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 dark:focus:ring-blue-500 dark:bg-gray-700 dark:text-white"
            />
            <p className="text-xs text-gray-500 dark:text-gray-400">
              {description}
            </p>
          </div>
        ) : type === "select" ? (
          <div className="space-y-1">
            <select
              value={value}
              onChange={(e) => updateOption(key, parseInt(e.target.value))}
              className="w-full px-3 py-2 border border-gray-300 dark:border-gray-600 rounded-md shadow-sm focus:outline-none focus:ring-2 focus:ring-blue-500 dark:focus:ring-blue-500 dark:bg-gray-700 dark:text-white"
            >
              <option value={1}>
                All - Preserve all disconnected structures
              </option>
              <option value={2}>
                Heavy Atom Count - Remove structures below threshold
              </option>
              <option value={3}>
                Molecular Weight - Remove structures below threshold
              </option>
            </select>
            <p className="text-xs text-gray-500 dark:text-gray-400">
              {description}
            </p>
          </div>
        ) : null}
      </div>
    );
  };

  // Render collapsible section
  const renderSection = (sectionKey, title, children) => {
    const isExpanded = expandedSections[sectionKey];

    return (
      <div className="border border-gray-200 dark:border-gray-700 rounded-lg overflow-hidden">
        <button
          type="button"
          onClick={() => toggleSection(sectionKey)}
          className="w-full px-4 py-3 bg-gray-50 dark:bg-gray-800 hover:bg-gray-100 dark:hover:bg-gray-750 flex items-center justify-between transition-colors"
        >
          <span className="text-sm font-semibold text-gray-700 dark:text-gray-300">
            {title}
          </span>
          {isExpanded ? (
            <HiChevronUp className="h-5 w-5 text-gray-500 dark:text-gray-400" />
          ) : (
            <HiChevronDown className="h-5 w-5 text-gray-500 dark:text-gray-400" />
          )}
        </button>

        {isExpanded && (
          <div className="px-4 py-4 space-y-4 bg-white dark:bg-gray-800">
            {children}
          </div>
        )}
      </div>
    );
  };

  return (
    <div className="space-y-6 p-4 md:p-6">
      {/* Header with Logo */}
      <div className="bg-gradient-to-br from-white via-blue-50 to-white dark:from-gray-800 dark:via-gray-800 dark:to-gray-900 p-8 rounded-xl shadow-xl dark:shadow-2xl relative overflow-hidden">
        {/* Decorative background elements */}
        <div className="absolute top-0 left-0 w-64 h-64 bg-blue-100 dark:bg-blue-900 rounded-full opacity-20 blur-3xl -ml-32 -mt-32"></div>
        <div className="absolute bottom-0 right-0 w-48 h-48 bg-purple-100 dark:bg-purple-900 rounded-full opacity-20 blur-3xl -mr-24 -mb-24"></div>

        <div className="flex items-center justify-between relative z-10">
          {/* Left side - Logo */}
          <div className="flex-shrink-0 mr-8">
            <div className="relative">
              {/* Logo with strong glow in dark mode */}
              <img
                src="/SDS.png"
                alt="Sugar Removal Logo"
                className="relative h-52 w-52 object-contain transform hover:scale-105 transition-transform duration-300"
                style={{
                  filter:
                    "drop-shadow(0 0 20px rgba(96, 165, 250, 0.8)) drop-shadow(0 0 40px rgba(56, 189, 248, 0.6)) drop-shadow(0 0 60px rgba(59, 130, 246, 0.4))",
                }}
              />
            </div>
          </div>

          {/* Right side - Reference */}
          <div className="flex-1">
            <div className="bg-white/60 dark:bg-gray-900/60 backdrop-blur-sm p-6 rounded-lg">
              <div className="flex items-start space-x-3">
                <div className="flex-shrink-0 w-1 h-16 bg-gradient-to-b from-blue-500 to-purple-500 rounded-full"></div>
                <div>
                  <h3 className="text-lg font-bold text-gray-800 dark:text-blue-300 mb-2">
                    <HiOutlineBookOpen
                      className="inline text-2xl mr-2 text-blue-600 dark:text-blue-400"
                      aria-hidden="true"
                    />
                    Reference
                  </h3>
                  <p className="text-sm text-gray-700 dark:text-gray-300 leading-relaxed">
                    <strong className="text-blue-600 dark:text-blue-400">
                      Schaub, J., Zielesny, A., Steinbeck, C., Sorokina, M.
                    </strong>
                    <br />
                    <em className="text-gray-600 dark:text-gray-400">
                      Too sweet: cheminformatics for deglycosylation in natural
                      products.
                    </em>
                    <br />
                    <span className="text-gray-500 dark:text-gray-500">
                      J Cheminform 12, 67 (2020).
                    </span>
                  </p>
                  <a
                    href="https://doi.org/10.1186/s13321-020-00467-y"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="inline-flex items-center mt-3 text-sm font-semibold text-blue-600 dark:text-blue-400 hover:text-blue-700 dark:hover:text-blue-300 transition-colors group"
                  >
                    <span>View Publication</span>
                    <svg
                      className="w-4 h-4 ml-1 group-hover:translate-x-1 transition-transform"
                      fill="none"
                      stroke="currentColor"
                      viewBox="0 0 24 24"
                    >
                      <path
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        strokeWidth={2}
                        d="M13 7l5 5m0 0l-5 5m5-5H6"
                      />
                    </svg>
                  </a>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* MORTAR Note */}
        <div className="mt-6 relative z-10">
          <div className="bg-gradient-to-r from-blue-50 via-indigo-50 to-purple-50 dark:from-blue-900/30 dark:via-indigo-900/30 dark:to-purple-900/30 border-l-4 border-blue-500 dark:border-blue-400 p-4 rounded-r-lg shadow-md backdrop-blur-sm">
            <div className="flex items-start">
              <div className="flex-shrink-0 mr-3">
                <HiOutlineInformationCircle className="h-6 w-6 text-blue-600 dark:text-blue-400 animate-pulse" />
              </div>
              <div className="text-sm">
                <p className="text-gray-800 dark:text-gray-200 font-medium mb-1">
                  <HiOutlineLightBulb className="inline h-5 w-5 mr-2 text-blue-600 dark:text-blue-400" />
                  <span className="text-blue-600 dark:text-blue-400 font-bold">
                    Batch Processing
                  </span>
                </p>
                <p className="text-gray-700 dark:text-gray-300">
                  For processing multiple molecules, use{" "}
                  <a
                    href="https://felixbaensch.github.io/MORTAR"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="font-bold text-blue-600 dark:text-blue-400 hover:text-blue-700 dark:hover:text-blue-300 underline decoration-2 underline-offset-2 hover:decoration-blue-500 transition-all inline-flex items-center"
                  >
                    MORTAR
                    <svg
                      className="w-3 h-3 ml-1"
                      fill="currentColor"
                      viewBox="0 0 20 20"
                    >
                      <path d="M11 3a1 1 0 100 2h2.586l-6.293 6.293a1 1 0 101.414 1.414L15 6.414V9a1 1 0 102 0V4a1 1 0 00-1-1h-5z"></path>
                      <path d="M5 5a2 2 0 00-2 2v8a2 2 0 002 2h8a2 2 0 002-2v-3a1 1 0 10-2 0v3H5V7h3a1 1 0 000-2H5z"></path>
                    </svg>
                  </a>{" "}
                  (MOlecule fragmenTation fRamework) - a comprehensive tool for
                  systematic molecule fragmentation and analysis, or the{" "}
                  <a
                    href="https://github.com/JonasSchaub/SugarRemoval"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="font-bold text-blue-600 dark:text-blue-400 hover:text-blue-700 dark:hover:text-blue-300 underline decoration-2 underline-offset-2 hover:decoration-blue-500 transition-all inline-flex items-center"
                  >
                    SugarRemoval CMD app
                    <svg
                      className="w-3 h-3 ml-1"
                      fill="currentColor"
                      viewBox="0 0 20 20"
                    >
                      <path d="M11 3a1 1 0 100 2h2.586l-6.293 6.293a1 1 0 101.414 1.414L15 6.414V9a1 1 0 102 0V4a1 1 0 00-1-1h-5z"></path>
                      <path d="M5 5a2 2 0 00-2 2v8a2 2 0 002 2h8a2 2 0 002-2v-3a1 1 0 10-2 0v3H5V7h3a1 1 0 000-2H5z"></path>
                    </svg>
                  </a>{" "}
                  for command-line batch processing.
                </p>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Input Section */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <h3 className="text-lg font-semibold text-gray-800 dark:text-gray-200 mb-4">
          Input Molecule
        </h3>

        <SMILESInput
          value={smiles}
          onChange={setSmiles}
          label="SMILES String"
          placeholder="Enter SMILES (e.g., C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)"
          required
        />
      </div>

      {/* Operation Mode Selection */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <h3 className="text-lg font-semibold text-gray-800 dark:text-gray-200 mb-4">
          Operation Mode
        </h3>

        <div className="space-y-4">
          {/* Mode selection radio buttons */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <label
              className={`flex items-center p-4 border-2 rounded-lg cursor-pointer transition-all ${
                operationMode === "detect"
                  ? "border-blue-500 bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20"
                  : "border-gray-300 dark:border-gray-600 hover:border-blue-300 dark:hover:border-blue-700"
              }`}
            >
              <input
                type="radio"
                value="detect"
                checked={operationMode === "detect"}
                onChange={(e) => setOperationMode(e.target.value)}
                className="h-4 w-4 text-blue-600 border-gray-300 focus:ring-blue-500"
              />
              <div className="ml-3">
                <HiOutlineSearch className="h-6 w-6 text-blue-600 dark:text-blue-400 mb-1" />
                <span className="block text-sm font-medium text-gray-900 dark:text-gray-100">
                  Detect
                </span>
                <span className="block text-xs text-gray-500 dark:text-gray-400">
                  Identify sugar moieties
                </span>
              </div>
            </label>

            <label
              className={`flex items-center p-4 border-2 rounded-lg cursor-pointer transition-all ${
                operationMode === "remove"
                  ? "border-red-500 bg-red-50 dark:bg-red-900 dark:bg-opacity-20"
                  : "border-gray-300 dark:border-gray-600 hover:border-red-300 dark:hover:border-red-700"
              }`}
            >
              <input
                type="radio"
                value="remove"
                checked={operationMode === "remove"}
                onChange={(e) => setOperationMode(e.target.value)}
                className="h-4 w-4 text-red-600 border-gray-300 focus:ring-red-500"
              />
              <div className="ml-3">
                <HiOutlineTrash className="h-6 w-6 text-red-600 dark:text-red-400 mb-1" />
                <span className="block text-sm font-medium text-gray-900 dark:text-gray-100">
                  Remove
                </span>
                <span className="block text-xs text-gray-500 dark:text-gray-400">
                  Remove sugar moieties
                </span>
              </div>
            </label>

            <label
              className={`flex items-center p-4 border-2 rounded-lg cursor-pointer transition-all ${
                operationMode === "extract"
                  ? "border-green-500 bg-green-50 dark:bg-green-900 dark:bg-opacity-20"
                  : "border-gray-300 dark:border-gray-600 hover:border-green-300 dark:hover:border-green-700"
              }`}
            >
              <input
                type="radio"
                value="extract"
                checked={operationMode === "extract"}
                onChange={(e) => setOperationMode(e.target.value)}
                className="h-4 w-4 text-green-600 border-gray-300 focus:ring-green-500"
              />
              <div className="ml-3">
                <HiOutlineBeaker className="h-6 w-6 text-green-600 dark:text-green-400 mb-1" />
                <span className="block text-sm font-medium text-gray-900 dark:text-gray-100">
                  Extract
                </span>
                <span className="block text-xs text-gray-500 dark:text-gray-400">
                  Extract aglycone & sugars
                </span>
              </div>
            </label>
          </div>

          {/* Remove sub-options */}
          {operationMode === "remove" && (
            <div className="pt-4 border-t border-gray-200 dark:border-gray-700">
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-3">
                Sugar Types to Remove
              </label>
              <div className="space-y-2">
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={removeCircular}
                    onChange={(e) => setRemoveCircular(e.target.checked)}
                    className="h-4 w-4 text-red-600 border-gray-300 rounded focus:ring-red-500 dark:border-gray-600 dark:bg-gray-700"
                  />
                  <span className="text-sm text-gray-700 dark:text-gray-300">
                    Remove Circular Sugars
                  </span>
                </label>
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={removeLinear}
                    onChange={(e) => setRemoveLinear(e.target.checked)}
                    className="h-4 w-4 text-red-600 border-gray-300 rounded focus:ring-red-500 dark:border-gray-600 dark:bg-gray-700"
                  />
                  <span className="text-sm text-gray-700 dark:text-gray-300">
                    Remove Linear Sugars
                  </span>
                </label>
              </div>
            </div>
          )}

          {/* Extract sub-options */}
          {operationMode === "extract" && (
            <div className="pt-4 border-t border-gray-200 dark:border-gray-700">
              <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-3">
                Sugar Types to Extract
              </label>
              <div className="space-y-2">
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={options.extract_circular_sugars}
                    onChange={(e) =>
                      updateOption("extract_circular_sugars", e.target.checked)
                    }
                    className="h-4 w-4 text-green-600 border-gray-300 rounded focus:ring-green-500 dark:border-gray-600 dark:bg-gray-700"
                  />
                  <span className="text-sm text-gray-700 dark:text-gray-300">
                    Extract Circular Sugars
                  </span>
                </label>
                <label className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={options.extract_linear_sugars}
                    onChange={(e) =>
                      updateOption("extract_linear_sugars", e.target.checked)
                    }
                    className="h-4 w-4 text-green-600 border-gray-300 rounded focus:ring-green-500 dark:border-gray-600 dark:bg-gray-700"
                  />
                  <span className="text-sm text-gray-700 dark:text-gray-300">
                    Extract Linear Sugars
                  </span>
                </label>
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Settings Section */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <div className="flex items-center justify-between mb-4">
          <h3 className="text-lg font-semibold text-gray-800 dark:text-gray-200">
            Detection & Processing Options
          </h3>
          <button
            type="button"
            onClick={resetAllOptions}
            className="flex items-center space-x-2 px-4 py-2 text-sm font-medium text-blue-600 dark:text-blue-400 hover:bg-blue-50 dark:hover:bg-blue-900 dark:hover:bg-opacity-20 rounded-lg transition-colors"
          >
            <HiOutlineRefresh className="h-4 w-4" />
            <span>Reset All to Default</span>
          </button>
        </div>

        <div className="space-y-4">
          {/* Basic Options Section */}
          {renderSection(
            "basic",
            "Basic Options (Both Sugar Types)",
            <>
              {renderOption(
                "only_terminal",
                "Remove Only Terminal Sugars",
                "boolean",
                null,
                null,
                null,
                "Only remove sugars at the terminal positions of the molecule"
              )}

              {renderOption(
                "preservation_mode",
                "Preservation Mode",
                "select",
                1,
                3,
                1,
                "Determines which disconnected structures to preserve after sugar removal"
              )}

              {renderOption(
                "preservation_threshold",
                "Preservation Threshold",
                "number",
                0,
                100,
                1,
                "Threshold value for the selected preservation mode (e.g., minimum heavy atoms)"
              )}

              {renderOption(
                "mark_attach_points",
                "Mark Attachment Points",
                "boolean",
                null,
                null,
                null,
                "Mark the attachment points of removed sugars with a dummy atom"
              )}
            </>
          )}

          {/* Circular Sugar Options Section */}
          {renderSection(
            "circular",
            "Circular Sugar Options",
            <>
              {renderOption(
                "gly_bond",
                "Detect Only with O-Glycosidic Bonds",
                "boolean",
                null,
                null,
                null,
                "Consider only circular sugars with glycosidic bonds"
              )}

              {renderOption(
                "oxygen_atoms",
                "Require Sufficient Exocyclic Oxygen Atoms",
                "boolean",
                null,
                null,
                null,
                "Consider only circular sugars with enough exocyclic oxygen atoms (see threshold below)"
              )}

              {renderOption(
                "oxygen_atoms_threshold",
                "Exocyclic Oxygen Atoms Ratio Threshold",
                "number",
                0.0,
                1.0,
                0.1,
                "Minimum ratio of exocyclic oxygen atoms to ring atoms (e.g., 0.5 means a 6-membered ring needs ≥3 exocyclic oxygens)"
              )}

              {renderOption(
                "spiro_sugars",
                "Include Spiro Sugars",
                "boolean",
                null,
                null,
                null,
                "Include rings that share one atom with another cycle (spiro rings)"
              )}

              {renderOption(
                "keto_sugars",
                "Detect Keto Sugars",
                "boolean",
                null,
                null,
                null,
                "Detect circular sugars with keto groups"
              )}
            </>
          )}

          {/* Linear Sugar Options Section */}
          {renderSection(
            "linear",
            "Linear Sugar Options",
            <>
              {renderOption(
                "linear_sugars_in_rings",
                "Detect Linear Sugars in Rings",
                "boolean",
                null,
                null,
                null,
                "Consider linear sugar patterns that are part of ring structures"
              )}

              {renderOption(
                "linear_sugars_min_size",
                "Linear Sugars Minimum Size",
                "number",
                0,
                20,
                1,
                "Minimum carbon chain length to be considered a linear sugar"
              )}

              {renderOption(
                "linear_sugars_max_size",
                "Linear Sugars Maximum Size",
                "number",
                1,
                20,
                1,
                "Maximum carbon chain length to be considered a linear sugar"
              )}

              {renderOption(
                "linear_acidic_sugars",
                "Detect Linear Acidic Sugars",
                "boolean",
                null,
                null,
                null,
                "Include linear sugar moieties with acidic functional groups"
              )}
            </>
          )}

          {/* Extract-specific options (only show when in extract mode) */}
          {operationMode === "extract" &&
            renderSection(
              "extract",
              "Advanced Extract Options",
              <>
                {renderOption(
                  "post_process_sugars",
                  "Post-process Extracted Sugars",
                  "boolean",
                  null,
                  null,
                  null,
                  "Split extracted sugars by breaking glycosidic, ether, ester, and peroxide bonds"
                )}

                {renderOption(
                  "limit_post_process_by_size",
                  "Limit Post-processing by Size",
                  "boolean",
                  null,
                  null,
                  null,
                  "Only post-process sugar structures larger than the preservation threshold"
                )}
              </>
            )}
        </div>
      </div>

      {/* Execute Button */}
      <div className="flex justify-center">
        <button
          onClick={handleExecute}
          disabled={!smiles.trim() || isLoading}
          className={`px-8 py-3 rounded-lg text-white font-semibold text-lg flex items-center justify-center space-x-3 transition-all shadow-lg ${
            !smiles.trim() || isLoading
              ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
              : operationMode === "detect"
              ? "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600"
              : operationMode === "remove"
              ? "bg-red-600 hover:bg-red-700 dark:bg-red-500 dark:hover:bg-red-600"
              : "bg-green-600 hover:bg-green-700 dark:bg-green-500 dark:hover:bg-green-600"
          }`}
        >
          {operationMode === "detect" && (
            <HiOutlineSearch className="h-6 w-6" />
          )}
          {operationMode === "remove" && <HiOutlineTrash className="h-6 w-6" />}
          {operationMode === "extract" && (
            <HiOutlineBeaker className="h-6 w-6" />
          )}
          <span>
            {isLoading
              ? "Processing..."
              : operationMode === "detect"
              ? "Detect Sugars"
              : operationMode === "remove"
              ? "Remove Sugars"
              : "Extract Aglycone & Sugars"}
          </span>
        </button>
      </div>

      {/* Loading State */}
      {isLoading && (
        <LoadingScreen
          text={`${
            operationMode === "detect"
              ? "Detecting"
              : operationMode === "remove"
              ? "Removing"
              : "Extracting"
          } sugar moieties...`}
        />
      )}

      {/* Error Display */}
      {error && !isLoading && (
        <div
          className="p-4 rounded-md flex items-start shadow bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700"
          role="alert"
        >
          <HiOutlineExclamationCircle
            className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
            aria-hidden="true"
          />
          <span>{error}</span>
        </div>
      )}

      {/* Results Display */}
      {results && !isLoading && !error && (
        <div className="space-y-6">
          {/* Detect Results */}
          {results.type === "detect" && (
            <div className="space-y-4">
              <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 shadow">
                <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-2">
                  Detection Result
                </h3>
                <p className="text-gray-700 dark:text-gray-200 font-medium">
                  {results.message}
                </p>
              </div>

              {/* Toggle for sugar highlighting - only show if both types are present */}
              {results.hasCircular && results.hasLinear && (
                <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow border border-gray-200 dark:border-gray-700">
                  <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-3">
                    Highlight Sugar Type
                  </label>
                  <div className="flex gap-3">
                    <button
                      onClick={() => setHighlightMode("circular")}
                      className={`flex-1 px-4 py-2 rounded-lg font-medium transition-all ${
                        highlightMode === "circular"
                          ? "bg-blue-600 text-white shadow-md"
                          : "bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-600"
                      }`}
                    >
                      Circular Sugars
                    </button>
                    <button
                      onClick={() => setHighlightMode("linear")}
                      className={`flex-1 px-4 py-2 rounded-lg font-medium transition-all ${
                        highlightMode === "linear"
                          ? "bg-blue-600 text-white shadow-md"
                          : "bg-gray-100 dark:bg-gray-700 text-gray-700 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-600"
                      }`}
                    >
                      Linear Sugars
                    </button>
                  </div>
                </div>
              )}

              <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
                <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-4">
                  Molecule with Highlighted Sugars
                </h4>
                <SMILESDisplay smiles={results.originalSmiles} label="SMILES" />
                <div className="mt-4">
                  <MoleculeCard
                    smiles={results.originalSmiles}
                    title="Analyzed Molecule"
                    description={
                      results.hasCircular && results.hasLinear
                        ? `Showing ${highlightMode} sugars highlighted`
                        : results.hasCircular
                        ? "Circular sugars highlighted"
                        : results.hasLinear
                        ? "Linear sugars highlighted"
                        : results.message
                    }
                    showActions={true}
                    substructures={
                      results.hasCircular && results.hasLinear
                        ? highlightMode === "circular"
                          ? results.circularSugars
                          : results.linearSugars
                        : results.hasCircular
                        ? results.circularSugars
                        : results.hasLinear
                        ? results.linearSugars
                        : []
                    }
                  />
                </div>
              </div>
            </div>
          )}

          {/* Remove Results */}
          {results.type === "remove" && (
            <div className="space-y-4">
              <h3 className="text-lg font-semibold text-gray-800 dark:text-gray-200">
                Sugar Removal Result (
                {results.removeCircular && results.removeLinear
                  ? "Circular & Linear Sugars"
                  : results.removeCircular
                  ? "Circular Sugars Only"
                  : "Linear Sugars Only"}
                )
              </h3>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
                  <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-3">
                    Original Molecule
                  </h4>
                  <SMILESDisplay
                    smiles={results.originalSmiles}
                    label="Original SMILES"
                  />
                  <div className="mt-4">
                    <MoleculeCard
                      smiles={results.originalSmiles}
                      title="Before Removal"
                      description="Original molecule with sugars"
                      showActions={false}
                    />
                  </div>
                </div>

                <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
                  <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-3">
                    Sugar-Free Molecule
                  </h4>
                  <SMILESDisplay
                    smiles={results.resultSmiles}
                    label="Aglycone SMILES"
                  />
                  <div className="mt-4">
                    <MoleculeCard
                      smiles={results.resultSmiles}
                      title="After Removal"
                      description={`${
                        results.removeCircular && results.removeLinear
                          ? "Circular & linear sugars"
                          : results.removeCircular
                          ? "Circular sugars"
                          : "Linear sugars"
                      } removed`}
                      showActions={false}
                    />
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* Extract Results */}
          {results.type === "extract" && (
            <div className="space-y-6">
              <h3 className="text-lg font-semibold text-gray-800 dark:text-gray-200">
                Extraction Results
              </h3>

              {/* Aglycone */}
              <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
                <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-3 flex items-center">
                  <span className="bg-green-100 dark:bg-green-900 text-green-800 dark:text-green-200 px-2 py-1 rounded text-xs font-semibold mr-2">
                    AGLYCONE
                  </span>
                  Core Structure (Sugar-Free)
                </h4>
                <SMILESDisplay
                  smiles={results.aglycone}
                  label="Aglycone SMILES"
                />
                <div className="mt-4">
                  <MoleculeCard
                    smiles={results.aglycone}
                    title="Aglycone"
                    description="Core structure after sugar removal"
                    showActions={true}
                  />
                </div>
              </div>

              {/* Extracted Sugars */}
              {results.sugars && results.sugars.length > 0 && (
                <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
                  <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-4">
                    Extracted Sugar Moieties ({results.sugars.length})
                  </h4>

                  <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                    {results.sugars.map((sugarSmiles, index) => (
                      <div
                        key={index}
                        className="border border-gray-200 dark:border-gray-700 rounded-lg p-4"
                      >
                        <div className="mb-2">
                          <span className="bg-blue-100 dark:bg-blue-900 text-blue-800 dark:text-blue-200 px-2 py-1 rounded text-xs font-semibold">
                            SUGAR {index + 1}
                          </span>
                        </div>
                        <SMILESDisplay
                          smiles={sugarSmiles}
                          label={`Sugar ${index + 1} SMILES`}
                        />
                        <div className="mt-3">
                          <MoleculeCard
                            smiles={sugarSmiles}
                            title={`Sugar ${index + 1}`}
                            description="Extracted sugar moiety"
                            showActions={true}
                            size="md"
                          />
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              )}

              {(!results.sugars || results.sugars.length === 0) && (
                <div className="p-4 rounded-md flex items-start shadow bg-yellow-50 dark:bg-yellow-900 dark:bg-opacity-20 text-yellow-700 dark:text-yellow-200 border border-yellow-300 dark:border-yellow-700">
                  <HiOutlineInformationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5" />
                  <span>
                    No sugar moieties were extracted from the molecule.
                  </span>
                </div>
              )}
            </div>
          )}
        </div>
      )}

      {/* Initial State / Help */}
      {!results && !isLoading && !error && (
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 flex items-start space-x-4 shadow">
          <HiOutlineInformationCircle className="h-6 w-6 text-blue-600 dark:text-blue-400 flex-shrink-0 mt-0.5" />
          <div>
            <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-2">
              Getting Started
            </h3>
            <div className="text-gray-700 dark:text-gray-300 space-y-2 text-sm">
              <p>
                <strong>1. Enter a SMILES string</strong> representing the
                molecule you want to analyze.
              </p>
              <p>
                <strong>2. Select an operation mode:</strong>
              </p>
              <ul className="list-disc list-inside ml-4 space-y-1">
                <li>
                  <strong>Detect:</strong> Identify if the molecule contains
                  linear and/or circular sugar moieties
                </li>
                <li>
                  <strong>Remove:</strong> Remove sugar moieties and obtain the
                  aglycone (core structure)
                </li>
                <li>
                  <strong>Extract:</strong> Separate the aglycone and individual
                  sugar moieties as distinct structures
                </li>
              </ul>
              <p>
                <strong>3. Configure detection options</strong> in the
                collapsible sections below (or use defaults).
              </p>
              <p>
                <strong>4. Click the Execute button</strong> to process your
                molecule.
              </p>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default SugarRemovalView;
