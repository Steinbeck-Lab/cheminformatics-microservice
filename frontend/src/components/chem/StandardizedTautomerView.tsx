// Description: This component allows users to input a SMILES string and calculate molecular descriptors using a selected toolkit (RDKit or CDK). It handles loading states, error messages, and displays results in a table format. The component is styled for both light and dark modes.
import React, { useState } from "react";
// Ensure all used icons are imported
import {
  HiOutlineRefresh,
  HiOutlineInformationCircle,
  HiOutlineExclamationCircle, // Added for error display
} from "react-icons/hi";
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import LoadingScreen from "../common/LoadingScreen";
// Assuming this service is configured correctly
import { generateStandardizedTautomer } from "../../services/chemService"; // Assuming this service exists

const StandardizedTautomerView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [standardizedSmiles, setStandardizedSmiles] = useState(null); // Store the result SMILES

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setStandardizedSmiles(null); // Clear previous results
      return;
    }

    setLoading(true);
    setError(null);
    setStandardizedSmiles(null); // Clear previous results before fetching

    try {
      // Call the service function
      const result = await generateStandardizedTautomer(trimmedSmiles);
      setStandardizedSmiles(result); // Store the result

      // Set an informational message instead of error if no different tautomer found
      if (!result || result === trimmedSmiles) {
        setError("Input is already the canonical tautomer or no different tautomers were found.");
        setStandardizedSmiles(trimmedSmiles); // Still show the original if it's canonical
      }
    } catch (err) {
      console.error("Tautomer standardization error:", err); // Log the error
      setError(`Error standardizing tautomer: ${err.message || "An unknown error occurred."}`);
      setStandardizedSmiles(null); // Ensure result is null on error
    } finally {
      setLoading(false);
    }
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Standardized Tautomer Generator
        </h2>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          {/* SMILES Input */}
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            required
            // Pass theme props if needed
          />

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              disabled={!smiles.trim() || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !smiles.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm"
              }`}
            >
              <HiOutlineRefresh className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? "Standardizing..." : "Generate Standardized Tautomer"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Standardizing tautomer..." />}

      {/* Error/Info Display */}
      {/* Show error/info only if not loading */}
      {error && !loading && (
        <div
          className={`p-4 rounded-md flex items-start shadow ${
            error.startsWith("Input is already") // Check if it's the informational message
              ? "bg-blue-50 dark:bg-blue-900 dark:bg-opacity-30 text-blue-700 dark:text-blue-200 border border-blue-300 dark:border-blue-700"
              : "bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700"
          }`}
          role={error.startsWith("Input is already") ? "status" : "alert"}
        >
          {error.startsWith("Input is already") ? (
            <HiOutlineInformationCircle
              className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-blue-500 dark:text-blue-400"
              aria-hidden="true"
            />
          ) : (
            <HiOutlineExclamationCircle
              className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
              aria-hidden="true"
            />
          )}
          <span>{error}</span>
        </div>
      )}

      {/* Results Display Section */}
      {/* Show only if standardizedSmiles exists and not loading */}
      {/* We show results even if error is set (for the 'already canonical' case) */}
      {standardizedSmiles && !loading && (
        <div className="space-y-4">
          {/* Main Results Card */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-4 border-b border-gray-200 dark:border-gray-700 pb-3">
              Tautomer Results
            </h3>

            {/* Grid for Original vs Standardized */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              {/* Original Structure Column */}
              <div>
                <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-2">
                  Original Structure
                </h4>
                {/* Ensure MoleculeCard is theme-aware */}
                <MoleculeCard
                  smiles={smiles} // Show original input SMILES
                  title="Original"
                  // description="Input molecule" // Optional description
                  size="md" // Adjust size as needed
                />
                {/* Display original SMILES string */}
                <div className="mt-2 bg-gray-100 dark:bg-gray-900 p-2 rounded-md border border-gray-200 dark:border-gray-700 font-mono text-xs overflow-x-auto shadow-sm">
                  <pre className="text-gray-700 dark:text-gray-300 whitespace-pre-wrap break-all">
                    {smiles}
                  </pre>
                </div>
              </div>

              {/* Standardized Structure Column */}
              <div>
                <h4 className="text-md font-medium text-blue-700 dark:text-blue-300 mb-2">
                  Standardized Tautomer
                </h4>
                {/* Ensure MoleculeCard is theme-aware */}
                <MoleculeCard
                  smiles={standardizedSmiles} // Show standardized SMILES
                  title="Standardized"
                  // description="Canonical tautomer" // Optional description
                  size="md"
                />
                {/* Display standardized SMILES string */}
                <div className="mt-2 bg-gray-100 dark:bg-gray-900 p-2 rounded-md border border-gray-200 dark:border-gray-700 font-mono text-xs overflow-x-auto shadow-sm">
                  <pre className="text-gray-700 dark:text-gray-300 whitespace-pre-wrap break-all">
                    {standardizedSmiles}
                  </pre>
                </div>
              </div>
            </div>
          </div>

          {/* Informational Box about Tautomers */}
          <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 shadow">
            <h4 className="text-md font-medium text-blue-800 dark:text-blue-300 mb-2">
              About Tautomers
            </h4>
            <p className="text-sm text-gray-700 dark:text-gray-300">
              Tautomers are constitutional isomers that readily interconvert, typically involving
              proton migration and a shift in double bonds (e.g., keto-enol).
            </p>
            <p className="text-sm text-gray-700 dark:text-gray-300 mt-2">
              The standardized tautomer is the canonical form chosen by RDKit's TautomerEnumerator
              for consistent representation.
            </p>
          </div>
        </div>
      )}

      {/* Initial State / About Box */}
      {/* Show only if no results, not loading, and no error */}
      {!standardizedSmiles && !loading && !error && (
        <div
          className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 flex items-start space-x-4 shadow"
          role="complementary"
        >
          <HiOutlineInformationCircle
            className="h-6 w-6 text-blue-600 dark:text-blue-400 flex-shrink-0 mt-0.5"
            aria-hidden="true"
          />
          <div>
            <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-2">
              About Tautomer Standardization
            </h3>
            <p className="text-gray-700 dark:text-gray-300">
              Tautomers are structural isomers (like keto-enol forms) that interconvert. This tool
              finds a consistent, canonical tautomeric form using RDKit.
            </p>
            {/* Example Section */}
            <div className="mt-4 grid grid-cols-1 sm:grid-cols-2 gap-4">
              {/* Keto Example */}
              <div className="bg-gray-100 dark:bg-gray-900 p-3 rounded-lg border border-gray-200 dark:border-gray-700 shadow-sm">
                <p className="text-gray-500 dark:text-gray-400 text-sm mb-1">
                  Example (Keto form):
                </p>
                <pre className="text-xs text-gray-700 dark:text-gray-300 font-mono">CC(=O)CC</pre>
              </div>
              {/* Enol Example */}
              <div className="bg-gray-100 dark:bg-gray-900 p-3 rounded-lg border border-gray-200 dark:border-gray-700 shadow-sm">
                <p className="text-gray-500 dark:text-gray-400 text-sm mb-1">
                  Example (Enol form):
                </p>
                <pre className="text-xs text-gray-700 dark:text-gray-300 font-mono">CC(O)=CC</pre>
              </div>
            </div>
            <p className="mt-3 text-xs text-gray-500 dark:text-gray-400">
              Both input forms typically yield the same standardized tautomer.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default StandardizedTautomerView;
