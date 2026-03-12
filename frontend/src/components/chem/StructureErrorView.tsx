// Description: This component provides a user interface for validating and standardizing chemical structures using SMILES notation. It includes input handling, error checking, and displays results with appropriate messaging and styling for both light and dark themes.
import React, { useState } from "react";
// Ensure all used icons are imported
import {
  HiOutlineCheck,
  HiOutlineShieldCheck,
  HiOutlinePencilAlt,
  HiOutlineExclamationCircle, // Added for error display
  HiOutlineInformationCircle, // Added for about section
} from "react-icons/hi";
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import LoadingScreen from "../common/LoadingScreen";
// Assuming this service is configured correctly
import { checkStructureErrors } from "../../services/chemService"; // Assuming this service exists

const StructureErrorView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null); // Stores { original: { smi, messages }, standardized: { smi, messages } }
  const [fix, setFix] = useState(true); // Option to attempt fixing issues

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setResult(null); // Clear previous results
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null); // Clear previous results before fetching

    try {
      // Call the service function
      const validationResult = await checkStructureErrors(trimmedSmiles, fix);
      setResult(validationResult);

      // Optional: Check if standardized structure is same as original even if messages exist
      // if (fix && validationResult?.standardized?.smi === validationResult?.original?.smi) {
      //    // Could set an info message here if desired
      // }
    } catch (err) {
      console.error("Structure check error:", err); // Log the error
      setError(`Error checking structure: ${err.message || "An unknown error occurred."}`);
      setResult(null); // Ensure result is null on error
    } finally {
      setLoading(false);
    }
  };

  // Render validation messages with appropriate styling
  const renderMessages = (messages) => {
    // Handle null or empty messages array
    if (!messages || messages.length === 0) {
      // Optionally return a default "No messages" state or null
      return (
        <div className="flex items-center text-gray-500 dark:text-gray-400 text-sm p-3 bg-gray-50 dark:bg-gray-900 border border-gray-200 dark:border-gray-700 rounded-lg shadow-sm">
          <HiOutlineInformationCircle className="mr-2 h-5 w-5 flex-shrink-0" />
          No validation messages reported.
        </div>
      );
    }

    // Handle the specific "No Errors Found" message
    if (messages.length === 1 && messages[0] === "No Errors Found") {
      return (
        // Success message styling
        <div className="flex items-center text-green-700 dark:text-green-400 font-medium p-3 bg-green-50 dark:bg-green-900 dark:bg-opacity-20 border border-green-300 dark:border-green-700 rounded-lg shadow-sm">
          <HiOutlineCheck className="mr-2 h-5 w-5 flex-shrink-0" aria-hidden="true" />
          {messages[0]}
        </div>
      );
    }

    // Handle actual error/warning messages
    return (
      // Warning/Error message styling (using amber for visibility)
      <div className="p-3 bg-amber-50 dark:bg-amber-900 dark:bg-opacity-20 border border-amber-300 dark:border-amber-700 rounded-lg shadow-sm">
        <h4 className="font-medium text-amber-700 dark:text-amber-300 mb-2 flex items-center">
          <HiOutlineExclamationCircle className="h-5 w-5 mr-2 flex-shrink-0" aria-hidden="true" />
          Structure Issues Found
        </h4>
        <ul className="list-disc list-inside space-y-1 text-sm text-amber-700 dark:text-amber-200 pl-5">
          {messages.map((msg, index) => (
            <li key={index}>{msg}</li>
          ))}
        </ul>
      </div>
    );
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Structure Validation & Standardization
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

          {/* Fix Option Checkbox */}
          <div className="flex items-center pt-1">
            <input
              id="fix-structure"
              type="checkbox"
              checked={fix}
              onChange={(e) => setFix(e.target.checked)}
              // Checkbox styling for light/dark mode
              className="h-4 w-4 rounded bg-gray-50 dark:bg-gray-700 border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 focus:ring-blue-500 dark:focus:ring-offset-gray-800 shadow-sm"
            />
            <label
              htmlFor="fix-structure"
              className="ml-3 block text-sm text-gray-700 dark:text-gray-300"
            >
              Attempt to fix issues by standardizing (using ChEMBL pipeline)
            </label>
          </div>

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
              <HiOutlineShieldCheck className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? "Checking..." : "Check Structure"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Checking chemical structure..." />}

      {/* Error Display */}
      {error &&
        !loading && ( // Show error only if not loading
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

      {/* Results Display Section */}
      {/* Show only if result object exists and not loading */}
      {result && !loading && (
        <div className="space-y-6">
          {" "}
          {/* Increased spacing */}
          {/* Original Structure Section */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-4">
              Original Structure
            </h3>
            <div className="space-y-4">
              {/* Original SMILES Display */}
              <div className="flex items-center space-x-3 text-gray-700 dark:text-gray-300">
                <HiOutlinePencilAlt
                  className="h-5 w-5 text-blue-600 dark:text-blue-400 flex-shrink-0"
                  aria-hidden="true"
                />
                <span className="font-mono text-sm bg-gray-100 dark:bg-gray-900 px-3 py-1 rounded border border-gray-200 dark:border-gray-700 break-all">
                  {result.original?.smi || smiles} {/* Show original SMILES from result or input */}
                </span>
              </div>
              {/* Original Structure Messages */}
              {renderMessages(result.original?.messages)}
            </div>
          </div>
          {/* Standardized Structure Section (if fix was true and data exists) */}
          {fix && result.standardized && (
            <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
              <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-4">
                Standardized Structure
              </h3>
              <div className="space-y-4">
                {/* Standardized SMILES Display */}
                <div className="flex items-center space-x-3 text-gray-700 dark:text-gray-300">
                  {/* Use check icon if standardized messages indicate no errors */}
                  {result.standardized.messages?.includes("No Errors Found") ? (
                    <HiOutlineCheck
                      className="h-5 w-5 text-green-600 dark:text-green-400 flex-shrink-0"
                      aria-hidden="true"
                    />
                  ) : (
                    <HiOutlinePencilAlt
                      className="h-5 w-5 text-green-600 dark:text-green-400 flex-shrink-0"
                      aria-hidden="true"
                    />
                  )}
                  <span className="font-mono text-sm bg-gray-100 dark:bg-gray-900 px-3 py-1 rounded border border-gray-200 dark:border-gray-700 break-all">
                    {result.standardized.smi}
                  </span>
                </div>
                {/* Standardized Structure Messages */}
                {renderMessages(result.standardized.messages)}

                {/* Molecule Card Comparison (only if standardization occurred) */}
                {/* Ensure MoleculeCard is theme-aware */}
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4 pt-4 border-t border-gray-200 dark:border-gray-700 mt-4">
                  <MoleculeCard
                    smiles={result.original?.smi || smiles}
                    title="Original"
                    size="sm" // Smaller cards for comparison
                  />
                  <MoleculeCard
                    smiles={result.standardized.smi}
                    title="Standardized / Fixed"
                    size="sm"
                  />
                </div>
              </div>
            </div>
          )}
        </div>
      )}

      {/* Initial State / About Box */}
      {/* Show only if no result, not loading, and no error */}
      {!result && !loading && !error && (
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
              About Structure Validation
            </h3>
            <p className="text-gray-700 dark:text-gray-300">
              This tool checks chemical structures for potential issues (e.g., incorrect valences,
              invalid aromaticity) using the ChEMBL structure curation pipeline.
            </p>
            <p className="mt-2 text-gray-700 dark:text-gray-300">
              Enable the "Attempt to fix" option to apply standardization rules (normalization,
              aromaticity perception, etc.) to resolve common problems and generate a standardized
              representation.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default StructureErrorView;
