// Description: This component allows users to input a SMILES string and calculate molecular descriptors using different toolkits (RDKit, CDK).
import React, { useState } from "react";
// Ensure all used icons are imported
import {
  HiOutlineDocumentSearch,
  HiOutlineClipboard,
  HiOutlineExclamationCircle, // Added for error display
} from "react-icons/hi";
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import LoadingScreen from "../common/LoadingScreen";
// Assuming this service is configured correctly
import { generateHOSECodes } from "../../services/chemService"; // Assuming this service exists

const HOSECodeView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [hoseCodes, setHoseCodes] = useState([]); // Initialize as empty array
  const [spheres, setSpheres] = useState(2); // Default spheres
  const [toolkit, setToolkit] = useState("rdkit"); // Default toolkit
  const [ringsize, setRingsize] = useState(false); // Default ringsize option
  const [copied, setCopied] = useState(false); // State for copy button feedback

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setHoseCodes([]); // Clear previous results
      return;
    }

    setLoading(true);
    setError(null);
    setHoseCodes([]); // Clear previous results before fetching

    try {
      // Call the service function
      const result = await generateHOSECodes(trimmedSmiles, spheres, toolkit, ringsize);
      // Ensure result is always an array
      setHoseCodes(Array.isArray(result) ? result : []);
    } catch (err) {
      console.error("HOSE code generation error:", err); // Log the error
      setError(`Error generating HOSE codes: ${err.message || "An unknown error occurred."}`);
      setHoseCodes([]); // Ensure codes are empty on error
    } finally {
      setLoading(false);
    }
  };

  // Copy all HOSE codes to clipboard and provide feedback
  const copyAllHoseCodes = () => {
    if (hoseCodes.length === 0 || !navigator.clipboard) return;

    const text = hoseCodes.join("\n"); // Join codes with newlines
    navigator.clipboard
      .writeText(text)
      .then(() => {
        setCopied(true); // Set copied state to true
        setTimeout(() => setCopied(false), 2000); // Reset after 2 seconds
      })
      .catch((err) => {
        console.error("Failed to copy HOSE codes:", err);
        setError("Failed to copy HOSE codes to clipboard."); // Show error to user
      });
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          HOSE Code Generator
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

          {/* Options Grid */}
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 pt-2">
            {/* Spheres Selection */}
            <div>
              <label
                htmlFor="spheres-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Number of Spheres
              </label>
              <select
                id="spheres-select"
                value={spheres}
                onChange={(e) => setSpheres(Number(e.target.value))}
                // Select styling for light/dark mode
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                <option value={1}>1</option>
                <option value={2}>2</option>
                <option value={3}>3</option>
                <option value={4}>4</option>
              </select>
            </div>

            {/* Toolkit Selection */}
            <div>
              <label
                htmlFor="toolkit-select-hose"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Toolkit
              </label>
              <select
                id="toolkit-select-hose"
                value={toolkit}
                onChange={(e) => setToolkit(e.target.value)}
                // Select styling for light/dark mode
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                <option value="rdkit">RDKit</option>
                <option value="cdk">CDK</option>
              </select>
            </div>

            {/* Ringsize Checkbox */}
            {/* Align checkbox with label vertically */}
            <div className="flex items-end pb-1">
              <div className="flex items-center">
                <input
                  id="ringsize-checkbox"
                  type="checkbox"
                  checked={ringsize}
                  onChange={(e) => setRingsize(e.target.checked)}
                  // Checkbox styling for light/dark mode
                  className="h-4 w-4 rounded bg-gray-50 dark:bg-gray-700 border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 focus:ring-blue-500 dark:focus:ring-offset-gray-800 shadow-sm"
                />
                <label
                  htmlFor="ringsize-checkbox"
                  className="ml-3 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Include ring size
                </label>
              </div>
            </div>
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
              <HiOutlineDocumentSearch className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? "Generating..." : "Generate HOSE Codes"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Generating HOSE codes..." />}

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
      {/* Show only if results exist and not loading and no error */}
      {hoseCodes.length > 0 && !loading && !error && (
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
          {/* Results Header */}
          <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center mb-4 gap-3 border-b border-gray-200 dark:border-gray-700 pb-3">
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white">
              HOSE Codes{" "}
              <span className="text-sm font-normal text-gray-500 dark:text-gray-400">
                ({hoseCodes.length} atoms)
              </span>
            </h3>
            {/* Copy Button */}
            <button
              onClick={copyAllHoseCodes}
              className={`px-3 py-1.5 text-sm rounded-md flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                copied
                  ? "bg-green-100 dark:bg-green-700 text-green-700 dark:text-green-200"
                  : "bg-gray-100 hover:bg-gray-200 dark:bg-gray-700 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600"
              }`}
              title="Copy all HOSE codes to clipboard"
            >
              <HiOutlineClipboard className="mr-1.5 h-4 w-4" aria-hidden="true" />
              {copied ? "Copied!" : "Copy All"}
            </button>
          </div>

          {/* Table Container */}
          {/* Added max-height and overflow for scrolling */}
          <div className="bg-gray-50 dark:bg-gray-900 rounded-lg p-4 border border-gray-200 dark:border-gray-700 shadow-sm overflow-y-auto max-h-96">
            <table className="w-full border-collapse">
              {/* Table Header */}
              <thead className="sticky top-0 bg-gray-100 dark:bg-gray-700 z-10">
                {/* Added sticky header */}
                <tr>
                  <th
                    scope="col"
                    className="px-4 py-2 text-left text-xs font-medium uppercase tracking-wider text-gray-500 dark:text-gray-400 border-b border-gray-200 dark:border-gray-600"
                  >
                    Atom Index
                  </th>
                  <th
                    scope="col"
                    className="px-4 py-2 text-left text-xs font-medium uppercase tracking-wider text-gray-500 dark:text-gray-400 border-b border-gray-200 dark:border-gray-600"
                  >
                    HOSE Code
                  </th>
                </tr>
              </thead>
              {/* Table Body */}
              <tbody className="divide-y divide-gray-100 dark:divide-gray-800">
                {hoseCodes.map((code, index) => (
                  <tr
                    key={index}
                    className="hover:bg-gray-100 dark:hover:bg-gray-800 transition-colors duration-150"
                  >
                    <td className="px-4 py-2 text-sm text-gray-600 dark:text-gray-300 whitespace-nowrap">
                      {index}
                    </td>
                    <td className="px-4 py-2 font-mono text-xs text-gray-700 dark:text-gray-300 break-all">
                      {code}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          {/* Informational Note */}
          <div className="mt-4 text-sm text-gray-600 dark:text-gray-400">
            <p>
              HOSE codes describe the chemical environment around each atom. Spheres indicate
              distance from the central atom.
            </p>
          </div>
        </div>
      )}

      {/* Initial State / About Box */}
      {/* Show only if no results, not loading, and no error */}
      {!hoseCodes.length && !loading && !error && (
        <div
          className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 shadow"
          role="complementary"
        >
          <h4 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-2">
            About HOSE Codes
          </h4>
          <p className="text-gray-700 dark:text-gray-300">
            HOSE (Hierarchically Ordered Spherical Environment) codes describe the chemical
            environment around each atom in a molecule. They are useful in NMR shift prediction,
            QSAR/QSPR studies, and structural elucidation.
          </p>
          <p className="mt-2 text-gray-700 dark:text-gray-300">
            The "spheres" parameter controls the depth of the environment considered (distance in
            bonds). More spheres yield more specific codes.
          </p>
        </div>
      )}
    </div>
  );
};

export default HOSECodeView;
