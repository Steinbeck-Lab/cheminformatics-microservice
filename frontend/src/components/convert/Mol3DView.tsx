// Description: This component allows users to input a SMILES string and generate 3D coordinates in MOL format using different cheminformatics toolkits. It includes error handling, loading states, and options for copying and downloading the generated data.
import React, { useState } from "react";
// Ensure all used icons are imported
import {
  HiOutlineCube,
  HiOutlineClipboard,
  HiOutlineDownload,
  HiOutlineCheck, // Added for copy success state
  HiOutlineExclamationCircle, // Added for error display
} from "react-icons/hi";
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import LoadingScreen from "../common/LoadingScreen";
// Assuming this service is configured correctly
import { generate3DCoordinates } from "../../services/convertService"; // Assuming this service exists

// Toolkit options configuration (adjust if different for 3D)
const TOOLKIT_OPTIONS_3D = [
  { id: "openbabel", label: "Open Babel" },
  { id: "rdkit", label: "RDKit" },
  // Add other toolkits your service supports for 3D generation
];

const Mol3DView = () => {
  const [smiles, setSmiles] = useState("");
  const [toolkit, setToolkit] = useState("openbabel"); // Default toolkit for 3D
  const [isLoading, setIsLoading] = useState(false);
  const [molblock, setMolblock] = useState(null); // Store the resulting 3D molblock string
  const [error, setError] = useState(null);
  const [copied, setCopied] = useState(false); // State for copy button feedback

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setMolblock(null); // Clear previous result
      return;
    }

    setIsLoading(true);
    setError(null);
    setMolblock(null); // Clear previous result before fetching

    try {
      // Call the service function for 3D coordinates
      const data = await generate3DCoordinates(trimmedSmiles, toolkit);
      // Ensure result is a non-empty string
      if (typeof data === "string" && data.trim()) {
        setMolblock(data);
      } else {
        // Handle cases where API might return empty or invalid data
        throw new Error("Received empty or invalid 3D Molblock data from the server.");
      }
    } catch (err) {
      console.error("3D Molblock generation error:", err); // Log the actual error
      setError(err.message || "Failed to generate 3D coordinates");
      setMolblock(null); // Ensure molblock is null on error
    } finally {
      setIsLoading(false);
    }
  };

  // Handle copying the molblock to clipboard
  const handleCopy = () => {
    if (!molblock || !navigator.clipboard) return;

    navigator.clipboard
      .writeText(molblock)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000); // Reset copied state
      })
      .catch((err) => {
        console.error("Failed to copy 3D molblock:", err);
        setError("Failed to copy 3D Molblock to clipboard."); // Show error to user
      });
  };

  // Handle downloading the molblock as a .mol file
  const downloadMolblock = () => {
    if (!molblock) return;

    try {
      const blob = new Blob([molblock], {
        type: "chemical/x-mdl-molfile;charset=utf-8",
      });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      // Create a filename from SMILES if possible, otherwise default
      const filenameBase = smiles
        ? smiles.replace(/[^a-z0-9]/gi, "_").substring(0, 30)
        : "molecule";
      a.download = `${filenameBase}_3d.mol`; // Indicate 3D in filename
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url); // Clean up object URL
    } catch (err) {
      console.error("Error creating download link:", err);
      setError("Could not create file for download.");
    }
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          3D Coordinate Generation (Molblock)
        </h2>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          {/* SMILES Input Component */}
          {/* Assuming SMILESInput handles its own dark/light mode */}
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            placeholder="Enter SMILES string (e.g., CN1C=NC2=C1C(=O)N(C(=O)N2C)C)"
            label="Input SMILES"
            required
          />

          {/* Toolkit Selection */}
          <div>
            <label
              htmlFor="toolkit-select-mol3d"
              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
            >
              Cheminformatics Toolkit
            </label>
            {/* Select styling */}
            <select
              id="toolkit-select-mol3d"
              value={toolkit}
              onChange={(e) => setToolkit(e.target.value)}
              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
            >
              {TOOLKIT_OPTIONS_3D.map(
                (
                  option // Use specific options for 3D
                ) => (
                  <option key={option.id} value={option.id}>
                    {option.label}
                  </option>
                )
              )}
            </select>
            {/* Hint text styling */}
            <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
              Different toolkits use different 3D conformation generation algorithms.
            </p>
          </div>

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              // Button styling
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !smiles.trim() || isLoading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed" // Disabled state
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm" // Enabled state
              }`}
              disabled={!smiles.trim() || isLoading}
            >
              <HiOutlineCube className="mr-2 h-5 w-5" aria-hidden="true" />
              {isLoading ? "Generating..." : "Generate 3D Coordinates"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {isLoading && <LoadingScreen text="Generating 3D coordinates..." />}

      {/* Error Display */}
      {error &&
        !isLoading && ( // Show error only if not loading
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

      {/* Results Display Section */}
      {/* Show only if molblock exists and not loading */}
      {molblock && !isLoading && (
        // Results card styling
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
          {/* Results Header with Action Buttons */}
          <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center mb-4 gap-3 border-b border-gray-200 dark:border-gray-700 pb-3">
            <h3 className="text-lg font-semibold text-gray-900 dark:text-blue-100">
              3D Molblock (MDL Molfile V2000)
            </h3>
            {/* Action Buttons */}
            <div className="flex space-x-2 flex-shrink-0">
              {/* Copy Button */}
              <button
                onClick={handleCopy}
                className={`px-3 py-1.5 text-sm rounded-md flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                  copied
                    ? "bg-green-100 dark:bg-green-700 text-green-700 dark:text-green-200" // Copied state
                    : "bg-gray-100 hover:bg-gray-200 dark:bg-gray-700 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600" // Default state
                }`}
                title="Copy Molblock to clipboard"
                aria-label={copied ? "Molblock Copied" : "Copy Molblock"}
              >
                {copied ? (
                  <HiOutlineCheck className="mr-1.5 h-4 w-4" /> // Use Check icon
                ) : (
                  <HiOutlineClipboard className="mr-1.5 h-4 w-4" />
                )}
                {copied ? "Copied!" : "Copy"}
              </button>
              {/* Download Button */}
              <button
                onClick={downloadMolblock}
                className="px-3 py-1.5 text-sm bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 text-white rounded-md flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500"
                title="Download as MOL file"
                aria-label="Download 3D Molblock as MOL file"
              >
                <HiOutlineDownload className="mr-1.5 h-4 w-4" />
                Download (.mol)
              </button>
            </div>
          </div>

          {/* Molblock Display Area */}
          <div className="bg-gray-100 dark:bg-gray-900 border border-gray-200 dark:border-gray-700 rounded-lg p-4 overflow-auto max-h-96 shadow-sm">
            {/* Molblock text styling */}
            <pre className="text-gray-700 dark:text-gray-300 font-mono text-xs sm:text-sm whitespace-pre">
              {molblock}
            </pre>
          </div>

          {/* Informational Note */}
          <div className="mt-4 text-xs text-gray-500 dark:text-gray-400">
            <p>
              This MOL file contains 3D coordinates (x, y, z) for each atom, representing a likely
              low-energy conformer. Generated using{" "}
              {TOOLKIT_OPTIONS_3D.find((opt) => opt.id === toolkit)?.label || toolkit}.
            </p>
          </div>
        </div>
      )}

      {/* Initial State Message */}
      {/* Show only if no molblock, not loading, and no error */}
      {!molblock && !isLoading && !error && (
        // Initial state card styling
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
          <div className="text-center text-gray-500 dark:text-gray-400 py-8">
            <HiOutlineCube className="h-12 w-12 mx-auto mb-3 text-gray-400 dark:text-gray-500" />
            <p>
              Enter a SMILES string and select a toolkit to generate 3D coordinates in MOL format.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default Mol3DView;
