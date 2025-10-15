// Description: A reusable React component for uploading and converting MOL/SDF files to SMILES
import React, { useState, useRef } from "react";
import {
  HiOutlineUpload,
  HiOutlineDocumentText,
  HiOutlineX,
  HiOutlineCheckCircle,
  HiOutlineExclamationCircle,
} from "react-icons/hi";
import convertService from "../../services/convertService";

/**
 * MolFileUpload Component
 *
 * A reusable component that allows users to upload MOL/SDF files and convert them to SMILES using CDK or RDKit.
 *
 * @param {Object} props
 * @param {Function} props.onConversionSuccess - Callback when conversion is successful, receives (smiles, molblock, filename)
 * @param {Function} props.onConversionError - Callback when conversion fails, receives (error)
 * @param {string} props.toolkit - Toolkit to use for conversion: "cdk" (default) or "rdkit"
 * @param {string} props.className - Additional CSS classes for the container
 * @param {boolean} props.showMolblock - Whether to display the MOL block content (default: false)
 * @param {boolean} props.allowMultipleMolecules - Whether to allow SDF files with multiple molecules (default: false)
 */
const MolFileUpload = ({
  onConversionSuccess,
  onConversionError,
  toolkit = "cdk",
  className = "",
  showMolblock = false,
  allowMultipleMolecules = false,
}) => {
  const [molblock, setMolblock] = useState("");
  const [filename, setFilename] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [smiles, setSmiles] = useState("");
  const fileInputRef = useRef(null);

  // Helper function to format molblock for backend
  const formatMolblockForBackend = (molblock) => {
    // Normalize line endings first
    let formatted = molblock.replace(/\r\n/g, "\n").replace(/\r/g, "\n");

    // Clean up - remove trailing spaces from each line and trim
    formatted = formatted
      .split("\n")
      .map((line) => line.trimEnd())
      .join("\n")
      .trim();

    // Ensure proper molblock format - starts with newline
    if (!formatted.startsWith("\n")) {
      formatted = "\n" + formatted;
    }

    // Ensure it ends with M  END followed by a newline
    if (!formatted.endsWith("M  END\n")) {
      if (formatted.endsWith("M  END")) {
        formatted += "\n";
      }
    }

    return formatted;
  };

  // Handle file upload and read content
  const handleFileUpload = async (e) => {
    const file = e.target.files[0];
    if (!file) return;

    // Reset state
    setError(null);
    setSmiles("");
    setMolblock("");
    setFilename("");

    // Validate file extension
    const fileName = file.name.toLowerCase();
    if (!fileName.endsWith(".mol") && !fileName.endsWith(".sdf")) {
      const errorMsg = "Please upload a valid .mol or .sdf file.";
      setError(errorMsg);
      if (onConversionError) onConversionError(errorMsg);
      return;
    }

    setFilename(file.name);
    setLoading(true);

    const reader = new FileReader();

    reader.onload = async (event) => {
      try {
        let fileContent = event.target.result;

        // Normalize line endings and trim
        fileContent = fileContent.replace(/\r\n/g, "\n").replace(/\r/g, "\n").trim();

        // Validate molblock format
        if (!fileContent.includes("M  END")) {
          throw new Error(
            'Invalid file format: The uploaded file does not appear to be a valid MOL/SDF file (missing "M  END").'
          );
        }

        // Check for multiple molecules if not allowed
        if (!allowMultipleMolecules) {
          const molBlockCount = (fileContent.match(/M {2}END/g) || []).length;
          if (molBlockCount > 1) {
            throw new Error(
              "Multiple molecules detected in file. This component only supports single molecules. Please upload a file containing only one molecule."
            );
          }
        }

        // Handle SDF terminator - extract first molecule
        if (fileContent.includes("$$$$")) {
          fileContent = fileContent.split("$$$$")[0].trim();
        }

        setMolblock(fileContent);

        // Format and convert to SMILES
        const formattedMolblock = formatMolblockForBackend(fileContent);
        const convertedSmiles = await convertService.molblockToSMILES(formattedMolblock, toolkit);

        // Clean up SMILES (remove quotes if present)
        const cleanSmiles = convertedSmiles.replace(/^"|"$/g, "").trim();

        setSmiles(cleanSmiles);
        setError(null);

        // Call success callback
        if (onConversionSuccess) {
          onConversionSuccess(cleanSmiles, fileContent, file.name);
        }
      } catch (err) {
        console.error("File conversion error:", err);
        const errorMsg = err.message || "Failed to convert MOL/SDF file to SMILES.";
        setError(errorMsg);
        if (onConversionError) onConversionError(errorMsg);
      } finally {
        setLoading(false);
      }
    };

    reader.onerror = () => {
      const errorMsg = "Failed to read the uploaded file. Please ensure it's a valid text file.";
      setError(errorMsg);
      setLoading(false);
      if (onConversionError) onConversionError(errorMsg);
    };

    reader.readAsText(file);
  };

  // Clear/reset the component
  const handleClear = () => {
    setMolblock("");
    setFilename("");
    setSmiles("");
    setError(null);
    if (fileInputRef.current) {
      fileInputRef.current.value = "";
    }
  };

  return (
    <div className={`space-y-4 ${className}`}>
      {/* File Upload Button */}
      <div>
        <label
          htmlFor="mol-file-upload"
          className="flex items-center justify-center px-4 py-3 border-2 border-dashed border-gray-300 dark:border-gray-600 rounded-lg cursor-pointer hover:border-blue-500 dark:hover:border-blue-400 transition-colors duration-200 bg-gray-50 dark:bg-gray-800/50"
        >
          <input
            ref={fileInputRef}
            id="mol-file-upload"
            type="file"
            accept=".mol,.sdf"
            onChange={handleFileUpload}
            className="hidden"
            disabled={loading}
          />
          <HiOutlineUpload className="h-5 w-5 text-gray-500 dark:text-gray-400 mr-2" />
          <span className="text-sm font-medium text-gray-700 dark:text-gray-300">
            {filename || "Upload MOL/SDF File"}
          </span>
        </label>
        <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
          Supports .mol and .sdf file formats
        </p>
      </div>

      {/* Loading Indicator */}
      {loading && (
        <div className="flex items-center justify-center py-4">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600 dark:border-blue-400"></div>
          <span className="ml-3 text-sm text-gray-600 dark:text-gray-400">
            Converting to SMILES...
          </span>
        </div>
      )}

      {/* Error Display */}
      {error && !loading && (
        <div className="flex items-start p-4 bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-lg">
          <HiOutlineExclamationCircle className="h-5 w-5 text-red-600 dark:text-red-400 mt-0.5 flex-shrink-0" />
          <div className="ml-3 flex-1">
            <p className="text-sm text-red-800 dark:text-red-200">{error}</p>
          </div>
          <button
            onClick={() => setError(null)}
            className="ml-2 text-red-500 hover:text-red-700 dark:text-red-400 dark:hover:text-red-300"
          >
            <HiOutlineX className="h-4 w-4" />
          </button>
        </div>
      )}

      {/* Success Display with SMILES */}
      {smiles && !loading && !error && (
        <div className="space-y-3">
          <div className="flex items-start p-4 bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800 rounded-lg">
            <HiOutlineCheckCircle className="h-5 w-5 text-green-600 dark:text-green-400 mt-0.5 flex-shrink-0" />
            <div className="ml-3 flex-1">
              <p className="text-sm font-medium text-green-800 dark:text-green-200 mb-1">
                Conversion Successful
              </p>
              <div className="bg-white dark:bg-gray-800 p-3 rounded border border-green-200 dark:border-green-700">
                <p className="text-xs text-gray-500 dark:text-gray-400 mb-1">SMILES:</p>
                <p className="text-sm font-mono text-gray-900 dark:text-gray-100 break-all">
                  {smiles}
                </p>
              </div>
            </div>
            <button
              onClick={handleClear}
              className="ml-2 text-green-500 hover:text-green-700 dark:text-green-400 dark:hover:text-green-300"
              title="Clear and upload new file"
            >
              <HiOutlineX className="h-4 w-4" />
            </button>
          </div>

          {/* Optional MOL Block Display */}
          {showMolblock && molblock && (
            <div className="bg-gray-50 dark:bg-gray-800/50 p-4 rounded-lg border border-gray-200 dark:border-gray-700">
              <div className="flex items-center mb-2">
                <HiOutlineDocumentText className="h-4 w-4 text-gray-500 dark:text-gray-400 mr-2" />
                <p className="text-xs font-medium text-gray-600 dark:text-gray-400">
                  MOL Block Content
                </p>
              </div>
              <pre className="text-xs font-mono text-gray-700 dark:text-gray-300 overflow-x-auto max-h-60 overflow-y-auto bg-white dark:bg-gray-900 p-3 rounded border border-gray-200 dark:border-gray-700">
                {molblock}
              </pre>
            </div>
          )}
        </div>
      )}

      {/* Filename Display */}
      {filename && !loading && (
        <div className="flex items-center text-sm text-gray-600 dark:text-gray-400">
          <HiOutlineDocumentText className="h-4 w-4 mr-2" />
          <span className="truncate">{filename}</span>
        </div>
      )}
    </div>
  );
};

export default MolFileUpload;
