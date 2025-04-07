// Description: StandardizeView component for molecule standardization
import React, { useState } from 'react';
// Ensure all used icons are imported
import {
  HiOutlineDocumentReport,
  HiOutlineClipboard,
  HiOutlineUpload,
  HiOutlineExclamationCircle // Added for error display
} from 'react-icons/hi';
// Assuming these components are correctly implemented and styled for dark/light mode
// import SMILESInput from '../common/SMILESInput'; // Removed unused import
import LoadingScreen from '../common/LoadingScreen';
import MoleculeCard from '../common/MoleculeCard';
// Assuming the API URL is configured correctly via environment variables or defaults
const API_BASE_URL = process.env.REACT_APP_API_URL || 'https://api.naturalproducts.net';

const StandardizeView = () => {
  const [molblock, setMolblock] = useState('');
  const [standardizedData, setStandardizedData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [copiedStates, setCopiedStates] = useState({}); // State for copy button feedback (multiple buttons)

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedMolblock = molblock.trim();
    if (!trimmedMolblock) {
      setError('Please paste or upload a molblock.');
      setStandardizedData(null); // Clear previous results
      return;
    }

    setLoading(true);
    setError(null);
    setStandardizedData(null); // Clear previous results before fetching

    try {
      const response = await fetch(`${API_BASE_URL}/chem/standardize`, {
        method: 'POST',
        headers: {
          'Content-Type': 'text/plain', // Sending molblock as plain text
          'Accept': 'application/json' // Expecting JSON response
        },
        body: trimmedMolblock // Send trimmed molblock
      });

      if (!response.ok) {
        // Try to get error details from response body
        let errorMsg = `Error ${response.status}: ${response.statusText}`;
        try {
          const errorData = await response.json();
          errorMsg = errorData.detail || errorMsg; // Use detail if available
        } catch (jsonError) {
          // Ignore if response is not JSON
        }
        throw new Error(errorMsg);
      }

      const data = await response.json();
      setStandardizedData(data);
    } catch (err) {
      console.error("Standardization error:", err); // Log the error
      setError(`Error standardizing molecule: ${err.message || 'An unknown error occurred.'}`);
      setStandardizedData(null); // Ensure data is null on error
    } finally {
      setLoading(false);
    }
  };

  // Handle copy action for different text types
  const handleCopy = (text, type) => {
    if (!navigator.clipboard) {
      setError('Clipboard API not available in this browser.');
      return;
    }
    navigator.clipboard.writeText(text)
      .then(() => {
        setCopiedStates(prev => ({ ...prev, [type]: true })); // Set specific copied state
        setTimeout(() => setCopiedStates(prev => ({ ...prev, [type]: false })), 2000); // Reset after 2 seconds
      })
      .catch(err => {
        console.error(`Failed to copy ${type}:`, err);
        setError(`Failed to copy ${type} to clipboard.`);
      });
  };

  // Handle file upload and read content into textarea
  const handleFileUpload = (e) => {
    const file = e.target.files[0];
    if (!file) return;

    // Basic check for file type (can be enhanced)
    if (!file.name.toLowerCase().endsWith('.mol') && !file.name.toLowerCase().endsWith('.sdf')) {
      setError('Please upload a valid .mol or .sdf file.');
      return;
    }

    const reader = new FileReader();
    reader.onload = (event) => {
      setMolblock(event.target.result);
      setError(null); // Clear error on successful load
    };
    reader.onerror = (event) => {
      console.error("File reading error:", event.target.error);
      setError("Failed to read the uploaded file.");
    };
    reader.readAsText(file);
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">Standardize Molecule (via Molblock)</h2>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          {/* Molblock Textarea */}
          <div>
            <label htmlFor="molblock-input" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
              Molblock (MDL MOL format)
            </label>
            <textarea
              id="molblock-input"
              value={molblock}
              onChange={(e) => setMolblock(e.target.value)}
              placeholder="Paste your molblock here, or upload a .mol/.sdf file below..."
              rows={10}
              // Textarea styling for light/dark mode
              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white font-mono text-xs sm:text-sm resize-y shadow-sm focus:ring-indigo-500 focus:border-indigo-500"
              required
              aria-required="true"
            />
            {/* File Upload and Info */}
            <div className="mt-2 flex flex-col sm:flex-row justify-between items-start sm:items-center gap-2">
              {/* File Upload Button */}
              <label className="inline-flex items-center px-4 py-2 bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-800 dark:text-gray-200 text-sm font-medium rounded-md cursor-pointer border border-gray-300 dark:border-gray-600 shadow-sm transition-colors duration-150">
                <HiOutlineUpload className="mr-2 h-5 w-5" aria-hidden="true" />
                <span>Upload MOL/SDF File</span>
                <input
                  type="file"
                  accept=".mol,.sdf,chemical/x-mdl-molfile,chemical/x-mdl-sdfile" // MIME types
                  onChange={handleFileUpload}
                  className="sr-only" // Hide the default file input
                />
              </label>
              {/* Info Text */}
              <div className="text-xs text-gray-500 dark:text-gray-400 text-right">
                Standardization uses the ChEMBL curation pipeline.
              </div>
            </div>
          </div>

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              disabled={!molblock.trim() || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${!molblock.trim() || loading
                ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm'
                }`}
            >
              <HiOutlineDocumentReport className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? 'Standardizing...' : 'Standardize Molecule'}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Standardizing molecule..." />}

      {/* Error Display */}
      {error && !loading && ( // Show error only if not loading
        <div className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow" role="alert">
          <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
          <span>{error}</span>
        </div>
      )}

      {/* Results Display Section */}
      {/* Show only if data exists and not loading */}
      {standardizedData && !loading && (
        <div className="space-y-6">
          {/* Main Results Card */}
          <div className="bg-white dark:bg-gray-800 rounded-lg overflow-hidden shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
            {/* Results Header */}
            <div className="p-4 bg-gray-50 dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700">
              <h3 className="text-lg font-semibold text-gray-900 dark:text-white">Standardization Results</h3>
            </div>

            {/* Grid for Molecule Card and Identifiers */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6 p-6">
              {/* Molecule Card */}
              {/* Ensure MoleculeCard can handle SMILES and is theme-aware */}
              <MoleculeCard
                smiles={standardizedData.canonical_smiles} // Use canonical SMILES for display
                title="Standardized Structure"
                showActions={true} // Assuming MoleculeCard handles actions like copy/download if true
                size="md"
              />

              {/* Identifiers List */}
              <div className="space-y-4">
                {/* Canonical SMILES */}
                <div>
                  <label className="block text-sm font-medium text-gray-500 dark:text-gray-400 mb-1">Canonical SMILES</label>
                  <div className="flex bg-gray-100 dark:bg-gray-900 rounded-lg p-2 items-center border border-gray-200 dark:border-gray-700 shadow-sm">
                    <div className="flex-1 font-mono text-xs sm:text-sm text-gray-800 dark:text-blue-300 break-all pr-2">
                      {standardizedData.canonical_smiles}
                    </div>
                    <button
                      onClick={() => handleCopy(standardizedData.canonical_smiles, 'smiles')}
                      className="ml-auto p-1 text-gray-400 hover:text-gray-700 dark:hover:text-white flex-shrink-0"
                      title="Copy SMILES"
                    >
                      <HiOutlineClipboard className="h-5 w-5" />
                    </button>
                  </div>
                </div>
                {/* InChI */}
                <div>
                  <label className="block text-sm font-medium text-gray-500 dark:text-gray-400 mb-1">InChI</label>
                  <div className="flex bg-gray-100 dark:bg-gray-900 rounded-lg p-2 items-center border border-gray-200 dark:border-gray-700 shadow-sm">
                    <div className="flex-1 font-mono text-xs sm:text-sm text-gray-800 dark:text-blue-300 break-all pr-2">
                      {standardizedData.inchi}
                    </div>
                    <button
                      onClick={() => handleCopy(standardizedData.inchi, 'inchi')}
                      className="ml-auto p-1 text-gray-400 hover:text-gray-700 dark:hover:text-white flex-shrink-0"
                      title="Copy InChI"
                    >
                      <HiOutlineClipboard className="h-5 w-5" />
                    </button>
                  </div>
                </div>
                {/* InChI Key */}
                <div>
                  <label className="block text-sm font-medium text-gray-500 dark:text-gray-400 mb-1">InChI Key</label>
                  <div className="flex bg-gray-100 dark:bg-gray-900 rounded-lg p-2 items-center border border-gray-200 dark:border-gray-700 shadow-sm">
                    <div className="flex-1 font-mono text-xs sm:text-sm text-gray-800 dark:text-blue-300 break-all pr-2">
                      {standardizedData.inchikey}
                    </div>
                    <button
                      onClick={() => handleCopy(standardizedData.inchikey, 'inchikey')}
                      className="ml-auto p-1 text-gray-400 hover:text-gray-700 dark:hover:text-white flex-shrink-0"
                      title="Copy InChI Key"
                    >
                      <HiOutlineClipboard className="h-5 w-5" />
                    </button>
                  </div>
                </div>
              </div>
            </div>

            {/* Standardized Molblock Section */}
            <div className="p-6 border-t border-gray-200 dark:border-gray-700">
              <h4 className="text-sm font-medium text-gray-500 dark:text-gray-400 mb-2">Standardized Molblock</h4>
              {/* Molblock Display Area */}
              <div className="bg-gray-100 dark:bg-gray-900 rounded-lg p-3 max-h-60 overflow-auto border border-gray-200 dark:border-gray-700 shadow-sm">
                <pre className="font-mono text-xs text-gray-700 dark:text-gray-300 whitespace-pre">
                  {standardizedData.standardized_mol}
                </pre>
              </div>
              {/* Copy Molblock Button */}
              <button
                onClick={() => handleCopy(standardizedData.standardized_mol, 'molblock')}
                className={`mt-3 px-3 py-1 text-sm rounded-md flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${copiedStates['molblock']
                  ? 'bg-green-100 dark:bg-green-700 text-green-700 dark:text-green-200'
                  : 'bg-gray-100 hover:bg-gray-200 dark:bg-gray-700 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600'
                  }`}
              >
                <HiOutlineClipboard className="mr-1.5 h-4 w-4" aria-hidden="true" />
                {copiedStates['molblock'] ? "Copied!" : "Copy Molblock"}
              </button>
            </div>
          </div>

          {/* Additional Properties Section (Optional) */}
          {/* Filter out known keys and check if any others remain */}
          {Object.keys(standardizedData).filter(key => !['standardized_mol', 'canonical_smiles', 'inchi', 'inchikey'].includes(key)).length > 0 && (
            <div className="bg-white dark:bg-gray-800 rounded-lg p-6 shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
              <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-4">Additional Properties</h3>
              {/* Grid for additional properties */}
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {Object.entries(standardizedData)
                  .filter(([key]) => !['standardized_mol', 'canonical_smiles', 'inchi', 'inchikey'].includes(key))
                  .map(([key, value]) => (
                    // Property Card
                    <div key={key} className="bg-gray-50 dark:bg-gray-900 rounded-lg p-3 border border-gray-200 dark:border-gray-700 shadow-sm">
                      <h4 className="text-sm font-medium text-gray-500 dark:text-gray-400 mb-1 capitalize">
                        {key.replace(/_/g, ' ')} {/* Format key */}
                      </h4>
                      <div className="font-mono text-sm text-gray-800 dark:text-blue-300 break-all">
                        {/* Display value appropriately */}
                        {typeof value === 'object' ? JSON.stringify(value) : String(value)}
                      </div>
                    </div>
                  ))}
              </div>
            </div>
          )}
        </div>
      )}

      {/* Initial State Message (Optional) */}
      {!standardizedData && !loading && !error && (
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 text-center shadow">
          <p className="text-gray-600 dark:text-gray-300">
            Paste a Molblock or upload a MOL/SDF file to standardize the molecule using the ChEMBL curation pipeline.
          </p>
        </div>
      )}

    </div>
  );
};

export default StandardizeView;
