// Description: StandardizeView component for molecule standardization
import React, { useState } from 'react';
// Ensure all used icons are imported
import {
  HiOutlineDocumentReport,
  HiOutlineClipboard,
  HiOutlineUpload,
  HiOutlineExclamationCircle, // Added for error display
  HiOutlineBeaker,
  HiOutlineDocumentText
} from 'react-icons/hi';
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from '../common/SMILESInput';
import LoadingScreen from '../common/LoadingScreen';
import MoleculeCard from '../common/MoleculeCard';
import { generate2DCoordinates } from '../../services/convertService';
// Assuming the API URL is configured correctly via environment variables or defaults
const API_BASE_URL = process.env.REACT_APP_API_URL || 'https://dev.api.naturalproducts.net/latest';

const StandardizeView = () => {
  const [inputMethod, setInputMethod] = useState('molblock'); // 'molblock' or 'smiles'
  const [molblock, setMolblock] = useState('');
  const [smiles, setSmiles] = useState('');
  const [standardizedData, setStandardizedData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [copiedStates, setCopiedStates] = useState({}); // State for copy button feedback (multiple buttons)

  // Helper function to ensure molblock is in proper format for backend
  const formatMolblockForBackend = (molblock) => {
    // Normalize line endings first
    let formatted = molblock
      .replace(/\r\n/g, '\n')
      .replace(/\r/g, '\n');
    
    // Clean up - remove trailing spaces from each line and trim
    formatted = formatted
      .split('\n')
      .map(line => line.trimEnd())
      .join('\n')
      .trim();
    
    // Ensure proper molblock format - starts with newline (like the test fixture)
    if (!formatted.startsWith('\n')) {
      formatted = '\n' + formatted;
    }
    
    // Ensure it ends with M  END followed by a newline
    if (!formatted.endsWith('M  END\n')) {
      if (formatted.endsWith('M  END')) {
        formatted += '\n';
      }
    }
    
    return formatted;
  };

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    
    let molblockToStandardize = '';
    
    if (inputMethod === 'smiles') {
      const trimmedSmiles = smiles.trim();
      if (!trimmedSmiles) {
        setError('Please enter a SMILES string.');
        setStandardizedData(null);
        return;
      }
      
      // Check if SMILES contains multiple molecules
      if (trimmedSmiles.includes('.') && trimmedSmiles.split('.').filter(s => s.trim().length > 0).length > 1) {
        setError('Please enter only one molecule at a time. Multiple molecules (separated by dots) are not supported.');
        setStandardizedData(null);
        return;
      }
      
      try {
        setLoading(true);
        setError(null);
        // Convert SMILES to molblock first
        const rawMolblock = await generate2DCoordinates(trimmedSmiles, 'rdkit');
        molblockToStandardize = formatMolblockForBackend(rawMolblock);
      } catch (conversionError) {
        console.error("SMILES to molblock conversion error:", conversionError);
        setError(`Error converting SMILES to molblock: ${conversionError.message || 'Invalid SMILES string.'}`);
        setStandardizedData(null);
        setLoading(false);
        return;
      }
    } else {
      const trimmedMolblock = molblock.trim();
      if (!trimmedMolblock) {
        setError('Please paste or upload a molblock.');
        setStandardizedData(null);
        return;
      }
      
      // Validate molblock format more thoroughly
      if (!trimmedMolblock.includes('M  END')) {
        setError('Invalid molblock format: Missing "M  END" terminator. Please ensure you have pasted a complete molblock.');
        setStandardizedData(null);
        return;
      }
      
      // Check for multiple molecules in molblock
      const molBlockCount = (trimmedMolblock.match(/M {2}END/g) || []).length;
      if (molBlockCount > 1) {
        setError('Please provide only one molecule at a time. Multiple molecules in a single input are not supported.');
        setStandardizedData(null);
        return;
      }
      
      // Normalize line endings and ensure proper formatting
      const normalizedMolblock = formatMolblockForBackend(trimmedMolblock);
      molblockToStandardize = normalizedMolblock;
      
      setLoading(true);
      setError(null);
    }

    setStandardizedData(null); // Clear previous results before fetching

    // Debug log to see what we're sending
    console.log('Sending molblock to standardize:', {
      inputMethod,
      molblockLength: molblockToStandardize.length,
      molblockPreview: molblockToStandardize.substring(0, 200) + '...',
      endsWithNewline: molblockToStandardize.endsWith('\n'),
      containsMEnd: molblockToStandardize.includes('M  END')
    });

    try {
      const response = await fetch(`${API_BASE_URL}/chem/standardize`, {
        method: 'POST',
        headers: {
          'Content-Type': 'text/plain; charset=utf-8', // Be explicit about charset
          'Accept': 'application/json' // Expecting JSON response
        },
        body: molblockToStandardize // Send molblock as UTF-8 text
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
      // Provide more specific error messages
      let errorMessage = err.message || 'An unknown error occurred.';
      
      if (errorMessage.includes('cannot access local variable')) {
        errorMessage = 'Error parsing molblock: The molblock format appears to be invalid or corrupted. Please ensure you have pasted a complete, valid molblock with proper structure and exactly one molecule.';
      } else if (errorMessage.includes('422') && inputMethod === 'molblock') {
        errorMessage = 'Invalid molblock format: Please ensure your molblock is in proper MDL MOL format and contains exactly one valid molecule structure.';
      } else if (errorMessage.includes('Unable to parse molblock') || errorMessage.includes('Failure parsing molblock')) {
        errorMessage = 'Unable to parse molblock: Please check that your molblock is in valid MDL MOL format and try again.';
      } else if (errorMessage.includes('500') || errorMessage.includes('Internal Server Error')) {
        errorMessage = 'Server processing error: There was an issue processing your molblock. Please verify the format is correct and contains a single valid molecule.';
      }
      
      setError(errorMessage);
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
      let fileContent = event.target.result;
      
      // Normalize line endings and trim
      fileContent = fileContent
        .replace(/\r\n/g, '\n')  // Convert Windows line endings
        .replace(/\r/g, '\n')    // Convert Mac line endings
        .trim();
      
      // Basic validation for molblock format
      if (!fileContent.includes('M  END')) {
        setError('Invalid file format: The uploaded file does not appear to be a valid molblock (missing "M  END").');
        return;
      }
      
      // Check for multiple molecules
      const molBlockCount = (fileContent.match(/M {2}END/g) || []).length;
      if (molBlockCount > 1) {
        setError('Multiple molecules detected in file. Please upload a file containing only one molecule.');
        return;
      }
      
      setMolblock(fileContent);
      setError(null); // Clear error on successful load
    };
    reader.onerror = (event) => {
      console.error("File reading error:", event.target.error);
      setError("Failed to read the uploaded file. Please ensure it's a valid text file.");
    };
    reader.readAsText(file);
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">Standardize Molecule</h2>
        
        {/* Input Method Selection */}
        <div className="mb-4">
          <div className="flex rounded-lg bg-gray-100 dark:bg-gray-700 p-1">
            <button
              type="button"
              onClick={() => setInputMethod('smiles')}
              className={`flex-1 px-3 py-2 text-sm font-medium rounded-md transition-colors duration-200 flex items-center justify-center ${
                inputMethod === 'smiles'
                  ? 'bg-blue-600 text-white shadow-sm'
                  : 'text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-300'
              }`}
            >
              <HiOutlineBeaker className="mr-2 h-4 w-4" />
              SMILES Input
            </button>
            <button
              type="button"
              onClick={() => setInputMethod('molblock')}
              className={`flex-1 px-3 py-2 text-sm font-medium rounded-md transition-colors duration-200 flex items-center justify-center ${
                inputMethod === 'molblock'
                  ? 'bg-blue-600 text-white shadow-sm'
                  : 'text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-300'
              }`}
            >
              <HiOutlineDocumentText className="mr-2 h-4 w-4" />
              Molblock Input
            </button>
          </div>
        </div>

        {/* Important Notice */}
        <div className="mb-4 p-3 bg-yellow-50 dark:bg-yellow-900 dark:bg-opacity-20 border border-yellow-200 dark:border-yellow-800 rounded-md">
          <div className="flex">
            <HiOutlineExclamationCircle className="h-5 w-5 text-yellow-400 mr-2 flex-shrink-0 mt-0.5" />
            <div className="text-sm text-yellow-800 dark:text-yellow-200">
              <strong>Important:</strong> This tool accepts only one molecule at a time. Multiple molecules (separated by dots in SMILES or multiple records in molblock) will cause errors.
            </div>
          </div>
        </div>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          {/* Conditional Input Fields */}
          {inputMethod === 'smiles' ? (
            /* SMILES Input */
            <div>
              <SMILESInput
                value={smiles}
                onChange={setSmiles}
                placeholder="Enter SMILES string (e.g., CN1C=NC2=C1C(=O)N(C(=O)N2C)C)"
                label="SMILES String"
                required
              />
              <p className="mt-2 text-xs text-gray-500 dark:text-gray-400">
                Enter a single SMILES string. The system will convert it to molblock format for standardization.
              </p>
            </div>
          ) : (
            /* Molblock Input */
            <div>
              <label htmlFor="molblock-input" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                Molblock (MDL MOL format)
              </label>
              <textarea
                id="molblock-input"
                value={molblock}
                onChange={(e) => setMolblock(e.target.value)}
                placeholder={`Paste your molblock here, or upload a .mol/.sdf file below...

Example format (simple molecule):

  CDK     08302311362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.9743    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    1.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    2.8125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6238    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  END`}
                rows={12}
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
                  Must be in MDL MOL format with "M  END" terminator.
                </div>
              </div>
              {/* Format Hint */}
              <p className="mt-2 text-xs text-gray-500 dark:text-gray-400">
                Ensure the molblock is in proper MDL MOL format and ends with "M  END". The system expects exactly one complete molecule structure.
              </p>
            </div>
          )}

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              disabled={(inputMethod === 'smiles' ? !smiles.trim() : !molblock.trim()) || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${(inputMethod === 'smiles' ? !smiles.trim() : !molblock.trim()) || loading
                ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm'
                }`}
            >
              <HiOutlineDocumentReport className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? 'Standardizing...' : 'Standardize Molecule'}
            </button>
          </div>
          
          {/* Additional Info */}
          <div className="text-xs text-gray-500 dark:text-gray-400 mt-2">
            Standardization uses the ChEMBL curation pipeline to ensure consistent molecular representation.
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
          <div className="space-y-3">
            <HiOutlineBeaker className="h-12 w-12 mx-auto text-blue-500 dark:text-blue-400" />
            <h3 className="text-lg font-medium text-gray-900 dark:text-white">Molecule Standardization</h3>
            <div className="text-gray-600 dark:text-gray-300 space-y-2">
              <p>
                Use either SMILES or Molblock input to standardize your molecule using the ChEMBL curation pipeline.
              </p>
              <p className="text-sm">
                <strong>Note:</strong> Only one molecule per submission is supported. Multiple molecules will result in an error.
              </p>
            </div>
          </div>
        </div>
      )}

    </div>
  );
};

export default StandardizeView;
