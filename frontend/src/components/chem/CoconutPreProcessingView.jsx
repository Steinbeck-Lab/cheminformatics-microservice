// Description: This component handles the COCONUT pre-processing view, including SMILES input, processing options, and displaying results.
import React, { useState } from 'react'; // Removed useEffect, useCallback
// Ensure all used icons are imported, including HiOutlineExclamationCircle for errors
import {
  HiOutlineDocumentReport,
  HiOutlineBeaker,
  HiOutlineInformationCircle,
  HiOutlineExclamationCircle // Added for error display
} from 'react-icons/hi';
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from '../common/SMILESInput';
import LoadingScreen from '../common/LoadingScreen';
// Assuming this service is configured correctly
import { coconutPreprocessing } from '../../services/chemService'; // Assuming this service exists

const CoconutPreProcessingView = () => {
  const [smiles, setSmiles] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [preprocessData, setPreprocessData] = useState(null);
  const [generate3D, setGenerate3D] = useState(false);
  const [generateDescriptors, setGenerateDescriptors] = useState(false);
  const [activeTab, setActiveTab] = useState('original'); // Default tab

  // Handle form submission to preprocess the molecule
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError('Please enter a SMILES string');
      return;
    }

    setLoading(true);
    setError(null);
    setPreprocessData(null); // Clear previous results
    setActiveTab('original'); // Reset to default tab on new submission

    try {
      // Call the service function (ensure it handles errors appropriately)
      const data = await coconutPreprocessing(trimmedSmiles, generate3D, generateDescriptors);
      setPreprocessData(data);
    } catch (err) {
      console.error("Preprocessing error:", err); // Log the error
      // Set a user-friendly error message
      setError(`Error processing data: ${err.message || 'An unknown error occurred.'}`);
    } finally {
      setLoading(false);
    }
  };

  // Helper function to format keys (e.g., 'some_key' -> 'Some Key')
  const formatKey = (key) => {
    if (typeof key !== 'string') return key;
    return key
      .replace(/_/g, ' ') // Replace underscores with spaces
      .replace(/\b\w/g, l => l.toUpperCase()); // Capitalize first letter of each word
  };

  // Helper function to format values for display
  const formatValue = (value) => {
    if (typeof value === 'boolean') {
      return value ? 'Yes' : 'No';
    }
    if (typeof value === 'number') {
      // Format numbers to a reasonable precision
      return parseFloat(value.toFixed(4)).toString(); // Avoid trailing zeros like 1.2300
    }
    // Return other types as is (or handle strings specifically if needed)
    return value;
  };

  // Renders the content for the active tab
  const renderTabContent = (tabData) => {
    if (!tabData) return <p className="text-gray-500 dark:text-gray-400 p-4">No data available for this tab.</p>;

    // Separate molblocks from other representations for distinct display
    const molblocks = {};
    const otherReps = {};

    if (tabData.representations) {
      Object.entries(tabData.representations).forEach(([key, value]) => {
        // Assuming keys for molblocks contain 'MOL'
        if (key.includes('MOL')) {
          molblocks[key] = value;
        } else {
          otherReps[key] = value;
        }
      });
    }

    return (
      <div className="space-y-6 py-4">
        {/* Molblocks Section */}
        {Object.keys(molblocks).length > 0 && (
          <div className="bg-gray-50 dark:bg-gray-900 rounded-lg p-4 border border-gray-200 dark:border-gray-700 shadow-sm">
            <h4 className="text-lg font-medium text-gray-800 dark:text-blue-300 mb-4">Molecular Structures</h4>
            <div className="space-y-4">
              {Object.entries(molblocks).map(([key, value]) => (
                <div key={key} className="space-y-1">
                  <div className="text-sm font-medium text-gray-500 dark:text-gray-400">{formatKey(key)}</div>
                  {/* Preformatted block for MOL data */}
                  <div className="bg-gray-100 dark:bg-gray-800 p-3 rounded overflow-x-auto max-h-80 border border-gray-200 dark:border-gray-700">
                    <pre className="text-gray-700 dark:text-gray-300 text-xs sm:text-sm font-mono whitespace-pre">{value}</pre>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Representations and Properties Grid */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Other Representations Card */}
          <div className="bg-gray-50 dark:bg-gray-900 rounded-lg p-4 border border-gray-200 dark:border-gray-700 shadow-sm">
            <h4 className="text-lg font-medium text-gray-800 dark:text-blue-300 mb-3">Representations</h4>
            <div className="space-y-3">
              {Object.keys(otherReps).length > 0 ? Object.entries(otherReps).map(([key, value]) => (
                <div key={key} className="space-y-1">
                  <div className="text-sm font-medium text-gray-500 dark:text-gray-400">{formatKey(key)}</div>
                  {/* Preformatted block for SMILES, InChI, etc. */}
                  <div className="bg-gray-100 dark:bg-gray-800 p-2 rounded overflow-x-auto border border-gray-200 dark:border-gray-700">
                    <pre className="text-gray-700 dark:text-gray-300 text-xs sm:text-sm font-mono break-words whitespace-pre-wrap">{value}</pre>
                  </div>
                </div>
              )) : <p className="text-sm text-gray-500 dark:text-gray-400">No other representations available.</p>}
            </div>
          </div>

          {/* Properties and Errors Column */}
          <div className="space-y-6">
            {/* Properties Card */}
            <div className="bg-gray-50 dark:bg-gray-900 rounded-lg p-4 border border-gray-200 dark:border-gray-700 shadow-sm">
              <h4 className="text-lg font-medium text-gray-800 dark:text-blue-300 mb-3">Properties</h4>
              <div className="space-y-2 text-sm">
                <div className="flex items-center justify-between">
                  <span className="text-gray-600 dark:text-gray-400">Has Stereo</span>
                  {/* Conditional styling for boolean properties */}
                  <span className={`font-medium ${tabData.has_stereo ? 'text-green-600 dark:text-green-400' : 'text-red-600 dark:text-red-400'}`}>
                    {tabData.has_stereo ? 'Yes' : 'No'}
                  </span>
                </div>
                {/* Only show 'Has Stereo Defined' if the property exists */}
                {tabData.has_stereo_defined !== undefined && (
                  <div className="flex items-center justify-between">
                    <span className="text-gray-600 dark:text-gray-400">Has Stereo Defined</span>
                    <span className={`font-medium ${tabData.has_stereo_defined ? 'text-green-600 dark:text-green-400' : 'text-red-600 dark:text-red-400'}`}>
                      {tabData.has_stereo_defined ? 'Yes' : 'No'}
                    </span>
                  </div>
                )}
                {/* Add more properties here if available in tabData */}
              </div>
            </div>

            {/* Validation Errors Card (only if errors exist) */}
            {tabData.errors && (
              <div className="bg-red-50 dark:bg-red-900 dark:bg-opacity-20 rounded-lg p-4 border border-red-200 dark:border-red-700 shadow-sm">
                <h4 className="text-lg font-medium text-red-700 dark:text-red-300 mb-3">Validation Errors</h4>
                {Array.isArray(tabData.errors) && tabData.errors.length > 0 ? (
                  <ul className="list-disc list-inside space-y-1 text-sm text-red-600 dark:text-red-400">
                    {tabData.errors.map((error, idx) => (
                      <li key={idx}>{error}</li>
                    ))}
                  </ul>
                ) : (
                  // This case might not happen if the card only renders when errors exist,
                  // but included for completeness.
                  <p className="text-sm text-green-600 dark:text-green-400">No errors found for this structure.</p>
                )}
              </div>
            )}
          </div>
        </div>

        {/* Descriptors Section (only if requested and available) */}
        {generateDescriptors && tabData.descriptors && (
          <div className="bg-gray-50 dark:bg-gray-900 rounded-lg p-4 border border-gray-200 dark:border-gray-700 shadow-sm mt-6">
            <h4 className="text-lg font-medium text-gray-800 dark:text-blue-300 mb-3">Descriptors</h4>
            {typeof tabData.descriptors === 'object' && Object.keys(tabData.descriptors).length > 0 ? (
              // Grid layout for descriptors
              <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-x-6 gap-y-2 text-sm">
                {Object.entries(tabData.descriptors).map(([key, value]) => (
                  <div key={key} className="flex justify-between border-b border-gray-100 dark:border-gray-700 py-1">
                    <span className="text-gray-600 dark:text-gray-400 mr-2">{formatKey(key)}</span>
                    <span className="text-gray-800 dark:text-gray-300 font-mono text-right">{formatValue(value)}</span>
                  </div>
                ))}
              </div>
            ) : (
              <p className="text-sm text-gray-500 dark:text-gray-400">No descriptors were calculated or available.</p>
            )}
          </div>
        )}
      </div>
    );
  };

  // Renders the main results section including tabs
  const renderResults = (data) => {
    if (!data) return null;

    return (
      <div className="space-y-6">
        {/* Tab Navigation */}
        <div className="border-b border-gray-200 dark:border-gray-700">
          <nav className="-mb-px flex space-x-6" aria-label="Tabs">
            {/* Map over available keys in preprocessData to create tabs */}
            {Object.keys(data).map((key) => (
              <button
                key={key}
                onClick={() => setActiveTab(key)}
                className={`whitespace-nowrap py-3 px-1 border-b-2 font-medium text-sm focus:outline-none ${activeTab === key
                    ? 'border-blue-500 dark:border-blue-400 text-blue-600 dark:text-blue-400'
                    : 'border-transparent text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-200 hover:border-gray-300 dark:hover:border-gray-600'
                  }`}
                aria-current={activeTab === key ? 'page' : undefined}
              >
                {formatKey(key)} {/* Format tab titles */}
              </button>
            ))}
          </nav>
        </div>

        {/* Render content for the active tab */}
        {renderTabContent(data[activeTab])}
      </div>
    );
  };


  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">COCONUT Pre-processing</h2>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          {/* SMILES Input Component */}
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            required
          // Pass theme props if SMILESInput needs them
          />

          {/* Options Checkboxes */}
          <div className="flex flex-col sm:flex-row sm:space-x-6 space-y-3 sm:space-y-0 pt-2">
            <div className="flex items-center">
              <input
                id="generate-3d"
                type="checkbox"
                checked={generate3D}
                onChange={(e) => setGenerate3D(e.target.checked)}
                // Checkbox styling for light/dark mode
                className="h-4 w-4 rounded bg-gray-50 dark:bg-gray-700 border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 focus:ring-blue-500 dark:focus:ring-offset-gray-800 shadow-sm"
              />
              <label htmlFor="generate-3d" className="ml-3 block text-sm text-gray-700 dark:text-gray-300">
                Generate 3D Coordinates
              </label>
            </div>
            <div className="flex items-center">
              <input
                id="generate-descriptors"
                type="checkbox"
                checked={generateDescriptors}
                onChange={(e) => setGenerateDescriptors(e.target.checked)}
                // Checkbox styling for light/dark mode
                className="h-4 w-4 rounded bg-gray-50 dark:bg-gray-700 border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 focus:ring-blue-500 dark:focus:ring-offset-gray-800 shadow-sm"
              />
              <label htmlFor="generate-descriptors" className="ml-3 block text-sm text-gray-700 dark:text-gray-300">
                Calculate Descriptors
              </label>
            </div>
          </div>

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              disabled={!smiles.trim() || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${!smiles.trim() || loading
                  ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                  : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm'
                }`}
            >
              <HiOutlineBeaker className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? 'Processing...' : 'Process for COCONUT'}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Pre-processing for COCONUT..." />}

      {/* Error Display */}
      {error && (
        <div className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow" role="alert">
          <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
          <span>{error}</span>
        </div>
      )}

      {/* Results Display Section */}
      {/* Show only if data exists and not loading */}
      {preprocessData && !loading && (
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
          {/* Results Header */}
          <div className="flex items-center mb-4 border-b border-gray-200 dark:border-gray-700 pb-3">
            <HiOutlineDocumentReport className="h-6 w-6 text-blue-600 dark:text-blue-400 mr-3 flex-shrink-0" aria-hidden="true" />
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white">COCONUT Pre-processing Results</h3>
          </div>
          {/* Render Tabs and Content */}
          {renderResults(preprocessData)}
        </div>
      )}

      {/* Initial State / About Box */}
      {/* Show only if no data, not loading, and no error */}
      {!preprocessData && !loading && !error && (
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 flex items-start space-x-4 shadow" role="complementary">
          <HiOutlineInformationCircle className="h-6 w-6 text-blue-600 dark:text-blue-400 flex-shrink-0 mt-0.5" aria-hidden="true" />
          <div>
            <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-2">About COCONUT Pre-processing</h3>
            <p className="text-gray-700 dark:text-gray-300">
              This tool prepares molecule data for submission to the COCONUT database (COlleCtion
              of Open Natural prodUcTs), a collection of natural products and natural product-like structures.
            </p>
            <p className="mt-3 text-gray-700 dark:text-gray-300">
              The pre-processing pipeline provides:
            </p>
            {/* List styling */}
            <ul className="list-disc list-inside mt-1 space-y-1 text-gray-600 dark:text-gray-400 marker:text-blue-500 dark:marker:text-blue-400">
              <li>Original molecule information</li>
              <li>Standardized molecule using ChEMBL curation pipeline</li>
              <li>Parent structure (salt and solvent stripped)</li>
              <li>Various chemical representations (SMILES, InChI, etc.)</li>
              <li>Optional 3D coordinates and molecular descriptors</li>
            </ul>
          </div>
        </div>
      )}
    </div>
  );
};

export default CoconutPreProcessingView;
