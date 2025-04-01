// Description: This component provides a user interface for detecting and removing sugar moieties from molecules using SMILES strings. It includes input handling, API calls for sugar detection and removal, and displays results in a user-friendly format.
import React, { useState } from 'react';
// Ensure all used icons are imported
import {
  HiOutlineSearch,
  HiOutlineTrash,
  HiOutlineInformationCircle,
  HiOutlineExclamationCircle // Added for error display
} from 'react-icons/hi';
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from '../common/SMILESInput';
import MoleculeCard from '../common/MoleculeCard';
import LoadingScreen from '../common/LoadingScreen';
import SMILESDisplay from '../common/SMILESDisplay';
// Assuming this service is configured correctly
import { detectSugars, removeSugars, removeLinearSugars, removeCircularSugars } from '../../services/toolsService'; // Assuming this service exists

const SugarRemovalView = () => {
  const [smiles, setSmiles] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [sugarInfo, setSugarInfo] = useState(null); // Stores result from detectSugars
  const [removedSugarSmiles, setRemovedSugarSmiles] = useState(null); // Stores result from removeSugars
  const [removalType, setRemovalType] = useState('all'); // 'all', 'linear', 'circular'

  // Handle sugar detection
  const handleDetectSugars = async (e) => {
    if (e) e.preventDefault(); // Prevent form submission if called from button
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError('Please enter a SMILES string');
      setSugarInfo(null); setRemovedSugarSmiles(null);
      return;
    }

    setIsLoading(true);
    setError(null);
    setSugarInfo(null);
    setRemovedSugarSmiles(null); // Clear previous removal result

    try {
      // Check if service is available
      if (!detectSugars || typeof detectSugars !== 'function') {
        throw new Error("Sugar detection service is not available.");
      }
      const info = await detectSugars(trimmedSmiles);
      setSugarInfo(info || "No specific sugar information returned."); // Provide default message if null/undefined
    } catch (err) {
      console.error("Error detecting sugars:", err);
      setError(`Error detecting sugars: ${err.message || 'Unknown error'}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Handle sugar removal based on selected type
  const handleRemoveSugars = async () => {
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError('Please enter a SMILES string');
      setSugarInfo(null); setRemovedSugarSmiles(null);
      return;
    }

    setIsLoading(true);
    setError(null);
    setSugarInfo(null); // Clear previous detection result
    setRemovedSugarSmiles(null);

    try {
      let result;
      let removeFunction;

      // Select the appropriate service function based on removalType
      switch (removalType) {
        case 'linear': removeFunction = removeLinearSugars; break;
        case 'circular': removeFunction = removeCircularSugars; break;
        case 'all': default: removeFunction = removeSugars; break;
      }

      // Check if service function exists
      if (!removeFunction || typeof removeFunction !== 'function') {
        throw new Error(`Sugar removal service for type '${removalType}' is not available.`);
      }

      result = await removeFunction(trimmedSmiles);

      // Check if the result is a valid SMILES string (basic check)
      // The API might return messages like "No sugars found" instead of SMILES
      if (typeof result === 'string' && result.trim() && !result.toLowerCase().startsWith('no')) {
        setRemovedSugarSmiles(result);
        // Set an info message if the result is the same as input (no sugars removed)
        if (result === trimmedSmiles) {
          setError(`No ${removalType !== 'all' ? removalType + ' ' : ''}sugars found to remove.`);
        }
      } else {
        // API returned a message or empty result
        setError(result || `No ${removalType !== 'all' ? removalType + ' ' : ''}sugars found to remove.`);
        setRemovedSugarSmiles(null);
      }
    } catch (err) {
      console.error("Error removing sugars:", err);
      setError(`Error removing sugars: ${err.message || 'Unknown error'}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Determine if the error state contains an informational message vs. a real error
  const isInfoMessage = typeof error === 'string' && (error.toLowerCase().includes('no sugar') || error.toLowerCase().includes('already canonical'));

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Sugar Moiety Detection & Removal
        </h2>

        {/* Using handleDetectSugars for form submission */}
        <form onSubmit={handleDetectSugars} className="space-y-4">
          {/* SMILES Input */}
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            label="Input SMILES"
            required
          />

          {/* Removal Type Radio Buttons */}
          <fieldset className="pt-2">
            <legend className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
              Sugar Removal Type (for 'Remove Sugars' action)
            </legend>
            {/* Radio button group styling */}
            <div className="flex flex-col sm:flex-row flex-wrap gap-x-6 gap-y-2">
              {['all', 'linear', 'circular'].map((type) => (
                <label key={type} className="flex items-center space-x-2 cursor-pointer">
                  <input
                    type="radio"
                    value={type}
                    checked={removalType === type}
                    onChange={() => setRemovalType(type)}
                    // Radio button styling
                    className="h-4 w-4 border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 focus:ring-blue-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
                  />
                  {/* Label text styling */}
                  <span className="text-sm text-gray-700 dark:text-gray-300 capitalize">
                    {type === 'all' ? 'All Sugars' : `${type} Sugars Only`}
                  </span>
                </label>
              ))}
            </div>
          </fieldset>

          {/* Action Buttons */}
          <div className="flex flex-wrap gap-3 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Detect Button */}
            <button
              type="submit" // Triggers handleDetectSugars via form onSubmit
              disabled={!smiles.trim() || isLoading}
              // Button styling
              className={`px-5 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${!smiles.trim() || isLoading
                ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm'
                }`}
            >
              <HiOutlineSearch className="mr-2 h-5 w-5" />
              {isLoading ? 'Detecting...' : 'Detect Sugars'}
            </button>

            {/* Remove Button */}
            <button
              type="button" // Important: type="button" to prevent form submission
              onClick={handleRemoveSugars}
              disabled={!smiles.trim() || isLoading}
              // Button styling (Red color)
              className={`px-5 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-red-500 ${!smiles.trim() || isLoading
                ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                : 'bg-red-600 hover:bg-red-700 dark:bg-red-500 dark:hover:bg-red-600 shadow-sm'
                }`}
            >
              <HiOutlineTrash className="mr-2 h-5 w-5" />
              {isLoading ? 'Removing...' : 'Remove Sugars'}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {isLoading && <LoadingScreen text="Processing molecule..." />}

      {/* Error/Info Display */}
      {error && !isLoading && (
        // Use blue for info messages, red for errors
        <div className={`p-4 rounded-md flex items-start shadow ${isInfoMessage
          ? 'bg-blue-50 dark:bg-blue-900 dark:bg-opacity-30 text-blue-700 dark:text-blue-200 border border-blue-300 dark:border-blue-700'
          : 'bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700'
          }`} role={isInfoMessage ? 'status' : 'alert'}>
          {isInfoMessage
            ? <HiOutlineInformationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-blue-500 dark:text-blue-400" aria-hidden="true" />
            : <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
          }
          <span>{error}</span>
        </div>
      )}

      {/* Sugar Detection Result Display */}
      {sugarInfo && !isLoading && !error && ( // Show only if no error
        // Info box styling
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 shadow">
          <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-2">Sugar Detection Result</h3>
          {/* Result text styling */}
          <p className="text-gray-700 dark:text-gray-200">{sugarInfo}</p>
        </div>
      )}

      {/* Sugar Removal Result Display */}
      {removedSugarSmiles && !isLoading && !isInfoMessage && (
        <div className="space-y-4">
          <h3 className="text-lg font-medium text-gray-800 dark:text-gray-200">
            Sugar Removal Result
          </h3>

          {/* SMILES Displays - Added section showing SMILES with copy/download functionality */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
            <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
              <SMILESDisplay
                smiles={smiles}
                label="Original SMILES"
              />
            </div>
            <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
              <SMILESDisplay
                smiles={removedSugarSmiles}
                label={`Sugar-Free SMILES (${removalType})`}
              />
            </div>
          </div>

          {/* Molecule visualization cards */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
              <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-2">Original Molecule</h4>
              <MoleculeCard
                smiles={smiles}
                title="Original"
                description="Molecule with potential sugar moieties"
                showActions={false} // Hide actions as we now have SMILESDisplay
              />
            </div>

            <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
              <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-2">Sugar-Free Molecule</h4>
              <MoleculeCard
                smiles={removedSugarSmiles}
                title="Sugar-Free"
                description={`After ${removalType === 'all' ? 'all' : removalType} sugar removal`}
                showActions={false} // Hide actions as we now have SMILESDisplay
              />
            </div>
          </div>
        </div>
      )}

      {/* Initial State / About Box */}
      {/* Show only if no results, not loading, and no error */}
      {!sugarInfo && !removedSugarSmiles && !isLoading && !error && (
        // About box styling
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 flex items-start space-x-4 shadow" role="complementary">
          <HiOutlineInformationCircle className="h-6 w-6 text-blue-600 dark:text-blue-400 flex-shrink-0 mt-0.5" aria-hidden="true" />
          <div>
            <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-2">About Sugar Removal</h3>
            {/* Text styling */}
            <p className="text-gray-700 dark:text-gray-300">
              This tool detects and removes sugar moieties (linear and circular) from molecules, often useful for natural product analysis.
            </p>
            <p className="mt-2 text-gray-700 dark:text-gray-300">
              Enter a SMILES string to first detect sugars. Then, choose a removal type and click "Remove Sugars" to generate the deglycosylated structure.
            </p>
            <p className="mt-3 text-xs text-gray-500 dark:text-gray-400">
              Reference: Schaub, J., Zielesny, A., Steinbeck, C. et al. J Cheminform 12, 67 (2020).
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default SugarRemovalView;