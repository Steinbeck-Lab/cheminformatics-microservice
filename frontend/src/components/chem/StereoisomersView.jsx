//Description: This component allows users to input a SMILES string and generates its stereoisomers using a service. It handles loading states, errors, and displays results in a user-friendly manner.
import React, { useState } from 'react';
// Ensure all used icons are imported
import {
  HiOutlineCube,
  HiOutlineExclamationCircle,
  HiOutlineInformationCircle // Added for consistency in info boxes
} from 'react-icons/hi';
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from '../common/SMILESInput';
import MoleculeCard from '../common/MoleculeCard';
import LoadingScreen from '../common/LoadingScreen';
// Assuming this service is configured correctly
import { generateStereoisomers } from '../../services/chemService'; // Assuming this service exists

const StereoisomersView = () => {
  const [smiles, setSmiles] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [stereoisomers, setStereoisomers] = useState([]); // Initialize as empty array

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError('Please enter a SMILES string');
      setStereoisomers([]); // Clear previous results
      return;
    }

    setLoading(true);
    setError(null);
    setStereoisomers([]); // Clear previous results before fetching

    try {
      // Call the service function
      const result = await generateStereoisomers(trimmedSmiles);
      // Ensure result is always an array
      const isomersArray = Array.isArray(result) ? result : [];
      setStereoisomers(isomersArray);

      // Set informational message if no isomers found (or only input returned)
      if (isomersArray.length === 0 || (isomersArray.length === 1 && isomersArray[0] === trimmedSmiles)) {
        setError('No distinct stereoisomers found for this molecule.');
        // Keep the input SMILES if it was the only result
        if (isomersArray.length === 1 && isomersArray[0] === trimmedSmiles) {
          setStereoisomers(isomersArray);
        } else {
          setStereoisomers([]); // Clear if API returned empty
        }
      }
    } catch (err) {
      console.error("Stereoisomer generation error:", err); // Log the error
      setError(`Error generating stereoisomers: ${err.message || 'An unknown error occurred.'}`);
      setStereoisomers([]); // Ensure isomers are empty on error
    } finally {
      setLoading(false);
    }
  };

  // Determine if the results indicate "no distinct isomers found" based on error state
  const noDistinctIsomersFound = error === 'No distinct stereoisomers found for this molecule.';
  // Determine if there was a real error
  const hasError = error && !noDistinctIsomersFound;

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">Stereoisomer Generator</h2>

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
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${!smiles.trim() || loading
                  ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                  : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm'
                }`}
            >
              <HiOutlineCube className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? 'Generating...' : 'Generate Stereoisomers'}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Generating stereoisomers..." />}

      {/* Error/Info Display */}
      {/* Show error/info only if not loading */}
      {error && !loading && (
        <div className={`p-4 rounded-md flex items-start shadow ${noDistinctIsomersFound // Check if it's the informational message
            ? 'bg-blue-50 dark:bg-blue-900 dark:bg-opacity-30 text-blue-700 dark:text-blue-200 border border-blue-300 dark:border-blue-700'
            : 'bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700'
          }`} role={noDistinctIsomersFound ? 'status' : 'alert'}>
          {noDistinctIsomersFound
            ? <HiOutlineInformationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-blue-500 dark:text-blue-400" aria-hidden="true" />
            : <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
          }
          <span>{error}</span>
        </div>
      )}


      {/* Results Display Section */}
      {/* Show only if results exist and not loading and no actual error occurred */}
      {stereoisomers.length > 0 && !loading && !hasError && (
        <div className="space-y-4">
          {/* Main Results Card */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
            {/* Results Header */}
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
              {/* Adjust count if input SMILES is included in the list when no distinct isomers are found */}
              Found {noDistinctIsomersFound ? 1 : stereoisomers.length} Stereoisomer{stereoisomers.length !== 1 || noDistinctIsomersFound ? 's' : ''}
            </h3>
            <p className="text-sm text-gray-600 dark:text-gray-400 mb-4">
              {noDistinctIsomersFound
                ? "The input molecule has no distinct stereoisomers."
                : "These are all possible stereoisomers (enantiomers and diastereomers) based on the input."
              }
            </p>

            {/* Grid for Molecule Cards */}
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4">
              {stereoisomers.map((isomer, index) => (
                // Ensure MoleculeCard is theme-aware
                <MoleculeCard
                  key={index}
                  smiles={isomer}
                  title={`Stereoisomer ${index + 1}`}
                  size="sm" // Use smaller cards for potentially many isomers
                // Pass theme props if needed
                />
              ))}
            </div>
          </div>

          {/* Informational Box about Stereoisomers */}
          <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 text-sm shadow">
            <h4 className="font-semibold text-blue-800 dark:text-blue-300 mb-1">About Stereoisomers</h4>
            <p className="text-gray-700 dark:text-gray-300">
              Stereoisomers have the same molecular formula and atom connectivity but differ in the 3D arrangement of atoms.
            </p>
            <p className="mt-2 text-gray-700 dark:text-gray-300">
              <strong>Enantiomers:</strong> Non-superimposable mirror images.
            </p>
            <p className="mt-1 text-gray-700 dark:text-gray-300">
              <strong>Diastereomers:</strong> Stereoisomers that are not mirror images.
            </p>
          </div>
        </div>
      )}

      {/* Initial State Message (Optional) */}
      {/* Show only if no results, not loading, and no error */}
      {!stereoisomers.length && !loading && !error && (
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 text-center shadow">
          <p className="text-gray-600 dark:text-gray-300">
            Enter a SMILES string to generate its possible stereoisomers.
          </p>
        </div>
      )}

    </div>
  );
};

export default StereoisomersView;
