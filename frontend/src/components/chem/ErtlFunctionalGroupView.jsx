// Description: This component allows users to input a SMILES string and detects functional groups using the Ertl algorithm.
import React, { useState } from "react";
import {
  HiOutlineSearch,
  HiOutlineInformationCircle,
  HiOutlineExclamationCircle,
} from "react-icons/hi";
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import LoadingScreen from "../common/LoadingScreen";
import HighlightedMoleculeCard from "../common/HighlightedMoleculeCard";
// Assuming this service is configured correctly
import { generateFunctionalGroups } from "../../services/chemService";

const ErtlFunctionalGroupView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [functionalGroups, setFunctionalGroups] = useState([]); // Initialize as empty array
  const [selectedGroupIndex, setSelectedGroupIndex] = useState(null); // Track which group is selected for highlighting

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setFunctionalGroups([]); // Clear previous results
      return;
    }

    setLoading(true);
    setError(null);
    setFunctionalGroups([]); // Clear previous results before fetching

    try {
      const result = await generateFunctionalGroups(trimmedSmiles);
      // Ensure result is always an array
      const groupsArray = Array.isArray(result) ? result : [];
      setFunctionalGroups(groupsArray);

      // Check if the result indicates no groups found
      if (groupsArray.length === 0 || (groupsArray.length === 1 && groupsArray[0]?.None)) {
        // Don't set an error, just let the results section handle the display
        // setError('No functional groups found in this molecule.');
      }
    } catch (err) {
      console.error("Functional group detection error:", err); // Log the error
      setError(`Error detecting functional groups: ${err.message || "An unknown error occurred."}`);
      setFunctionalGroups([]); // Ensure groups are empty on error
    } finally {
      setLoading(false);
    }
  };

  // Function to format an IFG object (or string) to a readable string
  const formatFunctionalGroup = (fg) => {
    // Handle simple cases first
    if (typeof fg === "string") return fg;
    if (!fg || typeof fg !== "object") return "Invalid group data"; // Handle null/undefined/non-object
    if (fg.None) return "No functional groups found"; // Specific check for 'None' object

    // Handle new structured format from backend
    if (fg.description) {
      return fg.description; // Use the full description for display
    }

    // Try to extract standard properties (adjust based on actual API response)
    if (fg.type && fg.atoms) {
      return `${fg.type} (atoms: ${fg.atoms})`; // Example format
    }
    if (fg.type && fg.atomIds) {
      return `${fg.type} (atom IDs: ${fg.atomIds.join(", ")})`;
    }
    if (fg.type) {
      return fg.type; // Fallback if only type is present
    }

    // Fallback: Construct a string from available properties, excluding 'None'
    const props = Object.entries(fg)
      .filter(([key]) => key !== "None")
      .map(([key, value]) => `${key}: ${value}`)
      .join(", ");

    return props || JSON.stringify(fg); // Return formatted props or raw JSON stringify
  };

  // Determine if the results indicate "no groups found"
  const noGroupsFound =
    functionalGroups.length === 0 || (functionalGroups.length === 1 && functionalGroups[0]?.None);

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Functional Group Detection (Ertl)
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
              <HiOutlineSearch className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? "Detecting..." : "Detect Functional Groups"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Detecting functional groups..." />}

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
      {/* Show only if results exist (even if empty/None) and not loading and no error */}
      {functionalGroups.length > 0 && !loading && !error && (
        <div className="space-y-4">
          {/* Main Results Card */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
            {/* Results Summary Header */}
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-4 border-b border-gray-200 dark:border-gray-700 pb-2">
              {noGroupsFound
                ? "No Functional Groups Found"
                : `Found ${functionalGroups.length} Functional Group${functionalGroups.length !== 1 ? "s" : ""}`}
            </h3>

            {/* Grid for Molecule Card and List */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              {/* Molecule Card Column */}
              <div>
                <div className="mb-4">
                  {/* Use highlighted molecule card for better visualization */}
                  <HighlightedMoleculeCard
                    smiles={smiles}
                    title="Input Molecule"
                    size="md"
                    functionalGroups={functionalGroups}
                    highlightedGroupIndex={selectedGroupIndex}
                    showActions={true}
                  />
                </div>
                {selectedGroupIndex !== null &&
                  functionalGroups[selectedGroupIndex] &&
                  !functionalGroups[selectedGroupIndex].None && (
                    <div className="text-center">
                      <button
                        onClick={() => setSelectedGroupIndex(null)}
                        className="px-3 py-1 text-xs bg-gray-200 dark:bg-gray-700 text-gray-700 dark:text-gray-300 rounded-full hover:bg-gray-300 dark:hover:bg-gray-600 transition-colors"
                      >
                        Clear Highlighting
                      </button>
                    </div>
                  )}
              </div>

              {/* Functional Groups List Column */}
              <div>
                <h4 className="text-md font-medium text-gray-800 dark:text-blue-300 mb-2">
                  Identified Groups
                </h4>
                {/* Container for the list */}
                <div className="bg-gray-50 dark:bg-gray-900 rounded-lg p-4 border border-gray-200 dark:border-gray-700 shadow-sm min-h-[100px]">
                  {noGroupsFound ? (
                    <p className="text-gray-500 dark:text-gray-400 italic">
                      No functional groups detected.
                    </p>
                  ) : (
                    // List of functional groups
                    <ul className="space-y-2">
                      {functionalGroups.map((fg, index) => (
                        <li
                          key={index}
                          // List item styling with clickable highlighting
                          className={`p-2 border rounded-md text-sm cursor-pointer transition-all duration-200 ${
                            selectedGroupIndex === index
                              ? "bg-blue-100 dark:bg-blue-900 border-blue-300 dark:border-blue-700 text-blue-800 dark:text-blue-200"
                              : "bg-white dark:bg-gray-800 border-gray-200 dark:border-gray-700 text-gray-700 dark:text-gray-300 hover:bg-gray-50 dark:hover:bg-gray-750"
                          }`}
                          onClick={() => {
                            if (fg.None) return; // Don't allow selecting "None" groups
                            setSelectedGroupIndex(selectedGroupIndex === index ? null : index);
                          }}
                          title={
                            fg.None
                              ? undefined
                              : "Click to highlight this functional group in the structure"
                          }
                        >
                          <div className="flex items-center justify-between">
                            <span>{formatFunctionalGroup(fg)}</span>
                            {!fg.None && (
                              <div className="flex items-center space-x-2">
                                {fg.atomIds && fg.atomIds.length > 0 && (
                                  <span className="text-xs opacity-70">
                                    {fg.atomIds.length} atom
                                    {fg.atomIds.length !== 1 ? "s" : ""}
                                  </span>
                                )}
                                {selectedGroupIndex === index && (
                                  <span className="text-xs font-medium">Highlighted</span>
                                )}
                              </div>
                            )}
                          </div>
                          {/* Show additional details for structured groups */}
                          {fg.atoms && fg.type && !fg.None && (
                            <div className="mt-1 text-xs opacity-70">
                              Type: {fg.type} | Atoms: {fg.atoms}
                            </div>
                          )}
                        </li>
                      ))}
                    </ul>
                  )}
                </div>
              </div>
            </div>
          </div>

          {/* Informational Note about Ertl */}
          <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 text-sm text-gray-700 dark:text-gray-300 shadow">
            <p>
              Functional groups are identified using the algorithm proposed by Peter Ertl, analyzing
              molecular structure for common chemical groups.
            </p>
            {!noGroupsFound && (
              <p className="mt-2 text-xs text-blue-700 dark:text-blue-300">
                💡 <strong>Tip:</strong> Click on any functional group in the list to highlight it
                in the molecular structure above.
              </p>
            )}
            {/* Optionally add reference back here if desired */}
            {/* <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">Ref: Ertl, P. J Cheminform 9, 36 (2017).</p> */}
          </div>
        </div>
      )}

      {/* Initial State / About Box */}
      {/* Show only if no results, not loading, and no error */}
      {!functionalGroups.length && !loading && !error && (
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
              About Functional Group Detection (Ertl)
            </h3>
            <p className="text-gray-700 dark:text-gray-300">
              This tool identifies functional groups in organic molecules using the algorithm
              proposed by Peter Ertl. Functional groups are specific arrangements of atoms that give
              a molecule its characteristic chemical properties.
            </p>
            <p className="mt-3 text-gray-700 dark:text-gray-300">Examples include:</p>
            {/* List styling */}
            <ul className="list-disc list-inside mt-1 space-y-1 text-gray-600 dark:text-gray-400 marker:text-blue-500 dark:marker:text-blue-400">
              <li>Hydroxyl groups (-OH)</li>
              <li>Carbonyl groups (C=O)</li>
              <li>Carboxyl groups (-COOH)</li>
              <li>Amino groups (-NH₂)</li>
              <li>And many others</li>
            </ul>
            <p className="mt-3 text-xs text-gray-500 dark:text-gray-400">
              Reference: Ertl, Peter. "An algorithm to identify functional groups in organic
              molecules." Journal of Cheminformatics 9.1 (2017): 36.
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default ErtlFunctionalGroupView;
