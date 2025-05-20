// Description: This component provides a user interface for filtering chemical compounds based on various medicinal chemistry rules and properties. It allows users to input SMILES strings, apply filters, and view results in a table format. The component is designed to be responsive and theme-aware, adapting to light and dark modes.
import React, { useState } from 'react';
// Import HiOutlineExclamationCircle along with other icons
import {
  HiOutlineFilter,
  HiOutlineClipboard,
  HiOutlineDocumentDownload,
  HiOutlineInformationCircle,
  HiCheck,
  HiX,
  HiOutlineExclamationCircle, // Added missing icon import
  HiOutlineSwitchHorizontal // Added for the filter operator toggle
} from 'react-icons/hi';
import LoadingScreen from '../common/LoadingScreen'; // Assuming this component is theme-aware
import { useAppContext } from '../../context/AppContext'; // Assuming this provides theme-aware context if needed
import api from '../../services/api'; // Assuming api service is configured

const AllFiltersView = () => {
  const [smilesInput, setSmilesInput] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [results, setResults] = useState([]);
  const [copied, setCopied] = useState(false);
  const [filterOperator, setFilterOperator] = useState('AND'); // Default to AND logic
  const { addRecentMolecule } = useAppContext(); // Assuming context provides this function

  // Filter options state
  const [filterOptions, setFilterOptions] = useState({
    pains: true,
    lipinski: true,
    veber: true,
    reos: true,
    ghose: true,
    ruleofthree: true,
    // Note: The API likely expects boolean for checkboxes,
    // and specific string formats or numbers for ranges.
    // These initial values might need adjustment based on API requirements.
    qedscore: '0-1', // Example range string
    sascore: '1-10', // Example range string (Corrected typical range)
    nplikeness: '-5-5' // Example range string
  });

  // Function to apply filters via API
  const applyFilters = async (e) => {
    e.preventDefault();
    const trimmedInput = smilesInput.trim();
    if (!trimmedInput) {
      setError('Please enter at least one SMILES string.');
      return;
    }

    setLoading(true);
    setError(null);
    setResults([]); // Clear previous results

    try {
      // Prepare parameters for the API call
      // Adjust the params object construction based on your actual API requirements.
      const apiParams = { 
        ...filterOptions,
        filterOperator: filterOperator // Add filter operator (AND/OR) to the API params
      };

      // Assuming api.post sends text data correctly and handles params
      const response = await api.post('/chem/all_filters', trimmedInput, {
        params: apiParams, // Send filter options as query parameters
        headers: {
          'Content-Type': 'text/plain', // Sending SMILES list as plain text
          'Accept': 'application/json' // Expecting JSON response, adjust if needed
        }
      });

      // Assuming the response data is the array of result strings
      // Adjust based on actual API response structure
      if (Array.isArray(response.data)) {
        setResults(response.data);
      } else {
        // Handle unexpected response format
        console.warn("Received non-array response:", response.data);
        setResults([]); // Set to empty or handle appropriately
        setError("Received unexpected data format from the server.");
      }


      // Add molecules to recent list (only if results were successful)
      if (Array.isArray(response.data)) {
        const smilesList = trimmedInput.split(/[\n\s,;]+/).filter(s => s.trim()); // Split by various delimiters
        smilesList.forEach(smiles => {
          const trimmedSmiles = smiles.trim();
          if (trimmedSmiles) {
            // Extract just the SMILES part if input format includes names/IDs
            const smilesOnly = trimmedSmiles.split(/[\s:]+/)[0];
            addRecentMolecule({
              smiles: smilesOnly,
              name: `Filtered Molecule (${new Date().toLocaleTimeString()})`, // Add a generic name/timestamp
              timestamp: new Date().toISOString()
            });
          }
        });
      }

    } catch (err) {
      console.error('Filter API error:', err);
      // Extract more specific error message if available (e.g., from response)
      const errorMsg = err.response?.data?.detail || err.message || 'An unknown error occurred.';
      setError(`Error applying filters: ${errorMsg}`);
      setResults([]); // Clear results on error
    } finally {
      setLoading(false);
    }
  };

  // Handle changes to filter options (generic)
  const handleFilterChange = (name, value) => {
    setFilterOptions(prev => ({
      ...prev,
      [name]: value
    }));
  };

  // Handle checkbox change specifically
  const handleCheckboxChange = (e) => {
    const { name, checked } = e.target;
    handleFilterChange(name, checked);
  };

  // Handle select/range input change specifically
  const handleSelectChange = (e) => {
    const { name, value } = e.target;
    handleFilterChange(name, value);
  };

  // Copy results to clipboard
  const copyToClipboard = () => {
    if (results.length === 0 || !navigator.clipboard) return;

    const text = results.join('\n'); // Join results array into a string
    navigator.clipboard.writeText(text)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000); // Reset copied state
      })
      .catch(err => {
        console.error('Failed to copy results:', err);
        setError('Failed to copy results to clipboard.');
      });
  };

  // Download results as a text file
  const downloadResults = () => {
    if (results.length === 0) return;

    const text = results.join('\n');
    const blob = new Blob([text], { type: 'text/plain;charset=utf-8' }); // Specify charset
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'filtered_molecules_results.txt'; // More descriptive filename
    document.body.appendChild(a); // Append to body
    a.click(); // Trigger download
    document.body.removeChild(a); // Clean up
    URL.revokeObjectURL(url); // Release object URL
  };

  // Parse and format the filter result string for display in the table
  // Adjust this function based on the exact format returned by your API
  const parseFilterResult = (resultString) => {
    if (typeof resultString !== 'string' || !resultString.includes(':')) {
      // Handle unexpected format or non-string results
      console.warn("Parsing unexpected result format:", resultString);
      return { smiles: resultString || 'Invalid Data', results: Array(9).fill('N/A') }; // Return placeholder
    }

    // Example format: "CCO : T, F, T, T, F, T, 0.5, 3.2, 1.1"
    const parts = resultString.split(' : ');
    const smiles = parts[0].trim();
    const filterDataString = parts[1] || ''; // Handle cases with no data after colon

    // Split the filter data part by comma and trim whitespace
    const filterValues = filterDataString.split(',').map(item => item.trim());

    // Map 'T'/'F' to boolean for the first 6 filters, keep others as strings/numbers
    // Ensure we have enough values, pad with 'N/A' if necessary
    const expectedLength = 9; // Adjust if the number of filters changes
    const parsedResults = [];
    for (let i = 0; i < expectedLength; i++) {
      const value = filterValues[i];
      if (value === undefined || value === '') {
        parsedResults.push('N/A'); // Handle missing values
      } else if (i < 6) { // Assuming first 6 are boolean T/F filters
        parsedResults.push(value.toUpperCase() === 'T');
      } else {
        parsedResults.push(value); // Keep scores/other values as strings (or parse to number if needed)
      }
    }

    return { smiles, results: parsedResults };
  };


  // Define column headers - makes the table more maintainable
  const filterColumns = [
    { key: 'pains', label: 'PAINS' },
    { key: 'lipinski', label: 'Lipinski' },
    { key: 'veber', label: 'Veber' },
    { key: 'reos', label: 'REOS' },
    { key: 'ghose', label: 'Ghose' },
    { key: 'ruleofthree', label: 'Rule of 3' },
    { key: 'qedscore', label: 'QED' }, // Shortened label
    { key: 'sascore', label: 'SA Score' },
    { key: 'nplikeness', label: 'NP-like' }, // Shortened label
  ];


  return (
    // Main container with spacing
    <div className="space-y-6 p-4 md:p-6">
      {/* Input and Filters Card */}
      <div className="bg-white dark:bg-gray-800 p-4 md:p-6 rounded-lg shadow-md dark:shadow-lg">
        {/* Card Header */}
        <div className="flex justify-between items-start mb-4">
          <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400">Chemical Filters</h2>
          {/* Info button - consider adding a tooltip or modal onClick */}
          <button
            className="p-1 text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white rounded-full focus:outline-none focus:ring-2 focus:ring-blue-500"
            title="Information about filters (Not implemented)" // Add tooltip text
            onClick={() => alert('Filter information modal/tooltip not implemented yet.')} // Placeholder action
          >
            <HiOutlineInformationCircle className="h-5 w-5" aria-hidden="true" />
          </button>
        </div>

        {/* Form for input and filter options */}
        <form onSubmit={applyFilters} className="space-y-6">
          {/* SMILES Input Area */}
          <div>
            <label htmlFor="smilesInputArea" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
              Molecule List (SMILES)
            </label>
            <textarea
              id="smilesInputArea"
              value={smilesInput}
              onChange={(e) => setSmilesInput(e.target.value)}
              placeholder="Enter SMILES strings, one per line or separated by space/comma/semicolon..."
              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 h-32 shadow-sm"
              required
              aria-required="true"
            />
            <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
              Enter one or more SMILES strings. Separate multiple entries with new lines, spaces, commas, or semicolons.
            </p>
          </div>

          {/* Filter Options Grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 border-t border-gray-200 dark:border-gray-700 pt-6">
            {/* Checkbox Filters */}
            <div className="space-y-4">
              <h3 className="text-md font-medium text-gray-800 dark:text-gray-200">Rule-Based Filters</h3>
              <div className="space-y-3">
                {/* Dynamically create checkboxes */}
                {filterColumns.slice(0, 6).map(filter => (
                  <div key={filter.key} className="flex items-center">
                    <input
                      id={filter.key}
                      name={filter.key}
                      type="checkbox"
                      checked={filterOptions[filter.key]}
                      onChange={handleCheckboxChange}
                      className="h-4 w-4 rounded bg-gray-50 dark:bg-gray-700 border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 focus:ring-blue-500 dark:focus:ring-offset-gray-800 shadow-sm"
                    />
                    <label htmlFor={filter.key} className="ml-3 block text-sm text-gray-700 dark:text-gray-300">
                      {filter.label}
                    </label>
                  </div>
                ))}
              </div>
            </div>

            {/* Score Range Filters */}
            <div className="space-y-4">
              <h3 className="text-md font-medium text-gray-800 dark:text-gray-200">Score Filters</h3>
              {/* QED Score */}
              <div>
                <label htmlFor="qedscore" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  QED Druglikeness (0-1)
                </label>
                <select
                  id="qedscore"
                  name="qedscore"
                  value={filterOptions.qedscore}
                  onChange={handleSelectChange} // Use specific handler
                  className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                >
                  <option value="0-1">0-1 (All)</option>
                  <option value="0-0.3">0-0.3 (Low)</option>
                  <option value="0.3-0.6">0.3-0.6 (Medium)</option>
                  <option value="0.6-1">0.6-1 (High)</option>
                </select>
              </div>
              {/* SA Score */}
              <div>
                <label htmlFor="sascore" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  Synthetic Accessibility (1-10) {/* Corrected range usually 1-10 */}
                </label>
                <select
                  id="sascore"
                  name="sascore"
                  value={filterOptions.sascore}
                  onChange={handleSelectChange}
                  className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                >
                  <option value="1-10">1-10 (All)</option> {/* Adjusted range */}
                  <option value="1-4">1-4 (Easy synthesis)</option>
                  <option value="4-7">4-7 (Moderate synthesis)</option>
                  <option value="7-10">7-10 (Difficult synthesis)</option>
                </select>
              </div>
              {/* NP-likeness Score */}
              <div>
                <label htmlFor="nplikeness" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                  NP-likeness Score (-5 to 5)
                </label>
                <select
                  id="nplikeness"
                  name="nplikeness"
                  value={filterOptions.nplikeness}
                  onChange={handleSelectChange}
                  className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                >
                  <option value="-5-5">-5 to 5 (All)</option>
                  <option value="-5-0">-5 to 0 (Synthetic-like)</option>
                  <option value="0-2">0 to 2 (Moderate NP-like)</option>
                  <option value="2-5">2 to 5 (High NP-like)</option>
                </select>
              </div>
            </div>
          </div>

          {/* Filter Operator Toggle */}
          <div className="pt-4 pb-2 border-t border-gray-200 dark:border-gray-700">
            <div className="flex flex-col sm:flex-row sm:items-center space-y-2 sm:space-y-0 sm:space-x-4">
              <label className="flex items-center text-sm font-medium text-gray-700 dark:text-gray-300">
                <HiOutlineSwitchHorizontal className="mr-2 h-5 w-5 text-blue-600 dark:text-blue-400" />
                Filter Match Logic:
              </label>
              <div className="flex p-1 space-x-1 bg-gray-100 dark:bg-gray-700 rounded-lg">
                <button
                  type="button"
                  onClick={() => setFilterOperator('AND')}
                  className={`px-3 py-1.5 text-sm font-medium rounded-md transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-blue-500 ${
                    filterOperator === 'AND'
                      ? 'bg-white dark:bg-gray-600 text-blue-600 dark:text-blue-400 shadow-sm'
                      : 'text-gray-600 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-600'
                  }`}
                >
                  Match All Filters (AND)
                </button>
                <button
                  type="button"
                  onClick={() => setFilterOperator('OR')}
                  className={`px-3 py-1.5 text-sm font-medium rounded-md transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-blue-500 ${
                    filterOperator === 'OR'
                      ? 'bg-white dark:bg-gray-600 text-blue-600 dark:text-blue-400 shadow-sm'
                      : 'text-gray-600 dark:text-gray-300 hover:bg-gray-200 dark:hover:bg-gray-600'
                  }`}
                >
                  Match Any Filter (OR)
                </button>
              </div>
            </div>
          </div>

          {/* Submit Button */}
          <div className="pt-4">
            <button
              type="submit"
              disabled={loading || !smilesInput.trim()}
              className={`w-full md:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${!smilesInput.trim() || loading
                  ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                  : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm'
                }`}
            >
              <HiOutlineFilter className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? 'Applying Filters...' : 'Apply Filters'}
            </button>
          </div>
        </form>
      </div>

      {/* Loading Indicator */}
      {loading && <LoadingScreen text="Applying filters to molecules..." />}

      {/* Error Message */}
      {error && (
        // Use the imported HiOutlineExclamationCircle here
        <div className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow" role="alert">
          <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
          <div>
            <h4 className="font-medium">Error</h4>
            <p className="text-sm">{error}</p>
          </div>
        </div>
      )}

      {/* Results Table Section */}
      {/* Show only if not loading and results array has items */}
      {!loading && results.length > 0 && (
        <div className="bg-white dark:bg-gray-800 p-4 md:p-6 rounded-lg shadow-md dark:shadow-lg">
          {/* Results Header with Actions */}
          <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center mb-4 gap-3">
            <h3 className="text-lg font-medium text-gray-900 dark:text-white">
              Filter Results <span className="text-sm font-normal text-gray-500 dark:text-gray-400">({results.length} molecules passing selected criteria)</span>
            </h3>
            {/* Action Buttons */}
            <div className="flex space-x-2 flex-shrink-0">
              <button
                onClick={copyToClipboard}
                className={`px-3 py-1.5 text-sm rounded-md flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${copied ? 'bg-green-100 dark:bg-green-700 text-green-700 dark:text-green-200' : 'bg-gray-100 hover:bg-gray-200 dark:bg-gray-700 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600'}`}
                title="Copy results list to clipboard"
              >
                <HiOutlineClipboard className={`mr-1.5 h-4 w-4 ${copied ? '' : ''}`} aria-hidden="true" />
                {copied ? "Copied!" : "Copy"}
              </button>
              <button
                onClick={downloadResults}
                className="px-3 py-1.5 text-sm bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 text-white rounded-md flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500"
                title="Download results as a text file"
              >
                <HiOutlineDocumentDownload className="mr-1.5 h-4 w-4" aria-hidden="true" />
                Download
              </button>
            </div>
          </div>

          {/* Results Table */}
          <div className="overflow-x-auto border border-gray-200 dark:border-gray-700 rounded-md">
            <table className="w-full min-w-[800px] border-collapse"> {/* min-w ensures table scrolls */}
              <thead className="bg-gray-50 dark:bg-gray-700">
                <tr>
                  {/* Fixed Columns */}
                  <th scope="col" className="sticky left-0 z-10 bg-gray-50 dark:bg-gray-700 px-3 py-2 text-left text-xs font-medium uppercase tracking-wider text-gray-500 dark:text-gray-400 border-r border-gray-200 dark:border-gray-600">#</th>
                  <th scope="col" className="sticky left-[calc(theme(space.12))] z-10 bg-gray-50 dark:bg-gray-700 px-4 py-2 text-left text-xs font-medium uppercase tracking-wider text-gray-500 dark:text-gray-400 border-r border-gray-200 dark:border-gray-600">SMILES</th>
                  {/* Dynamic Filter Columns */}
                  {filterColumns.map((col) => (
                    <th key={col.key} scope="col" className="px-3 py-2 text-center text-xs font-medium uppercase tracking-wider text-gray-500 dark:text-gray-400 whitespace-nowrap">{col.label}</th>
                  ))}
                </tr>
              </thead>
              <tbody className="bg-white dark:bg-gray-800 divide-y divide-gray-200 dark:divide-gray-700">
                {results.map((resultString, index) => {
                  const { smiles, results: filterValues } = parseFilterResult(resultString);
                  return (
                    <tr key={index} className="hover:bg-gray-50 dark:hover:bg-gray-700/50 transition-colors duration-150">
                      {/* Row Number (Sticky) */}
                      <td className="sticky left-0 z-10 bg-white dark:bg-gray-800 px-3 py-2 text-sm text-gray-500 dark:text-gray-400 border-r border-gray-200 dark:border-gray-600">{index + 1}</td>
                      {/* SMILES String (Sticky) */}
                      <td className="sticky left-[calc(theme(space.12))] z-10 bg-white dark:bg-gray-800 px-4 py-2 text-sm text-blue-600 dark:text-blue-300 font-mono whitespace-nowrap border-r border-gray-200 dark:border-gray-600">
                        {smiles}
                      </td>
                      {/* Filter Results */}
                      {filterValues.map((value, filterIndex) => (
                        <td key={filterIndex} className="px-3 py-2 text-center text-sm">
                          {typeof value === 'boolean' ? (
                            value ? (
                              <HiCheck className="h-5 w-5 text-green-500 dark:text-green-400 inline-block" title="Pass" aria-label="Pass" />
                            ) : (
                              <HiX className="h-5 w-5 text-red-500 dark:text-red-400 inline-block" title="Fail" aria-label="Fail" />
                            )
                          ) : (
                            // Display scores/other values
                            <span className="text-gray-700 dark:text-gray-300">{value === 'N/A' ? '-' : value}</span>
                          )}
                        </td>
                      ))}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Initial State / Info Box */}
      {/* Show only if not loading, no error, and no results */}
      {/* Restored detailed about section with light/dark mode styling */}
      {!loading && !error && results.length === 0 && (
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 shadow" role="complementary">
          <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-3">About Chemical Filters</h3>
          <p className="text-gray-700 dark:text-gray-300 mb-4">
            This tool allows you to filter a list of molecules based on various medicinal chemistry and
            drug-like property filters.
          </p>
          <div className="space-y-2">
            <p className="text-sm text-gray-600 dark:text-gray-400"><strong>PAINS</strong>: Pan-Assay Interference compounds - structures known to interfere with biochemical assays</p>
            <p className="text-sm text-gray-600 dark:text-gray-400"><strong>Lipinski Rule of 5</strong>: Evaluates drug-likeness based on molecular weight, logP, H-bond donors/acceptors</p>
            <p className="text-sm text-gray-600 dark:text-gray-400"><strong>Veber</strong>: Filters for oral bioavailability based on rotatable bonds and polar surface area</p>
            <p className="text-sm text-gray-600 dark:text-gray-400"><strong>REOS</strong>: Rapid Elimination Of Swill - property filters for lead-like compounds</p>
            <p className="text-sm text-gray-600 dark:text-gray-400"><strong>Ghose</strong>: Filters based on logP, molecular weight, and number of atoms</p>
            <p className="text-sm text-gray-600 dark:text-gray-400"><strong>Rule of 3</strong>: Criteria for fragment-based drug discovery</p>
            <p className="mt-4 text-sm text-gray-600 dark:text-gray-400 border-t border-blue-100 dark:border-blue-800 pt-4">
              <strong>Filter Match Logic</strong>: Choose how filters are combined - "Match All Filters" requires molecules to pass all selected filters (AND logic), while "Match Any Filter" includes molecules that pass at least one filter (OR logic).
            </p>
          </div>
        </div>
      )}
    </div>
  );
};

export default AllFiltersView;
