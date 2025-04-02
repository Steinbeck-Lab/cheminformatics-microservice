// Description: A React component for inputting SMILES notation with example and recent molecule dropdowns.
import React, { useState, useRef, useEffect } from 'react';
// Import necessary icons
import {
  HiOutlineX,
  HiOutlineLightBulb,
  HiOutlineClipboard,
  HiOutlineClock // Replaced SVG for Recent
} from 'react-icons/hi';
import { useAppContext } from '../../context/AppContext'; // Assuming context provides recentMolecules

// Example molecules data (static, moved outside component)
const EXAMPLE_MOLECULES = [
  { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', description: 'Stimulant found in coffee' },
  { name: 'Aspirin', smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O', description: 'Analgesic & anti-inflammatory' },
  { name: 'Sucrose', smiles: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O', description: 'Common sugar' },
  { name: 'Cholesterol', smiles: 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C', description: 'Steroid in cell membranes' },
  { name: 'Paracetamol', smiles: 'CC(=O)NC1=CC=C(C=C1)O', description: 'Pain reliever & fever reducer' },
  { name: 'Ibuprofen', smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', description: 'Anti-inflammatory drug' },
];

// Dropdown List Item Component (for cleaner code)
const DropdownItem = ({ molecule, onClick }) => (
  <li
    className="px-3 py-2 hover:bg-gray-100 dark:hover:bg-gray-700 cursor-pointer rounded-md mx-1"
    onClick={() => onClick(molecule.smiles)}
    role="menuitem"
    tabIndex={-1} // Can be focused programmatically if needed
  >
    <p className="font-medium text-sm text-gray-800 dark:text-gray-100 truncate" title={molecule.name || `SMILES: ${molecule.smiles}`}>
      {molecule.name || 'Recent Molecule'}
    </p>
    <p className="text-xs text-gray-500 dark:text-gray-400 truncate" title={molecule.smiles}>
      {molecule.smiles}
    </p>
    {molecule.description && (
      <p className="text-xs text-gray-400 dark:text-gray-500 mt-1">{molecule.description}</p>
    )}
  </li>
);


const SMILESInput = ({
  value,
  onChange,
  placeholder = "Enter SMILES notation...",
  showExamples = true,
  showRecent = true,
  label = "SMILES",
  required = false
}) => {
  const [showExampleList, setShowExampleList] = useState(false);
  const [showRecentList, setShowRecentList] = useState(false);
  // Refs for detecting outside clicks
  const examplesRef = useRef(null);
  const recentRef = useRef(null);

  // Assuming context provides recentMolecules array [{ smiles, name?, timestamp? }]
  const { recentMolecules = [] } = useAppContext() || {}; // Default to empty array if context missing

  // Close dropdowns if clicked outside
  useEffect(() => {
    const handleClickOutside = (event) => {
      if (examplesRef.current && !examplesRef.current.contains(event.target)) {
        setShowExampleList(false);
      }
      if (recentRef.current && !recentRef.current.contains(event.target)) {
        setShowRecentList(false);
      }
    };
    document.addEventListener('mousedown', handleClickOutside);
    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, []);


  const handleClear = () => {
    onChange('');
  };

  const handleExampleClick = (smiles) => {
    onChange(smiles);
    setShowExampleList(false); // Close dropdown after selection
  };

  const handleRecentClick = (smiles) => {
    onChange(smiles);
    setShowRecentList(false); // Close dropdown after selection
  };

  const handlePaste = async () => {
    if (!navigator.clipboard?.readText) {
      console.error("Clipboard API not available or permission denied.");
      // Optionally show user feedback
      return;
    }
    try {
      const text = await navigator.clipboard.readText();
      onChange(text); // Update input value with pasted text
    } catch (err) {
      console.error('Failed to read clipboard:', err);
      // Optionally show user feedback
    }
  };

  return (
    <div className="space-y-2">
      {/* Label */}
      {label && (
        <label htmlFor="smiles-input" className="block text-sm font-medium text-gray-700 dark:text-gray-300">
          {label} {required && <span className="text-red-500 dark:text-red-400">*</span>}
        </label>
      )}

      {/* Input field with Clear button */}
      <div className="relative">
        <input
          id="smiles-input"
          type="text"
          value={value}
          onChange={(e) => onChange(e.target.value)}
          placeholder={placeholder}
          // Input styling with theme awareness
          className="w-full px-4 py-2 bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-700 rounded-md text-gray-900 dark:text-white shadow-sm focus:outline-none focus:border-indigo-500 dark:focus:border-blue-500 focus:ring-1 focus:ring-indigo-500 dark:focus:ring-blue-500 pr-10" // Added pr-10 for clear button space
          required={required}
          aria-required={required}
        />

        {/* Clear button */}
        {value && (
          <button
            type="button" // Prevent form submission
            onClick={handleClear}
            // Button styling
            className="absolute inset-y-0 right-0 pr-3 flex items-center text-gray-400 dark:text-gray-500 hover:text-gray-600 dark:hover:text-gray-300 focus:outline-none focus:ring-1 focus:ring-indigo-500 rounded-r-md"
            aria-label="Clear input"
          >
            <HiOutlineX className="h-5 w-5" />
          </button>
        )}
      </div>

      {/* Action Buttons and Dropdowns */}
      <div className="flex flex-wrap items-center gap-2 text-sm pt-1">
        {/* Paste from clipboard button */}
        <button
          type="button"
          onClick={handlePaste}
          // Button styling
          className="inline-flex items-center px-2.5 py-1 bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-300 rounded-md border border-gray-300 dark:border-gray-600 text-xs font-medium shadow-sm transition-colors duration-150 focus:outline-none focus:ring-1 focus:ring-indigo-500"
          title="Paste from clipboard"
        >
          <HiOutlineClipboard className="mr-1 h-4 w-4" />
          <span>Paste</span>
        </button>

        {/* Examples button and dropdown */}
        {showExamples && (
          <div className="relative" ref={examplesRef}>
            <button
              type="button"
              onClick={() => {
                setShowExampleList(!showExampleList);
                setShowRecentList(false); // Close other dropdown
              }}
              aria-haspopup="true"
              aria-expanded={showExampleList}
              aria-controls="examples-dropdown"
              // Button styling
              className="inline-flex items-center px-2.5 py-1 bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-300 rounded-md border border-gray-300 dark:border-gray-600 text-xs font-medium shadow-sm transition-colors duration-150 focus:outline-none focus:ring-1 focus:ring-indigo-500"
            >
              <HiOutlineLightBulb className="mr-1 h-4 w-4" />
              <span>Examples</span>
            </button>

            {/* Examples Dropdown */}
            {showExampleList && (
              <div
                id="examples-dropdown"
                // Dropdown styling
                className="absolute z-20 mt-1 w-64 bg-white dark:bg-gray-800 border border-gray-200 dark:border-gray-700 rounded-md shadow-lg ring-1 ring-black ring-opacity-5 focus:outline-none"
                role="menu"
                aria-orientation="vertical"
                aria-labelledby="examples-button" // Needs id on button
              >
                <ul className="py-1 max-h-64 overflow-y-auto" role="none">
                  {EXAMPLE_MOLECULES.map((molecule, index) => (
                    <DropdownItem key={index} molecule={molecule} onClick={handleExampleClick} />
                  ))}
                </ul>
              </div>
            )}
          </div>
        )}

        {/* Recent molecules button and dropdown */}
        {showRecent && recentMolecules.length > 0 && (
          <div className="relative" ref={recentRef}>
            <button
              type="button"
              onClick={() => {
                setShowRecentList(!showRecentList);
                setShowExampleList(false); // Close other dropdown
              }}
              aria-haspopup="true"
              aria-expanded={showRecentList}
              aria-controls="recent-dropdown"
              // Button styling
              className="inline-flex items-center px-2.5 py-1 bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-300 rounded-md border border-gray-300 dark:border-gray-600 text-xs font-medium shadow-sm transition-colors duration-150 focus:outline-none focus:ring-1 focus:ring-indigo-500"
            >
              <HiOutlineClock className="mr-1 h-4 w-4" /> {/* Replaced SVG */}
              <span>Recent</span>
            </button>

            {/* Recent Dropdown */}
            {showRecentList && (
              <div
                id="recent-dropdown"
                // Dropdown styling
                className="absolute z-20 mt-1 w-64 bg-white dark:bg-gray-800 border border-gray-200 dark:border-gray-700 rounded-md shadow-lg ring-1 ring-black ring-opacity-5 focus:outline-none"
                role="menu"
                aria-orientation="vertical"
                aria-labelledby="recent-button" // Needs id on button
              >
                <ul className="py-1 max-h-64 overflow-y-auto" role="none">
                  {recentMolecules.map((molecule, index) => (
                    <DropdownItem key={molecule.timestamp || index} molecule={molecule} onClick={handleRecentClick} />
                  ))}
                </ul>
              </div>
            )}
          </div>
        )}

        {/* Hint text */}
        <p className="text-xs text-gray-500 dark:text-gray-400 pl-1">
          Enter a valid SMILES string.
        </p>
      </div>
    </div>
  );
};

export default SMILESInput;
