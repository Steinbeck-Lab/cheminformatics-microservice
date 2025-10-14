// Description: A React component for inputting SMILES notation with example and recent molecule dropdowns.
import React, { useState, useRef, useEffect } from "react";
// Import necessary icons
import {
  HiOutlineX,
  HiOutlineLightBulb,
  HiOutlineClipboard,
  HiOutlineClock,
  HiOutlineExclamationCircle, // Added for paste modal
} from "react-icons/hi";
import { useAppContext } from "../../context/AppContext"; // Assuming context provides recentMolecules

// Example molecules data (static, moved outside component)
const EXAMPLE_MOLECULES = [
  {
    name: "Caffeine",
    smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    description: "Stimulant found in coffee",
  },
  {
    name: "Aspirin",
    smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
    description: "Analgesic & anti-inflammatory",
  },
  {
    name: "Sucrose",
    smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O",
    description: "Common sugar",
  },
  {
    name: "Cholesterol",
    smiles: "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
    description: "Steroid in cell membranes",
  },
  {
    name: "Paracetamol",
    smiles: "CC(=O)NC1=CC=C(C=C1)O",
    description: "Pain reliever & fever reducer",
  },
  {
    name: "Ibuprofen",
    smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    description: "Anti-inflammatory drug",
  },
];

// Dropdown List Item Component (for cleaner code)
const DropdownItem = ({ molecule, onClick }) => (
  <li
    className="px-3 py-2 hover:bg-gray-100 dark:hover:bg-gray-700 cursor-pointer rounded-md mx-1"
    onClick={() => onClick(molecule.smiles)}
    role="menuitem"
    tabIndex={-1} // Can be focused programmatically if needed
  >
    <p
      className="font-medium text-sm text-gray-800 dark:text-gray-100 truncate"
      title={molecule.name || `SMILES: ${molecule.smiles}`}
    >
      {molecule.name || "Recent Molecule"}
    </p>
    <p className="text-xs text-gray-500 dark:text-gray-400 truncate" title={molecule.smiles}>
      {molecule.smiles}
    </p>
    {molecule.description && (
      <p className="text-xs text-gray-400 dark:text-gray-500 mt-1">{molecule.description}</p>
    )}
  </li>
);

// Paste Modal Component
const PasteModal = ({ isOpen, onClose, onPaste }) => {
  const [text, setText] = useState("");
  const textareaRef = useRef(null);

  useEffect(() => {
    if (isOpen && textareaRef.current) {
      textareaRef.current.focus();
    }
  }, [isOpen]);

  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4">
      <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-2xl max-w-lg w-full animate-fadeIn">
        <div className="flex justify-between items-center mb-4">
          <h3 className="text-lg font-bold text-gray-900 dark:text-white flex items-center">
            <HiOutlineClipboard className="h-5 w-5 mr-2 text-indigo-500 dark:text-indigo-400" />
            Paste SMILES
          </h3>
          <button
            onClick={onClose}
            className="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200"
            aria-label="Close"
          >
            <HiOutlineX className="w-5 h-5" />
          </button>
        </div>
        <p className="mb-4 text-gray-700 dark:text-gray-300">
          Please paste your SMILES string in the area below:
        </p>
        <div className="mb-4">
          <textarea
            ref={textareaRef}
            value={text}
            onChange={(e) => setText(e.target.value)}
            className="w-full p-2 border border-gray-300 dark:border-gray-600 rounded font-mono text-sm bg-gray-50 dark:bg-gray-900 text-gray-900 dark:text-gray-100 h-24"
            placeholder="Paste SMILES here..."
            autoFocus
          />
        </div>
        <div className="flex justify-end gap-3">
          <button
            onClick={onClose}
            className="px-4 py-2 bg-gray-200 dark:bg-gray-700 text-gray-800 dark:text-gray-200 rounded hover:bg-gray-300 dark:hover:bg-gray-600"
          >
            Cancel
          </button>
          <button
            onClick={() => {
              onPaste(text);
              onClose();
            }}
            className="px-4 py-2 bg-indigo-600 text-white rounded hover:bg-indigo-700"
            disabled={!text.trim()}
          >
            Use SMILES
          </button>
        </div>
      </div>
    </div>
  );
};

// Add custom styles for animations
if (typeof document !== "undefined") {
  const styleSheet = document.createElement("style");
  styleSheet.type = "text/css";
  styleSheet.textContent = `
    @keyframes fadeIn {
      from { opacity: 0; transform: translateY(10px); }
      to { opacity: 1; transform: translateY(0); }
    }
    .animate-fadeIn {
      animation: fadeIn 0.3s ease-out forwards;
    }
  `;
  document.head.appendChild(styleSheet);
}

const SMILESInput = ({
  value,
  onChange,
  placeholder = "Enter SMILES notation...",
  showExamples = true,
  showRecent = true,
  label = "SMILES",
  required = false,
}) => {
  const [showExampleList, setShowExampleList] = useState(false);
  const [showRecentList, setShowRecentList] = useState(false);
  const [showPasteModal, setShowPasteModal] = useState(false);
  const [pasteError, setPasteError] = useState("");

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
    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, []);

  const handleClear = () => {
    onChange("");
    setPasteError("");
  };

  const handleExampleClick = (smiles) => {
    onChange(smiles);
    setShowExampleList(false); // Close dropdown after selection
    setPasteError("");
  };

  const handleRecentClick = (smiles) => {
    onChange(smiles);
    setShowRecentList(false); // Close dropdown after selection
    setPasteError("");
  };

  // Enhanced paste functionality with multiple approaches
  const handlePaste = async () => {
    setPasteError("");

    // Method 1: Try using navigator.clipboard API
    if (navigator.clipboard && navigator.clipboard.readText) {
      try {
        const text = await navigator.clipboard.readText();
        if (text) {
          onChange(text);
          return;
        }
      } catch (err) {
        console.debug("Clipboard readText failed:", err);
        // Fall through to alternatives
      }
    }

    // Method 2: Try using document.execCommand
    try {
      // Create a hidden textarea
      const textArea = document.createElement("textarea");
      textArea.style.position = "fixed";
      textArea.style.top = "-999px";
      textArea.style.left = "-999px";
      document.body.appendChild(textArea);

      // Focus and select content (may trigger paste permission in some browsers)
      textArea.focus();
      document.execCommand("paste");

      // Get pasted content
      const pastedText = textArea.value;
      document.body.removeChild(textArea);

      if (pastedText) {
        onChange(pastedText);
        return;
      }
    } catch (err) {
      console.debug("execCommand paste failed:", err);
      // Fall through to manual method
    }

    // Method 3: Show modal for manual paste
    setShowPasteModal(true);
  };

  // Handle paste from the modal
  const handleModalPaste = (text) => {
    if (text && text.trim()) {
      onChange(text.trim());
    }
  };

  return (
    <div className="space-y-2">
      {/* Label */}
      {label && (
        <label
          htmlFor="smiles-input"
          className="block text-sm font-medium text-gray-700 dark:text-gray-300"
        >
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

      {/* Paste Error Display */}
      {pasteError && (
        <div className="text-sm text-red-600 dark:text-red-400 flex items-start">
          <HiOutlineExclamationCircle className="h-4 w-4 mt-0.5 mr-1" />
          <span>{pasteError}</span>
        </div>
      )}

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
                    <DropdownItem
                      key={molecule.timestamp || index}
                      molecule={molecule}
                      onClick={handleRecentClick}
                    />
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

      {/* Manual Paste Modal */}
      <PasteModal
        isOpen={showPasteModal}
        onClose={() => setShowPasteModal(false)}
        onPaste={handleModalPaste}
      />
    </div>
  );
};

export default SMILESInput;
