// Description: A React component to display SMILES strings with options to copy to clipboard and download as a MOL file.
import React, { useState } from 'react';
import { HiOutlineClipboard, HiOutlineDownload, HiOutlineCheck } from 'react-icons/hi';
import { generate2DCoordinates } from '../../services/convertService';

const SMILESDisplay = ({ smiles, label = "SMILES", showDownload = true }) => {
  const [copied, setCopied] = useState(false);
  const [downloading, setDownloading] = useState(false);

  // Handle copy to clipboard
  const handleCopy = async () => {
    if (!smiles) return;

    try {
      await navigator.clipboard.writeText(smiles);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000); // Reset after 2 seconds
    } catch (err) {
      console.error('Failed to copy text:', err);
    }
  };

  // Handle download as MOL
  const handleDownload = async () => {
    if (!smiles) return;

    setDownloading(true);
    try {
      // Convert SMILES to MOL using the service
      const molData = await generate2DCoordinates(smiles, 'rdkit');

      // Create and trigger download
      const blob = new Blob([molData], { type: 'chemical/x-mdl-molfile' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'molecule.mol';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    } catch (err) {
      console.error('Failed to convert and download MOL:', err);
    } finally {
      setDownloading(false);
    }
  };

  return (
    <div className="w-full">
      <div className="flex justify-between items-center mb-1">
        <label className="block text-sm font-medium text-gray-700 dark:text-gray-300">
          {label}
        </label>
        <div className="flex space-x-2">
          <button
            onClick={handleCopy}
            disabled={!smiles}
            className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 dark:focus:ring-blue-500 ${!smiles ? 'text-gray-400 dark:text-gray-600 cursor-not-allowed' : 'text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-200 hover:bg-gray-100 dark:hover:bg-gray-700'
              }`}
            title="Copy SMILES to clipboard"
          >
            {copied ? (
              <HiOutlineCheck className="h-5 w-5 text-green-500 dark:text-green-400" />
            ) : (
              <HiOutlineClipboard className="h-5 w-5" />
            )}
          </button>
          {showDownload && (
            <button
              onClick={handleDownload}
              disabled={downloading || !smiles}
              className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-2 focus:ring-blue-500 dark:focus:ring-blue-500 ${downloading || !smiles ? 'text-gray-400 dark:text-gray-600 cursor-not-allowed' : 'text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-200 hover:bg-gray-100 dark:hover:bg-gray-700'
                }`}
              title="Download as MOL file"
            >
              <HiOutlineDownload className="h-5 w-5" />
            </button>
          )}
        </div>
      </div>
      <div className="w-full px-3 py-2 bg-gray-50 dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-gray-100 font-mono text-sm break-all">
        {smiles || "No SMILES data available"}
      </div>
    </div>
  );
};

export default SMILESDisplay;