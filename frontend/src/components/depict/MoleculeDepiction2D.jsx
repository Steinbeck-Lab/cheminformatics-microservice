// Description: A component for 2D molecule depiction
import React, { useState, useEffect, useCallback } from 'react';
import {
  HiOutlineDownload,
  HiOutlineRefresh
} from 'react-icons/hi';

// Assuming this service is configured correctly
import depictService from '../../services/depictService';

const MoleculeDepiction2D = ({ smiles, title, toolkit = 'rdkit' }) => {
  const [error, setError] = useState(null);
  const [imageUrl, setImageUrl] = useState('');
  const [rotation, setRotation] = useState(0);
  const [useUnicolor, setUseUnicolor] = useState(false);
  const [showCIP, setShowCIP] = useState(false);
  
  // Fixed width and height - not using state to avoid ESLint warnings
  const imageWidth = 400;
  const imageHeight = 350;

  // Generate the depiction URL - using useCallback for proper memoization
  const generateDepiction = useCallback(() => {
    if (!smiles || !depictService || typeof depictService.get2DDepictionUrl !== 'function') {
      setError("Cannot generate 2D depiction: missing SMILES or service unavailable");
      return;
    }

    try {
      // Prepare options
      const options = {
        toolkit,
        width: imageWidth,
        height: imageHeight,
        rotate: rotation,
        CIP: toolkit === 'cdk' ? showCIP : undefined,
        unicolor: useUnicolor
      };

      // Get the URL
      const url = depictService.get2DDepictionUrl(smiles, options);
      setImageUrl(url);
      setError(null);
    } catch (err) {
      console.error("Error generating 2D depiction:", err);
      setError(`Failed to generate 2D depiction: ${err.message}`);
    }
  }, [smiles, toolkit, rotation, useUnicolor, showCIP]);

  // Generate depiction when props or settings change
  useEffect(() => {
    if (smiles) {
      generateDepiction();
    }
  }, [smiles, generateDepiction]);

  // Handle rotation change
  const handleRotationChange = (e) => {
    setRotation(Number(e.target.value));
  };

  // Handle refresh button click - reset to original settings
  const handleRefresh = () => {
    // Reset all settings to their initial values
    setRotation(0);
    setUseUnicolor(false);
    setShowCIP(false);
    // The image will update due to dependency changes
  };

  // Handle download
  const handleDownload = () => {
    if (!imageUrl) return;
    
    // Determine file extension based on current URL
    const extension = imageUrl.includes('format=png') ? 'png' : 'svg';
    
    // Create filename from title or SMILES
    const filenameSafe = (title || 'molecule').replace(/[^a-z0-9]/gi, '_').toLowerCase();
    
    // Create the download link
    const a = document.createElement('a');
    a.href = imageUrl;
    a.download = `${filenameSafe}_2d.${extension}`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700 overflow-hidden h-full flex flex-col">
      {/* Header */}
      <div className="px-4 py-3 border-b border-gray-200 dark:border-gray-700 flex justify-between items-center bg-gradient-to-r from-gray-50 to-gray-100 dark:from-gray-800 dark:to-gray-750">
        <h3 className="font-medium text-gray-800 dark:text-white text-sm sm:text-base truncate">
          {title || "2D Structure"}
        </h3>
        
        {/* Download Button */}
        <button
          onClick={handleDownload}
          disabled={!imageUrl}
          className="p-1.5 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500"
          title="Download 2D structure"
        >
          <HiOutlineDownload className="h-5 w-5" />
        </button>
      </div>
      
      {/* Image Display */}
      <div className="p-4 bg-white flex-grow flex items-center justify-center" style={{ minHeight: "350px" }}>
        {imageUrl ? (
          <img
            src={imageUrl}
            alt={title || "2D Molecular Structure"}
            className="max-w-full max-h-full object-contain"
            onError={() => setError("Failed to load 2D depiction")}
          />
        ) : (
          <div className="text-center text-gray-500 dark:text-gray-400">
            {error || "Loading depiction..."}
          </div>
        )}
      </div>
      
      {/* Controls */}
      <div className="px-4 py-3 bg-gray-50 dark:bg-gray-900 border-t border-gray-200 dark:border-gray-700">
        <div className="flex flex-wrap gap-4">
          {/* Rotation Slider */}
          <div className="flex-1 min-w-[120px]">
            <label htmlFor="rotation-slider" className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
              Rotate: {rotation}°
            </label>
            <input
              id="rotation-slider"
              type="range"
              value={rotation}
              onChange={handleRotationChange}
              min="0" max="359" step="1"
              className="w-full h-2 bg-gray-200 dark:bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500 dark:accent-blue-400"
            />
          </div>
          
          {/* Unicolor Toggle */}
          <div className="flex items-center">
            <input
              id="unicolor-toggle"
              type="checkbox"
              checked={useUnicolor}
              onChange={(e) => setUseUnicolor(e.target.checked)}
              className="h-4 w-4 rounded border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
            />
            <label htmlFor="unicolor-toggle" className="ml-2 text-xs text-gray-700 dark:text-gray-300">
              Black & White
            </label>
          </div>
          
          {/* CIP Stereo Toggle (CDK only) */}
          {toolkit === 'cdk' && (
            <div className="flex items-center">
              <input
                id="cip-toggle"
                type="checkbox"
                checked={showCIP}
                onChange={(e) => setShowCIP(e.target.checked)}
                className="h-4 w-4 rounded border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
              />
              <label htmlFor="cip-toggle" className="ml-2 text-xs text-gray-700 dark:text-gray-300">
                Show R/S
              </label>
            </div>
          )}
          
          {/* Refresh Button */}
          <button
            onClick={handleRefresh}
            className="ml-auto p-1.5 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500"
            title="Reset to original view"
          >
            <HiOutlineRefresh className="h-5 w-5" />
          </button>
        </div>
      </div>
    </div>
  );
};

export default MoleculeDepiction2D;