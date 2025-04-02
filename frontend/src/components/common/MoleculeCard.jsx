// Description: A React component for displaying a molecule card with SMILES representation, image, and action buttons.
import React, { useState } from 'react';
import {
  HiOutlineClipboard,
  HiOutlineDownload,
  HiOutlineInformationCircle,
  HiOutlineCheckCircle // Used for copy success indication
} from 'react-icons/hi';
import { useAppContext } from '../../context/AppContext'; // Assuming context provides addRecentMolecule

// Default fallback SVG with a neutral gray fill color (#888888)
const FALLBACK_SVG_BASE64 = 'data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48dGV4dCB4PSI1MCIgeT0iNTAiIGRvbWluYW50LWJhc2VsaW5lPSJtaWRkbGUiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTAiIGZpbGw9IiM4ODg4ODg4Ij5FcnJvciBsb2FkaW5nIHN0cnVjdHVyZTwvdGV4dD48L3N2Zz4=';


const MoleculeCard = ({
  smiles,
  title = 'Molecule',
  description,
  showActions = true,
  size = 'xl', // sm, md, lg, xl
  onClick = null // Optional click handler for the whole card
}) => {
  const [copied, setCopied] = useState(false);
  const [showFullSmiles, setShowFullSmiles] = useState(false);

  // --- Context Handling ---
  const appContext = useAppContext();
  // Declare addRecentMolecule once with let, provide default no-op function
  let addRecentMolecule = () => {
    // console.warn("MoleculeCard: addRecentMolecule context function not available or not used.");
  };

  if (!appContext) {
    console.error("MoleculeCard: AppContext is not available.");
    // addRecentMolecule remains the default no-op function
  } else {
    // Assign the function from context if it exists
    if (typeof appContext.addRecentMolecule === 'function') {
      addRecentMolecule = appContext.addRecentMolecule;
    } else {
      console.warn("MoleculeCard: addRecentMolecule function not found in AppContext.");
      // addRecentMolecule remains the default no-op function
    }
  }
  // --- End Context Handling ---


  // Configuration for image container height based on size prop
  const sizeClasses = {
    sm: 'h-32',
    md: 'h-48',
    lg: 'h-64',
    xl: 'h-80'
  };

  // --- Reverted API URL Construction ---
  // Using the hardcoded base URL that was likely working in your original code.
  const baseUrl = 'https://api.naturalproducts.net'; // Hardcoded base URL
  const imageUrl = smiles // Only generate URL if SMILES is valid
    ? `${baseUrl}/latest/depict/2D?smiles=${encodeURIComponent(smiles)}&width=512&height=512&toolkit=cdk` // Reverted to original structure
    : FALLBACK_SVG_BASE64; // Use fallback if no SMILES provided initially

  // Handle copying the SMILES string
  const handleCopy = (e) => {
    e.stopPropagation(); // Prevent card click if clicking button
    if (!navigator.clipboard) {
      console.error('Clipboard API not available');
      return;
    }
    const smilesToCopy = typeof smiles === 'string' ? smiles : '';
    navigator.clipboard.writeText(smilesToCopy)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      })
      .catch(err => console.error('Failed to copy SMILES:', err));
  };

  // Handle downloading the structure as SVG
  const handleDownload = (e) => {
    e.stopPropagation(); // Prevent card click
    if (typeof smiles !== 'string' || !smiles.trim()) {
      console.error("Cannot download: Invalid SMILES string provided.");
      return;
    }
    // Prefer a dedicated SVG endpoint if available (using the same base URL)
    const svgUrl = `${baseUrl}/latest/depict/svg?smiles=${encodeURIComponent(smiles)}`;

    fetch(svgUrl)
      .then(response => {
        if (!response.ok) {
          console.warn(`SVG download failed (${response.statusText}), attempting image download from: ${imageUrl}`);
          throw new Error('SVG endpoint failed');
        }
        return response.text();
      })
      .then(svgText => {
        const blob = new Blob([svgText], { type: 'image/svg+xml;charset=utf-8' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        const filename = `${title.replace(/[^a-z0-9]/gi, '_').toLowerCase() || 'molecule'}.svg`;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
      })
      .catch(err => {
        console.error('SVG download failed, falling back to image URL:', err);
        const fallbackUrl = imageUrl;
        const a = document.createElement('a');
        a.href = fallbackUrl;
        a.download = `${title.replace(/[^a-z0-9]/gi, '_').toLowerCase() || 'molecule'}.png`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
      });
  };

  // Handle showing molecule details (placeholder)
  const handleShowDetails = (e) => {
    e.stopPropagation();
    console.log('Placeholder: Show molecule details for:', smiles);
  };


  // Handle click on the card itself (if onClick prop is provided)
  const handleCardClick = () => {
    if (onClick) {
      onClick(smiles);
    }
    addRecentMolecule({
      smiles,
      name: title !== 'Molecule' ? title : null,
      timestamp: new Date().toISOString()
    });
  };

  // Format SMILES for display (shorten if too long)
  const MAX_SMILES_LENGTH = 35;
  const formattedSmiles = showFullSmiles ? smiles : (
    smiles && smiles.length > MAX_SMILES_LENGTH ? `${smiles.substring(0, MAX_SMILES_LENGTH)}...` : smiles
  );

  return (
    <div
      // Card container styling with theme awareness
      className={`bg-white dark:bg-gray-800 rounded-lg overflow-hidden shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700
        hover:border-blue-500 dark:hover:border-blue-500 hover:shadow-lg transition-all duration-200 group ${onClick ? 'cursor-pointer focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 dark:focus:ring-offset-gray-800' : ''}`}
      onClick={onClick ? handleCardClick : undefined}
      onKeyPress={onClick ? (e) => { if (e.key === 'Enter' || e.key === ' ') handleCardClick(); } : undefined} // Keyboard accessibility
      role={onClick ? 'button' : 'figure'} // Role depends on interactivity
      tabIndex={onClick ? 0 : undefined} // Make clickable cards focusable
      aria-label={title}
    >
      {/* Title bar */}
      <div className="px-4 py-2 bg-gray-50 dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700 flex justify-between items-center">
        {/* Theme-aware title text */}
        <h3 className="font-medium text-sm text-gray-700 dark:text-gray-200 truncate" title={title}>
          {title}
        </h3>
      </div>

      {/* Molecule visualization area */}
      {/* Centering the image, theme-aware background */}
      <div className={`relative p-2 flex justify-center items-center bg-gray-50 dark:bg-gray-800 ${sizeClasses[size] || 'h-60'}`}>
        {smiles ? ( // Only render img tag if smiles is present
          <img
            key={imageUrl} // Add key to force re-render if URL changes significantly
            src={imageUrl}
            alt={`Chemical structure of ${title}`}
            // Fallback SVG text color is now hardcoded gray
            className="max-h-full max-w-full object-contain"
            // Error handler with improved fallback SVG and logging
            onError={(e) => {
              console.error(`Error loading image for ${title} from ${imageUrl}:`, e);
              e.target.onerror = null; // Prevent infinite loop if fallback fails
              e.target.src = FALLBACK_SVG_BASE64;
              e.target.classList.add('p-4'); // Add padding to fallback text
            }}
          />
        ) : (
          // Placeholder if no SMILES is provided
          <div className="text-center text-gray-400 dark:text-gray-500 text-xs p-4">
            No structure to display.
          </div>
        )}
      </div>

      {/* SMILES information section */}
      <div className="px-4 py-3 border-t border-gray-200 dark:border-gray-700 text-sm">
        <div className="flex items-start">
          {/* Label styling */}
          <span className="font-medium text-gray-500 dark:text-gray-400 mr-2 flex-shrink-0">SMILES:</span>
          {/* SMILES text container */}
          <div className="text-gray-700 dark:text-gray-300 flex-1 min-w-0"> {/* Ensure flex item can shrink */}
            <div
              className="inline-block cursor-pointer group/smiles" // Group for hover effect
              onClick={(e) => {
                e.stopPropagation(); // Prevent card click
                setShowFullSmiles(!showFullSmiles);
              }}
              title={showFullSmiles ? "Click to shorten" : "Click to expand"}
            >
              {/* SMILES text styling */}
              <span className="font-mono text-xs break-all group-hover/smiles:text-blue-600 dark:group-hover/smiles:text-blue-400">
                {formattedSmiles || '-'} {/* Show dash if smiles is empty */}
              </span>
              {/* Show more/less indicator */}
              {smiles && !showFullSmiles && smiles.length > MAX_SMILES_LENGTH && (
                <span className="ml-1 text-blue-500 dark:text-blue-400 text-xs font-medium opacity-70 group-hover/smiles:opacity-100">
                  more
                </span>
              )}
              {smiles && showFullSmiles && smiles.length > MAX_SMILES_LENGTH && ( // Only show 'less' if it was actually truncated
                <span className="ml-1 text-blue-500 dark:text-blue-400 text-xs font-medium opacity-70 group-hover/smiles:opacity-100">
                  less
                </span>
              )}
            </div>
          </div>
        </div>

        {/* Optional description */}
        {description && (
          <div className="mt-1 text-gray-500 dark:text-gray-400 text-xs">
            {description}
          </div>
        )}
      </div>

      {/* Action buttons (optional) */}
      {showActions && (
        // Actions bar styling
        <div className="px-3 py-1.5 bg-gray-50 dark:bg-gray-800 border-t border-gray-200 dark:border-gray-700 flex justify-end space-x-1">
          {/* Copy Button */}
          <button
            onClick={handleCopy}
            disabled={!smiles} // Disable if no SMILES
            // Button styling with theme awareness
            className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 ${!smiles ? 'text-gray-300 dark:text-gray-600 cursor-not-allowed' :
                copied ? 'text-green-500 dark:text-green-500' :
                  'text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700'
              }`}
            title={copied ? "Copied!" : "Copy SMILES to clipboard"}
            aria-label={copied ? "SMILES Copied" : "Copy SMILES"}
          >
            {copied ? (
              <HiOutlineCheckCircle className="h-5 w-5" />
            ) : (
              <HiOutlineClipboard className="h-5 w-5" />
            )}
          </button>

          {/* Download Button */}
          <button
            onClick={handleDownload}
            disabled={!smiles} // Disable if no SMILES
            className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 ${!smiles ? 'text-gray-300 dark:text-gray-600 cursor-not-allowed' :
                'text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700'
              }`}
            title="Download structure as SVG"
            aria-label="Download structure as SVG"
          >
            <HiOutlineDownload className="h-5 w-5" />
          </button>

          {/* Info Button (Example) */}
          <button
            onClick={handleShowDetails}
            disabled={!smiles} // Disable if no SMILES
            className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 ${!smiles ? 'text-gray-300 dark:text-gray-600 cursor-not-allowed' :
                'text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700'
              }`}
            title="View molecule details (placeholder)"
            aria-label="View molecule details"
          >
            <HiOutlineInformationCircle className="h-5 w-5" />
          </button>
        </div>
      )}
    </div>
  );
};

export default MoleculeCard;
