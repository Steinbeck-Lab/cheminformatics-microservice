import React, { useState, useRef, useMemo } from 'react';
import {
  HiOutlineClipboard,
  HiOutlineDownload,
  HiOutlineCheckCircle,
  HiOutlineX
} from 'react-icons/hi';
import { useAppContext } from '../../context/AppContext';

// Default fallback SVG
const FALLBACK_SVG_BASE64 = 'data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48dGV4dCB4PSI1MCIgeT0iNTAiIGRvbWluYW50LWJhc2VsaW5lPSJtaWRkbGUiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTAiIGZpbGw9IiM4ODg4ODg4Ij5FcnJvciBsb2FkaW5nIHN0cnVjdHVyZTwvdGV4dD48L3N2Zz4=';

const HighlightedMoleculeCard = ({
  smiles,
  title = 'Molecule',
  description,
  functionalGroups = [],
  highlightedGroupIndex = null, // Index of the functional group to highlight
  showActions = true,
  size = 'xl',
  onClick = null
}) => {
  const [copied, setCopied] = useState(false);
  const [showFullSmiles, setShowFullSmiles] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const copyTextRef = useRef(null);

  // Context handling
  const appContext = useAppContext();
  let addRecentMolecule = () => {};
  if (appContext && typeof appContext.addRecentMolecule === 'function') {
    addRecentMolecule = appContext.addRecentMolecule;
  }

  // Configuration for image container height based on size
  const sizeClasses = {
    sm: 'h-32',
    md: 'h-48',
    lg: 'h-64',
    xl: 'h-80'
  };

  const baseUrl = 'https://dev.api.naturalproducts.net';

  // Generate highlighting SMARTS pattern from functional groups
  const highlightPattern = useMemo(() => {
    if (highlightedGroupIndex === null || !functionalGroups[highlightedGroupIndex]) {
      return '';
    }

    const group = functionalGroups[highlightedGroupIndex];
    if (group.None) return '';

    // Try to use atom IDs to create a highlighting pattern
    if (group.atomIds && group.atomIds.length > 0) {
      // For now, we'll use a simple pattern based on the group type
      // This could be enhanced with more sophisticated SMARTS generation
      if (group.type) {
        // Use the type as a SMARTS pattern if it looks like one
        const typeStr = group.type.toString();
        if (typeStr.includes('[') || typeStr.includes('(')) {
          return typeStr;
        }
      }
      
      // Fallback patterns for common functional groups
      const commonPatterns = {
        'alcohol': '[OH]',
        'hydroxyl': '[OH]',
        'carbonyl': '[C]=[O]',
        'carboxyl': '[C](=[O])[OH]',
        'carboxylic': '[C](=[O])[OH]',
        'amino': '[NH2]',
        'amine': '[N]',
        'amide': '[C](=[O])[N]',
        'ester': '[C](=[O])[O][C]',
        'ether': '[O]([C])[C]',
        'ketone': '[C](=[O])[C]',
        'aldehyde': '[CHO]',
        'phenyl': 'c1ccccc1',
        'aromatic': 'c',
        'benzene': 'c1ccccc1',
        'nitro': '[N+](=O)[O-]',
        'sulfur': '[S]',
        'phosphorus': '[P]',
        'halogen': '[F,Cl,Br,I]',
        'fluorine': '[F]',
        'chlorine': '[Cl]',
        'bromine': '[Br]',
        'iodine': '[I]'
      };

      // Try to match group description to common patterns
      const description = (group.description || group.type || group.atoms || '').toLowerCase();
      for (const [name, pattern] of Object.entries(commonPatterns)) {
        if (description.includes(name)) {
          return pattern;
        }
      }

      // If we have a simple atom type, try to match it
      if (group.atoms) {
        const atomStr = group.atoms.toLowerCase();
        if (atomStr === 'o') return '[O]';
        if (atomStr === 'n') return '[N]';
        if (atomStr === 's') return '[S]';
        if (atomStr === 'p') return '[P]';
        if (atomStr === 'c') return '[C]';
      }
    }

    return '';
  }, [functionalGroups, highlightedGroupIndex]);

  // Generate image URL with highlighting
  const imageUrl = useMemo(() => {
    if (!smiles) return FALLBACK_SVG_BASE64;

    const params = new URLSearchParams({
      smiles: smiles,
      width: '512',
      height: '512',
      toolkit: 'cdk'
    });

    // Add highlighting if we have a pattern
    if (highlightPattern) {
      params.set('highlight', highlightPattern);
    }

    return `${baseUrl}/latest/depict/2D?${params.toString()}`;
  }, [smiles, highlightPattern]);

  // Copy functionality
  const handleCopy = (e) => {
    e.stopPropagation();
    
    const smilesToCopy = typeof smiles === 'string' ? smiles : '';
    if (!smilesToCopy) {
      console.error('No SMILES to copy');
      return;
    }
    
    try {
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(smilesToCopy)
          .then(() => {
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
          })
          .catch(() => {
            fallbackCopy();
          });
        return;
      }
      
      fallbackCopy();
      
    } catch (err) {
      console.error('Copy failed with all methods:', err);
      setShowCopyModal(true);
    }
  };

  const fallbackCopy = () => {
    try {
      const textArea = document.createElement('textarea');
      textArea.value = smiles;
      textArea.style.position = 'fixed';
      textArea.style.left = '-999999px';
      textArea.style.top = '-999999px';
      document.body.appendChild(textArea);
      textArea.select();
      
      const successful = document.execCommand('copy');
      document.body.removeChild(textArea);
      
      if (successful) {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
        return;
      } else {
        throw new Error('execCommand copy failed');
      }
    } catch (err) {
      console.debug('Fallback copy failed:', err);
      setShowCopyModal(true);
    }
  };

  // Auto-select text in copy modal
  React.useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      setTimeout(() => {
        copyTextRef.current.select();
      }, 100);
    }
  }, [showCopyModal]);

  // Handle downloading the structure
  const handleDownload = (e) => {
    e.stopPropagation();
    if (typeof smiles !== 'string' || !smiles.trim()) {
      console.error("Cannot download: Invalid SMILES string provided.");
      return;
    }

    const molUrl = `${baseUrl}/latest/convert/mol2D?smiles=${encodeURIComponent(smiles)}&toolkit=cdk`;
    
    fetch(molUrl)
      .then(response => {
        if (!response.ok) {
          throw new Error(`Failed to download structure: ${response.statusText}`);
        }
        return response.text();
      })
      .then(molData => {
        const blob = new Blob([molData], { type: 'chemical/x-mdl-molfile' });
        const url = URL.createObjectURL(blob);
        
        const a = document.createElement('a');
        a.href = url;
        const filename = `${title.replace(/[^a-z0-9]/gi, '_').toLowerCase() || 'molecule'}.sdf`;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
      })
      .catch(err => {
        console.error('Structure download failed:', err);
        alert('Failed to download structure. Please try again later.');
      });
  };

  // Handle card click
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

  // Format SMILES for display
  const MAX_SMILES_LENGTH = 35;
  const formattedSmiles = showFullSmiles ? smiles : (
    smiles && smiles.length > MAX_SMILES_LENGTH ? `${smiles.substring(0, MAX_SMILES_LENGTH)}...` : smiles
  );

  return (
    <>
      <div
        className={`bg-white dark:bg-gray-800 rounded-lg overflow-hidden shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700
          hover:border-blue-500 dark:hover:border-blue-500 hover:shadow-lg transition-all duration-200 group ${onClick ? 'cursor-pointer focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 dark:focus:ring-offset-gray-800' : ''}`}
        onClick={onClick ? handleCardClick : undefined}
        onKeyPress={onClick ? (e) => { if (e.key === 'Enter' || e.key === ' ') handleCardClick(); } : undefined}
        role={onClick ? 'button' : 'figure'}
        tabIndex={onClick ? 0 : undefined}
        aria-label={title}
      >
        {/* Title bar */}
        <div className="px-4 py-2 bg-gray-50 dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700 flex justify-between items-center">
          <h3 className="font-medium text-sm text-gray-700 dark:text-gray-200 truncate" title={title}>
            {title}
          </h3>
          {highlightedGroupIndex !== null && functionalGroups[highlightedGroupIndex] && !functionalGroups[highlightedGroupIndex].None && (
            <span className="text-xs bg-blue-100 dark:bg-blue-900 text-blue-800 dark:text-blue-200 px-2 py-1 rounded-full">
              Highlighted
            </span>
          )}
        </div>

        {/* Molecule visualization area */}
        <div className={`relative p-2 flex justify-center items-center bg-gray-50 dark:bg-gray-800 ${sizeClasses[size] || 'h-60'}`}>
          {smiles ? (
            <img
              key={imageUrl}
              src={imageUrl}
              alt={`Chemical structure of ${title}`}
              className="max-h-full max-w-full object-contain"
              onError={(e) => {
                console.error(`Error loading image for ${title} from ${imageUrl}:`, e);
                e.target.onerror = null;
                e.target.src = FALLBACK_SVG_BASE64;
                e.target.classList.add('p-4');
              }}
            />
          ) : (
            <div className="text-center text-gray-400 dark:text-gray-500 text-xs p-4">
              No structure to display.
            </div>
          )}
        </div>

        {/* SMILES information section */}
        <div className="px-4 py-3 border-t border-gray-200 dark:border-gray-700 text-sm">
          <div className="flex items-start">
            <span className="font-medium text-gray-500 dark:text-gray-400 mr-2 flex-shrink-0">SMILES:</span>
            <div className="text-gray-700 dark:text-gray-300 flex-1 min-w-0">
              <div
                className="inline-block cursor-pointer group/smiles"
                onClick={(e) => {
                  e.stopPropagation();
                  setShowFullSmiles(!showFullSmiles);
                }}
                title={showFullSmiles ? "Click to shorten" : "Click to expand"}
              >
                <span className="font-mono text-xs break-all group-hover/smiles:text-blue-600 dark:group-hover/smiles:text-blue-400">
                  {formattedSmiles || '-'}
                </span>
                {smiles && !showFullSmiles && smiles.length > MAX_SMILES_LENGTH && (
                  <span className="ml-1 text-blue-500 dark:text-blue-400 text-xs font-medium opacity-70 group-hover/smiles:opacity-100">
                    more
                  </span>
                )}
                {smiles && showFullSmiles && smiles.length > MAX_SMILES_LENGTH && (
                  <span className="ml-1 text-blue-500 dark:text-blue-400 text-xs font-medium opacity-70 group-hover/smiles:opacity-100">
                    less
                  </span>
                )}
              </div>
            </div>
          </div>
        </div>

        {/* Optional description */}
        {description && (
          <div className="px-4 pb-3 text-gray-500 dark:text-gray-400 text-xs">
            {description}
          </div>
        )}

        {/* Action buttons */}
        {showActions && (
          <div className="px-3 py-1.5 bg-gray-50 dark:bg-gray-800 border-t border-gray-200 dark:border-gray-700 flex justify-end space-x-1">
            {/* Copy Button */}
            <button
              onClick={handleCopy}
              disabled={!smiles}
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
              disabled={!smiles}
              className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 ${!smiles ? 'text-gray-300 dark:text-gray-600 cursor-not-allowed' :
                  'text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700'
                }`}
              title="Download structure as SDF"
              aria-label="Download structure as SDF"
            >
              <HiOutlineDownload className="h-5 w-5" />
            </button>
          </div>
        )}
      </div>

      {/* Manual Copy Modal */}
      {showCopyModal && (
        <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4">
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-2xl max-w-lg w-full">
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-bold text-gray-900 dark:text-white">Copy SMILES</h3>
              <button 
                onClick={() => setShowCopyModal(false)}
                className="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200"
              >
                <HiOutlineX className="w-5 h-5" />
              </button>
            </div>
            <p className="mb-4 text-gray-700 dark:text-gray-300">
              Automatic copying failed. Please select and copy this text manually:
            </p>
            <div className="mb-4">
              <input
                type="text"
                ref={copyTextRef}
                value={smiles || ''}
                readOnly
                className="w-full p-2 border border-gray-300 dark:border-gray-600 rounded font-mono text-sm bg-gray-50 dark:bg-gray-900"
                onClick={e => e.target.select()}
              />
            </div>
            <div className="flex justify-end gap-3">
              <button
                onClick={() => setShowCopyModal(false)}
                className="px-4 py-2 bg-gray-200 dark:bg-gray-700 text-gray-800 dark:text-gray-200 rounded hover:bg-gray-300 dark:hover:bg-gray-600"
              >
                Close
              </button>
            </div>
          </div>
        </div>
      )}
    </>
  );
};

export default HighlightedMoleculeCard;
