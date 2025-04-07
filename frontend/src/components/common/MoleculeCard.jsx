import React, { useState, useRef } from 'react';
import {
  HiOutlineClipboard,
  HiOutlineDownload,
  HiOutlineInformationCircle,
  HiOutlineCheckCircle, // Used for copy success indication
  HiOutlineExclamationCircle, // Used for error display
  HiOutlineX         // Used for closing modals
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
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [showDetailsModal, setShowDetailsModal] = useState(false);
  const [moleculeDetails, setMoleculeDetails] = useState(null);
  const [isLoadingDetails, setIsLoadingDetails] = useState(false);
  const copyTextRef = useRef(null);

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
  const baseUrl = 'https://dev.api.naturalproducts.net'; // Hardcoded base URL
  const imageUrl = smiles // Only generate URL if SMILES is valid
    ? `${baseUrl}/latest/depict/2D?smiles=${encodeURIComponent(smiles)}&width=512&height=512&toolkit=cdk` // Reverted to original structure
    : FALLBACK_SVG_BASE64; // Use fallback if no SMILES provided initially

  // Enhanced handleCopy function with multiple fallback methods for production environments
  const handleCopy = (e) => {
    e.stopPropagation(); // Prevent card click if clicking button
    
    const smilesToCopy = typeof smiles === 'string' ? smiles : '';
    if (!smilesToCopy) {
      console.error('No SMILES to copy');
      return;
    }
    
    // Try multiple clipboard copy methods in sequence
    try {
      // Method 1: Use the Clipboard API (modern browsers)
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(smilesToCopy)
          .then(() => {
            setCopied(true);
            setTimeout(() => setCopied(false), 2000);
          })
          .catch(err => {
            console.debug('Clipboard API failed, trying alternative method', err);
            fallbackCopy();
          });
        return;
      }
      
      // If Clipboard API not available, try fallback immediately
      fallbackCopy();
      
    } catch (err) {
      console.error('Copy failed with all methods:', err);
      // Show manual copy modal as last resort
      setShowCopyModal(true);
    }
  };
  
  // Fallback copy method using execCommand
  const fallbackCopy = () => {
    try {
      // Method 2: Use execCommand (older browsers)
      const textArea = document.createElement('textarea');
      textArea.value = typeof smiles === 'string' ? smiles : '';
      
      // Make the textarea out of viewport
      textArea.style.position = 'fixed';
      textArea.style.left = '-999999px';
      textArea.style.top = '-999999px';
      document.body.appendChild(textArea);
      
      // Select and copy
      textArea.focus();
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
      // Show manual copy modal as last resort
      setShowCopyModal(true);
    }
  };

  // Auto-select text in copy modal when it appears
  React.useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      setTimeout(() => {
        copyTextRef.current.select();
      }, 100);
    }
  }, [showCopyModal]);

  // Handle downloading the structure as Mol
  const handleDownload = (e) => {
    e.stopPropagation(); // Prevent card click
    if (typeof smiles !== 'string' || !smiles.trim()) {
      console.error("Cannot download: Invalid SMILES string provided.");
      return;
    }
    // Prefer a dedicated mol2D endpoint if available (using the same base URL)
    // Use the molfile endpoint to download structure as SDF (Mol) format
    const molUrl = `${baseUrl}/latest/convert/mol2D?smiles=${encodeURIComponent(smiles)}&toolkit=cdk`;
    
    fetch(molUrl)
      .then(response => {
      if (!response.ok) {
        throw new Error(`Failed to download structure: ${response.statusText}`);
      }
      return response.text();
      })
      .then(molData => {
      // Create blob with correct MIME type for mol/SDF file
      const blob = new Blob([molData], { type: 'chemical/x-mdl-molfile' });
      const url = URL.createObjectURL(blob);
      
      // Create and trigger download link
      const a = document.createElement('a');
      a.href = url;
      // Use SDF extension for the correct file format
      const filename = `${title.replace(/[^a-z0-9]/gi, '_').toLowerCase() || 'molecule'}.sdf`;
      a.download = filename;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
      })
      .catch(err => {
      console.error('Structure download failed:', err);
      // Alert user of failure
      alert('Failed to download structure. Please try again later.');
      });
    };


  // Fetch and display molecule details using the descriptor endpoint
  const fetchMoleculeDetails = async (smiles) => {
    if (!smiles) return null;
    
    setIsLoadingDetails(true);
    
    try {
      // Use the descriptors endpoint to fetch molecule properties
      const detailsUrl = `${baseUrl}/latest/chem/descriptors`;
      
      // Set up the request with proper params
      const params = new URLSearchParams();
      params.append('smiles', smiles);
      params.append('toolkit', 'rdkit'); // Use RDKit
      params.append('format', 'json'); // Request JSON format
      
      const response = await fetch(`${detailsUrl}?${params.toString()}`);
      
      if (!response.ok) {
        throw new Error(`Failed to fetch molecule details: ${response.statusText}`);
      }
      
      const data = await response.json();
      
      // Process the descriptors based on the provided caffeine example JSON
      const processedData = {
        smiles,
        physicalProperties: {
          molecularWeight: data.molecular_weight,
          exactMolecularWeight: data.exact_molecular_weight,
          alogp: data.alogp,
          rotatableBondCount: data.rotatable_bond_count,
          polarSurfaceArea: data.topological_polar_surface_area, 
          formalCharge: data.formal_charge,
          fractionCsp3: data.fractioncsp3,
          vanDerWaalsVolume: data.van_der_waals_volume
        },
        atomCounts: {
          atomCount: data.atom_count,
          heavyAtomCount: data.heavy_atom_count,
          aromaticRingsCount: data.aromatic_rings_count,
          minimumRingCount: data.number_of_minimal_rings
        },
        drugLikeProperties: {
          hydrogenBondAcceptors: data.hydrogen_bond_acceptors,
          hydrogenBondDonors: data.hydrogen_bond_donors,
          hydrogenBondAcceptorsLipinski: data.hydrogen_bond_acceptors_lipinski,
          hydrogenBondDonorsLipinski: data.hydrogen_bond_donors_lipinski,
          lipinskiViolations: data.lipinski_rule_of_five_violations,
          qedDrugLikeliness: data.qed_drug_likeliness,
          npLikeness: data.nplikeness,
          hasLinearSugars: data.linear_sugars,
          hasCircularSugars: data.circular_sugars,
          murckoFramework: data.murcko_framework
        },
        molecularFormula: data.molecular_formula
      };
      
      setMoleculeDetails(processedData);
      return processedData;
    } catch (error) {
      console.error('Error fetching molecule descriptors:', error);
      
      // Create fallback details with SMILES parsing
      const fallbackDetails = {
        smiles,
        estimatedProperties: {
          heavyAtomCount: (smiles.match(/[A-Z][a-z]?/g) || []).length,
          carbonCount: (smiles.match(/C/g) || []).length,
          oxygenCount: (smiles.match(/O/g) || []).length,
          nitrogenCount: (smiles.match(/N/g) || []).length,
          ringCount: (smiles.match(/[0-9]/g) || []).length / 2
        },
        error: error.message
      };
      
      setMoleculeDetails(fallbackDetails);
      return fallbackDetails;
    } finally {
      setIsLoadingDetails(false);
    }
  };

  // Handle showing molecule details
  const handleShowDetails = async (e) => {
    e.stopPropagation(); // Prevent card click
    
    if (!smiles) {
      console.warn('No SMILES available to show details');
      return;
    }
    
    setShowDetailsModal(true);
    
    // Fetch details if we don't have them yet
    if (!moleculeDetails) {
      await fetchMoleculeDetails(smiles);
    }
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

  // Render property sections in the details modal with formatted values
  const renderPropertySection = (sectionTitle, properties) => {
    if (!properties || Object.keys(properties).length === 0) return null;
    
    // Format values based on property type
    const formatValue = (key, value) => {
      if (value === null || value === undefined) return 'N/A';
      
      // Format numbers with appropriate precision
      if (typeof value === 'number') {
        // Different precisions for different property types
        if (key.toLowerCase().includes('weight') || key.toLowerCase().includes('mass')) {
          return value.toFixed(2) + ' g/mol';
        } else if (key.toLowerCase().includes('logp')) {
          return value.toFixed(2);
        } else if (key.toLowerCase().includes('surface') || key.toLowerCase().includes('area')) {
          return value.toFixed(2) + ' Å²';
        } else if (Number.isInteger(value)) {
          return value.toString();
        } else {
          return value.toFixed(2);
        }
      }
      
      return value.toString();
    };
    
    // Format property name for display
    const formatPropertyName = (key) => {
      // Special case handling for common abbreviations
      const specialCases = {
        'logP': 'LogP',
        'tpsa': 'TPSA',
        'hBondDonors': 'H-Bond Donors',
        'hBondAcceptors': 'H-Bond Acceptors',
        'molecularWeight': 'Molecular Weight',
        'polarSurfaceArea': 'Polar Surface Area',
        'rotatableBonds': 'Rotatable Bonds',
        'heavyAtomCount': 'Heavy Atom Count',
        'aromaticRingCount': 'Aromatic Ring Count'
      };
      
      if (specialCases[key]) return specialCases[key];
      
      // Default formatting: camelCase/PascalCase to Title Case with spaces
      return key
        .replace(/([A-Z])/g, ' $1') // Insert space before capital letters
        .replace(/^./, (str) => str.toUpperCase()) // Capitalize first letter
        .trim();
    };
    
    return (
      <div className="mb-4">
        <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">{sectionTitle}</h4>
        <div className="bg-gray-50 dark:bg-gray-900 rounded-md p-3">
          {Object.entries(properties)
            .filter(([_, value]) => value !== null) // Filter out null values
            .map(([key, value]) => (
              <div key={key} className="flex justify-between py-1 text-sm border-b border-gray-100 dark:border-gray-800 last:border-0">
                <span className="text-gray-600 dark:text-gray-400">{formatPropertyName(key)}</span>
                <span className="font-mono text-gray-800 dark:text-gray-200">{formatValue(key, value)}</span>
              </div>
          ))}
        </div>
      </div>
    );
  };

  return (
    <>
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

            {/* Info Button (Implemented) */}
            <button
              onClick={handleShowDetails}
              disabled={!smiles} // Disable if no SMILES
              className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 ${!smiles ? 'text-gray-300 dark:text-gray-600 cursor-not-allowed' :
                  'text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700'
                }`}
              title="View molecule details"
              aria-label="View molecule details"
            >
              <HiOutlineInformationCircle className="h-5 w-5" />
            </button>
          </div>
        )}
      </div>

      {/* Manual Copy Modal (Shown when automatic copy fails) */}
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

      {/* Molecule Details Modal */}
      {showDetailsModal && (
        <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4">
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-2xl max-w-2xl w-full max-h-[90vh] overflow-y-auto">
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-bold text-gray-900 dark:text-white">
                {title !== 'Molecule' ? title : 'Molecule Details'}
              </h3>
              <button 
                onClick={() => setShowDetailsModal(false)}
                className="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200"
              >
                <HiOutlineX className="w-5 h-5" />
              </button>
            </div>

            {/* Loading state */}
            {isLoadingDetails && (
              <div className="py-8 text-center">
                <div className="inline-block h-8 w-8 animate-spin rounded-full border-4 border-solid border-blue-500 border-r-transparent"></div>
                <p className="mt-4 text-gray-600 dark:text-gray-400">Loading molecule information...</p>
              </div>
            )}

            {/* Error state */}
            {!isLoadingDetails && moleculeDetails?.error && (
              <div className="py-4">
                <div className="bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-md p-4 mb-4">
                  <p className="text-red-700 dark:text-red-400">
                    Could not fetch complete molecule details. Showing limited information.
                  </p>
                </div>
              </div>
            )}

            {/* Content when details are loaded */}
                          {!isLoadingDetails && moleculeDetails && (
              <div className="space-y-4">
                {/* Structure visualization */}
                <div className="bg-gray-50 dark:bg-gray-900 rounded-lg p-4">
                  <div className="flex justify-between items-center mb-2">
                    <span className="text-xs text-gray-500 dark:text-gray-400">
                      Structure visualization
                    </span>
                  </div>
                  <div className="flex justify-center">
                    <img 
                      src={imageUrl}
                      alt={`Chemical structure of ${title}`}
                      className="max-h-64 max-w-full object-contain"
                      onError={(e) => {
                        e.target.onerror = null;
                        e.target.src = FALLBACK_SVG_BASE64;
                      }}
                    />
                  </div>
                </div>

                {/* Basic Information */}
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  {/* SMILES */}
                  <div className="bg-gray-50 dark:bg-gray-900 rounded-md p-3">
                    <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">SMILES</h4>
                    <div className="font-mono text-xs break-all text-gray-800 dark:text-gray-200 p-2 bg-white dark:bg-gray-800 rounded border border-gray-200 dark:border-gray-700">
                      {smiles || 'N/A'}
                    </div>
                  </div>

                  {/* Identification */}
                  {moleculeDetails.inchi && (
                    <div className="bg-gray-50 dark:bg-gray-900 rounded-md p-3">
                      <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">InChI</h4>
                      <div className="font-mono text-xs overflow-x-auto whitespace-nowrap text-gray-800 dark:text-gray-200 p-2 bg-white dark:bg-gray-800 rounded border border-gray-200 dark:border-gray-700">
                        {moleculeDetails.inchi || 'N/A'}
                      </div>
                    </div>
                  )}
                </div>

                {/* Physical Properties */}
                {renderPropertySection('Physical Properties', moleculeDetails.physicalProperties || {})}

                {/* Atom Counts */}
                {renderPropertySection('Atom Count', moleculeDetails.atomCounts || {})}
                
                {/* Drug-likeness Properties */}
                {moleculeDetails.drugLikeProperties && renderPropertySection('Drug-likeness Properties', {
                  // Only include key drug-likeness properties
                  hydrogenBondAcceptors: moleculeDetails.drugLikeProperties.hydrogenBondAcceptors,
                  hydrogenBondDonors: moleculeDetails.drugLikeProperties.hydrogenBondDonors,
                  lipinskiViolations: moleculeDetails.drugLikeProperties.lipinskiViolations,
                  qedDrugLikeliness: moleculeDetails.drugLikeProperties.qedDrugLikeliness,
                  npLikeness: moleculeDetails.drugLikeProperties.npLikeness
                })}
                
                {/* QED Score Visualization */}
                {moleculeDetails.drugLikeProperties?.qedDrugLikeliness !== undefined && (
                  <div className="mb-4">
                    <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">QED Score</h4>
                    <div className="bg-gray-50 dark:bg-gray-900 rounded-md p-3">
                      <div className="mb-2 text-sm text-gray-600 dark:text-gray-400">
                        Quantitative Estimate of Drug-likeness (0-1, higher is better)
                      </div>
                      {/* QED Score progress bar */}
                      <div className="w-full bg-gray-200 dark:bg-gray-700 rounded-full h-4">
                        <div 
                          className="bg-blue-600 dark:bg-blue-500 h-4 rounded-full"
                          style={{ 
                            width: `${(moleculeDetails.drugLikeProperties.qedDrugLikeliness * 100).toFixed(0)}%`,
                            background: `linear-gradient(90deg, 
                              ${moleculeDetails.drugLikeProperties.qedDrugLikeliness < 0.3 ? '#EF4444' : 
                                moleculeDetails.drugLikeProperties.qedDrugLikeliness < 0.6 ? '#F59E0B' : 
                                '#10B981'
                              } 0%, 
                              ${moleculeDetails.drugLikeProperties.qedDrugLikeliness < 0.3 ? '#F87171' : 
                                moleculeDetails.drugLikeProperties.qedDrugLikeliness < 0.6 ? '#FBBF24' : 
                                '#34D399'
                              } 100%)`
                          }}
                        />
                      </div>
                      <div className="flex justify-between text-xs text-gray-600 dark:text-gray-400 mt-1">
                        <span>0.0</span>
                        <span className="font-medium text-blue-600 dark:text-blue-400">{moleculeDetails.drugLikeProperties.qedDrugLikeliness.toFixed(2)}</span>
                        <span>1.0</span>
                      </div>
                    </div>
                  </div>
                )}
                
                {/* Identifiers from RDKit */}
                {(moleculeDetails.inchi || moleculeDetails.inchiKey) && (
                  <div className="mb-4">
                    <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">Chemical Identifiers</h4>
                    <div className="bg-gray-50 dark:bg-gray-900 rounded-md p-3 space-y-3">
                      {/* InChI */}
                      {moleculeDetails.inchi && (
                        <div className="overflow-hidden text-sm">
                          <div className="font-medium text-gray-600 dark:text-gray-400 mb-1">InChI</div>
                          <div className="font-mono text-xs break-all bg-white dark:bg-gray-800 p-2 rounded border border-gray-200 dark:border-gray-700">
                            {moleculeDetails.inchi}
                          </div>
                        </div>
                      )}
                      
                      {/* InChI Key */}
                      {moleculeDetails.inchiKey && (
                        <div className="overflow-hidden text-sm">
                          <div className="font-medium text-gray-600 dark:text-gray-400 mb-1">InChI Key</div>
                          <div className="font-mono text-xs bg-white dark:bg-gray-800 p-2 rounded border border-gray-200 dark:border-gray-700">
                            {moleculeDetails.inchiKey}
                          </div>
                        </div>
                      )}
                    </div>
                  </div>
                )}

                {/* Molecular Formula */}
                {moleculeDetails.molecularFormula && (
                  <div className="mb-4">
                    <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">Molecular Formula</h4>
                    <div className="bg-gray-50 dark:bg-gray-900 rounded-md p-3">
                      <div className="font-mono text-lg text-center text-gray-800 dark:text-gray-200 p-2">
                        {moleculeDetails.molecularFormula}
                      </div>
                    </div>
                  </div>
                )}

                {/* Fallback information if API call failed */}
                {moleculeDetails.error && (
                  <div className="bg-gray-50 dark:bg-gray-900 rounded-md p-3">
                    <div className="flex items-center mb-2">
                      <HiOutlineExclamationCircle className="h-5 w-5 text-orange-500 dark:text-orange-400 mr-2" />
                      <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300">Estimated Information</h4>
                    </div>
                    <div className="space-y-2">
                      <p className="text-sm text-gray-600 dark:text-gray-400">
                        Unable to fetch complete molecule data. These values are rough estimates parsed from the SMILES string:
                      </p>
                      <ul className="list-disc pl-5 text-sm text-gray-700 dark:text-gray-300">
                        {Object.entries(moleculeDetails.estimatedProperties || {}).map(([key, value]) => (
                          <li key={key}>
                            <span className="font-medium">{key.replace(/([A-Z])/g, ' $1').trim()}:</span> {value}
                          </li>
                        ))}
                      </ul>
                    </div>
                    <div className="mt-3 pt-3 border-t border-gray-200 dark:border-gray-700">
                      <button
                        onClick={handleShowDetails} // Re-attempt fetch
                        className="text-sm text-blue-600 dark:text-blue-400 hover:text-blue-800 dark:hover:text-blue-300 flex items-center"
                      >
                        <svg className="h-4 w-4 mr-1" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
                        </svg>
                        Retry fetching details
                      </button>
                    </div>
                  </div>
                )}
              </div>
            )}

            {/* Action buttons */}
            <div className="mt-6 flex justify-end gap-3">
              {/* Download button in modal */}
              <button
                onClick={handleDownload}
                className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2"
              >
                Download Structure
              </button>
              
              <button
                onClick={() => setShowDetailsModal(false)}
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

export default MoleculeCard;