import React, { useState, useRef, useEffect, useCallback } from 'react';
import {
  HiOutlineCamera,
  HiOutlineRefresh, // Changed from HiOutlineColorSwatch to appropriate reset icon
  HiOutlineDownload,
  HiOutlineExclamationCircle
} from 'react-icons/hi';

// Assuming this service is configured correctly
import convertService from '../../services/convertService';

// Visualization styles and configurations
const VISUALIZATION_STYLES = [
  { id: 'stick', label: 'Stick' },
  { id: 'line', label: 'Line' },
  { id: 'sphere', label: 'Sphere' }
];

const BACKGROUND_COLORS = [
  { id: 'white', label: 'White', value: '#ffffff' },
  { id: 'black', label: 'Black', value: '#000000' },
  { id: 'gray', label: 'Dark Gray', value: '#2d3748' }
];

const COLOR_SCHEMES = [
  { id: 'default', label: 'Element (Default)' },
  { id: 'cyanCarbon', label: 'Cyan Carbon' },
  { id: 'greenCarbon', label: 'Green Carbon' }
];

const MoleculeDepiction3D = ({ smiles, title, toolkit = 'openbabel' }) => {
  // State for 3D viewer
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [style, setStyle] = useState('stick');
  const [backgroundColor, setBackgroundColor] = useState('white');
  const [colorScheme, setColorScheme] = useState('default');
  const [showLabels, setShowLabels] = useState(false);
  const [spin, setSpin] = useState(false);
  const [molData, setMolData] = useState(null);
  const [viewerInitialized, setViewerInitialized] = useState(false);
  const [viewerMounted, setViewerMounted] = useState(false);

  // Refs
  const viewerRef = useRef(null);
  const containerRef = useRef(null);
  const modelRef = useRef(null);
  const mountTimerRef = useRef(null);
  const initialViewStateRef = useRef(null); // New ref to store initial view state

  // Generate 3D structure from SMILES
  const generateStructure = useCallback(async (smilesStr) => {
    if (!smilesStr) return;
    if (!convertService || typeof convertService.generate3DCoordinates !== 'function') {
      console.error("convertService.generate3DCoordinates is not available.");
      setError("3D Conversion service is not configured correctly.");
      return;
    }
    
    setLoading(true);
    setError(null);
    setMolData(null);
    
    try {
      const molblockResult = await convertService.generate3DCoordinates(smilesStr, toolkit);
      if (typeof molblockResult === 'string' && molblockResult.trim() !== '') {
        setMolData(molblockResult);
      } else {
        throw new Error('Invalid or empty molecule data received from API.');
      }
    } catch (err) {
      console.error("Error generating 3D structure:", err);
      setError(`Error generating 3D structure: ${err.message || 'Unknown error'}`);
      setMolData(null);
    } finally {
      setLoading(false);
    }
  }, [toolkit]);

  // Destroy and cleanup viewer
  const destroyViewer = useCallback(() => {
    if (viewerRef.current) {
      try {
        viewerRef.current.removeAllModels();
        viewerRef.current.removeAllLabels();
        viewerRef.current.clear();
      } catch (e) { 
        console.warn("Error cleaning up viewer:", e); 
      } finally {
        viewerRef.current = null; 
        modelRef.current = null;
        initialViewStateRef.current = null; // Clear the initial view state
        setViewerInitialized(false);
      }
    }
  }, []);

  // Create new viewer instance
  const createViewer = useCallback(() => {
    if (!containerRef.current || !window.$3Dmol) return null;
    
    const containerWidth = containerRef.current.clientWidth;
    const containerHeight = containerRef.current.clientHeight;
    
    if (containerWidth <= 0 || containerHeight <= 0) return null;
    
    try {
      const bgColorHex = BACKGROUND_COLORS.find(c => c.id === backgroundColor)?.value || '#ffffff';
      destroyViewer();
      
      const viewer = window.$3Dmol.createViewer(
        containerRef.current,
        { 
          backgroundColor: bgColorHex, 
          width: '100%', 
          height: '100%', 
          antialias: true, 
          defaultcolors: window.$3Dmol.elementColors.rasmol 
        }
      );
      
      if (!viewer) throw new Error("$3Dmol.createViewer returned null or undefined.");
      viewer.resize();
      return viewer;
    } catch (err) {
      console.error("Failed to create $3Dmol viewer:", err);
      setError(`Failed to initialize 3D viewer: ${err.message}`);
      return null;
    }
  }, [backgroundColor, destroyViewer]);

  // Apply visual style settings
  const applyStyle = useCallback(() => {
    if (!viewerRef.current || !modelRef.current) return false;
    const viewer = viewerRef.current;
    
    try {
      // Reset styles and labels
      viewer.setStyle({}, {});
      viewer.removeAllLabels();
      
      // Apply color scheme
      const styleColorScheme = colorScheme === 'default' ? undefined : colorScheme;
      const styleConfig = { colorscheme: styleColorScheme };
      
      // Apply visualization style
      switch (style) {
        case 'stick':
          viewer.setStyle({}, { stick: { radius: 0.15, ...styleConfig } });
          break;
        case 'line':
          viewer.setStyle({}, { line: { linewidth: 2, ...styleConfig } });
          break;
        case 'sphere':
          viewer.setStyle({}, { sphere: { scale: 0.3, ...styleConfig } });
          break;
        default:
          viewer.setStyle({}, { stick: { ...styleConfig } });
      }
      
      // Add atom labels if enabled
      if (showLabels && modelRef.current) {
        try {
          const model = modelRef.current;
          const atoms = model.atoms || [];
          
          // Determine label font color based on background
          let labelFontColor;
          switch (backgroundColor) {
            case 'white': labelFontColor = 'black'; break;
            case 'black':
            case 'gray': labelFontColor = 'white'; break;
            default: labelFontColor = 'black';
          }
          
          // Label offset
          const labelOffset = { x: 0.1, y: 0.1, z: 0 };
          
          // Add labels to non-hydrogen atoms
          for (let i = 0; i < atoms.length; i++) {
            const atom = atoms[i];
            if (atom && atom.elem && atom.elem !== 'H' && 
                typeof atom.x === 'number' && 
                typeof atom.y === 'number' && 
                typeof atom.z === 'number') {
              viewer.addLabel(atom.elem, {
                position: {
                  x: atom.x + labelOffset.x,
                  y: atom.y + labelOffset.y,
                  z: atom.z + labelOffset.z
                },
                useScreen: false,
                fontColor: labelFontColor,
                fontSize: 14,
                showBackground: false,
                alignment: 'center'
              });
            }
          }
        } catch (labelError) {
          console.warn("Could not add atom labels:", labelError);
        }
      }
      
      // Apply spin animation
      viewer.spin(spin);
      viewer.render();
      return true;
    } catch (e) {
      console.error("Error applying style:", e);
      setError(`Failed to apply style: ${e.message}`);
      return false;
    }
  }, [style, colorScheme, showLabels, spin, backgroundColor]);

  // Render molecule
  const renderMolecule = useCallback(() => {
    if (!viewerRef.current || !molData) return;
    const viewer = viewerRef.current;
    
    try {
      viewer.removeAllModels(); 
      modelRef.current = null;
      
      if (typeof molData !== 'string' || !molData.includes('M  END')) {
        throw new Error('Invalid Molblock data format.');
      }
      
      const model = viewer.addModel(molData, "mol");
      if (!model) throw new Error('$3Dmol.addModel failed to return a model.');
      
      modelRef.current = model;
      const styleApplied = applyStyle();
      
      if (styleApplied) {
        viewer.zoomTo();
        viewer.render();
        
        // Store the initial view state after first render
        if (!initialViewStateRef.current && viewer.getView) {
          try {
            initialViewStateRef.current = viewer.getView();
          } catch (viewErr) {
            console.warn("Could not store initial view state:", viewErr);
          }
        }
      } else {
        viewer.render();
      }
    } catch (e) {
      console.error("Error rendering molecule:", e);
      setError(`Failed to render molecule: ${e.message}. Try refreshing.`);
      destroyViewer();
    }
  }, [molData, applyStyle, destroyViewer]);

  // Initialize viewer
  const initializeViewer = useCallback(() => {
    if (viewerInitialized || !containerRef.current) return false;
    
    const viewer = createViewer();
    if (viewer) {
      viewerRef.current = viewer;
      setViewerInitialized(true);
      if (molData) setTimeout(() => renderMolecule(), 100);
      return true;
    } else {
      setViewerInitialized(false);
      return false;
    }
  }, [createViewer, molData, renderMolecule, viewerInitialized]);

  // Effect: Initialize on mount
  useEffect(() => {
    if (!viewerMounted && containerRef.current) {
      mountTimerRef.current = setTimeout(() => {
        setViewerMounted(true);
        initializeViewer();
      }, 150);
    }
    
    return () => {
      if (mountTimerRef.current) clearTimeout(mountTimerRef.current);
    };
  }, [viewerMounted, initializeViewer]);

  // Effect: Generate 3D structure when SMILES changes
  useEffect(() => {
    if (smiles) {
      generateStructure(smiles);
    }
  }, [smiles, toolkit, generateStructure]);

  // Effect: Render molecule when molData changes
  useEffect(() => {
    if (molData && viewerInitialized) renderMolecule();
    else if (!molData && viewerRef.current) {
      viewerRef.current.removeAllModels();
      modelRef.current = null;
      viewerRef.current.render();
    }
  }, [molData, viewerInitialized, renderMolecule]);

  // Effect: Apply style changes
  useEffect(() => {
    if (viewerInitialized && viewerRef.current && modelRef.current) applyStyle();
  }, [style, colorScheme, showLabels, spin, viewerInitialized, applyStyle]);

  // Effect: Change background color
  useEffect(() => {
    if (viewerRef.current && viewerInitialized) {
      try {
        const bgColorHex = BACKGROUND_COLORS.find(c => c.id === backgroundColor)?.value || '#ffffff';
        viewerRef.current.setBackgroundColor(bgColorHex);
        viewerRef.current.render();
      } catch (e) {
        console.warn("Error changing background color:", e);
      }
    }
  }, [backgroundColor, viewerInitialized]);

  // Effect: Clean up on unmount
  useEffect(() => {
    return () => { destroyViewer(); };
  }, [destroyViewer]);

  // Effect: Handle window resize
  useEffect(() => {
    const handleResize = () => {
      if (viewerRef.current && containerRef.current && viewerInitialized) {
        try {
          viewerRef.current.resize();
          viewerRef.current.render();
        } catch (e) {
          console.warn("Error handling resize:", e);
        }
      }
    };
    
    window.addEventListener('resize', handleResize);
    return () => window.removeEventListener('resize', handleResize);
  }, [viewerInitialized]);

  // Reset view handler - IMPROVED VERSION
  const handleResetView = () => {
    if (!viewerRef.current || !viewerInitialized) return;
    try {
      // First approach: Use the stored initial view if available
      if (initialViewStateRef.current) {
        try {
          viewerRef.current.setView(initialViewStateRef.current);
          viewerRef.current.render();
          return;
        } catch (viewErr) {
          console.warn("Could not restore initial view state, falling back to zoomTo:", viewErr);
        }
      }
      
      // Second approach: Use zoomTo with default parameters to reset
      viewerRef.current.zoomTo();
      
      // Third approach: If zoom doesn't work well, completely reset the view and rotation
      if (typeof viewerRef.current.resetView === 'function') {
        viewerRef.current.resetView();
      } else {
        // Manual reset if resetView is not available
        viewerRef.current.setRotationCenter({x: 0, y: 0, z: 0});
        viewerRef.current.setCenter({x: 0, y: 0, z: 0});
      }
      
      // Stop any spin animation temporarily and restart if needed
      const wasSpin = spin;
      if (wasSpin) {
        viewerRef.current.spin(false);
      }
      
      // Force a reapply of the current style
      applyStyle();
      
      // Restart spin if it was on
      if (wasSpin) {
        setTimeout(() => {
          if (viewerRef.current) viewerRef.current.spin(true);
        }, 100);
      }
      
      viewerRef.current.render();
    } catch (e) {
      console.warn("Error resetting view:", e);
      setError("Failed to reset view. Try refreshing the page.");
    }
  };

  // Take screenshot
  const handleTakeScreenshot = () => {
    if (!viewerRef.current || !viewerInitialized) {
      setError("Viewer not ready.");
      return;
    }
    
    try {
      const pngDataUri = viewerRef.current.pngURI();
      if (!pngDataUri) throw new Error("Failed to generate PNG data URI.");
      
      const a = document.createElement('a');
      a.href = pngDataUri;
      a.download = 'molecule-3d-screenshot.png';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    } catch (err) {
      console.error("Error taking screenshot:", err);
      setError("Failed to take screenshot.");
    }
  };

  // Download molblock
  const downloadMolblock = () => {
    if (!molData) return;
    
    try {
      const blob = new Blob([molData], { type: 'chemical/x-mdl-molfile;charset=utf-8' });
      const url = URL.createObjectURL(blob);
      
      const a = document.createElement('a');
      a.href = url;
      const filenameBase = title ? title.replace(/[^a-z0-9]/gi, '_').substring(0, 30) : 'molecule';
      a.download = `${filenameBase}_3d.mol`;
      
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    } catch (err) {
      console.error("Error creating download link:", err);
      setError("Could not create file.");
    }
  };

  return (
    <div className="bg-white dark:bg-gray-800 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700 overflow-hidden h-full flex flex-col">
      {/* Header */}
      <div className="px-4 py-3 border-b border-gray-200 dark:border-gray-700 flex justify-between items-center bg-gradient-to-r from-gray-50 to-gray-100 dark:from-gray-800 dark:to-gray-750">
        <h3 className="font-medium text-gray-800 dark:text-white text-sm sm:text-base truncate">
          {title || "3D Structure"}
        </h3>
        
        {/* Action buttons */}
        <div className="flex items-center space-x-1">
          <button
            onClick={handleTakeScreenshot}
            disabled={!viewerInitialized || !molData}
            className="p-1.5 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
            title="Take screenshot"
          >
            <HiOutlineCamera className="h-5 w-5" />
          </button>
          
          <button
            onClick={handleResetView}
            disabled={!viewerInitialized || !molData}
            className="p-1.5 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 disabled:opacity-50 disabled:cursor-not-allowed"
            title="Reset view"
          >
            <HiOutlineRefresh className="h-5 w-5" />
          </button>
        </div>
      </div>
      
      {/* Viewer container */}
      <div
        className="relative flex-grow"
        style={{ 
          backgroundColor: BACKGROUND_COLORS.find(c => c.id === backgroundColor)?.value || '#ffffff',
          minHeight: "350px"
        }}
      >
        <div
          className="w-full h-full absolute"
          style={{ minHeight: "350px" }}
          ref={containerRef}
          data-testid="molecule-3d-container"
        >
          {/* $3Dmol viewer mounts here */}
        </div>

        {/* Loading overlay */}
        {(!viewerInitialized || !viewerMounted) && (
          <div className="absolute inset-0 flex flex-col items-center justify-center bg-gray-200/80 dark:bg-black/80 backdrop-blur-sm z-10">
            <div className="text-center p-6 rounded-lg bg-white/70 dark:bg-gray-800/80 shadow-lg border border-gray-300 dark:border-gray-700">
              <div className="animate-spin mb-4 h-10 w-10 border-4 border-blue-500 dark:border-blue-400 border-t-transparent rounded-full mx-auto"></div>
              <p className="text-gray-700 dark:text-white text-lg font-medium">
                Initializing 3D viewer...
              </p>
            </div>
          </div>
        )}
        
        {/* Loading indicator during generation */}
        {loading && (
          <div className="absolute inset-0 flex flex-col items-center justify-center bg-gray-200/80 dark:bg-black/80 backdrop-blur-sm z-10">
            <div className="text-center p-6 rounded-lg bg-white/70 dark:bg-gray-800/80 shadow-lg border border-gray-300 dark:border-gray-700">
              <div className="animate-spin mb-4 h-10 w-10 border-4 border-blue-500 dark:border-blue-400 border-t-transparent rounded-full mx-auto"></div>
              <p className="text-gray-700 dark:text-white text-lg font-medium">
                Generating 3D structure...
              </p>
            </div>
          </div>
        )}
      </div>
      
      {/* Controls */}
      <div className="px-4 py-3 bg-gray-50 dark:bg-gray-900 border-t border-gray-200 dark:border-gray-700">
        <div className="flex flex-wrap gap-4 mb-2">
          {/* Style Select */}
          <div className="flex-1 min-w-[120px]">
            <label htmlFor="style-select" className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
              Style
            </label>
            <select
              id="style-select"
              value={style}
              onChange={(e) => setStyle(e.target.value)}
              className="w-full text-xs bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 py-1 px-2"
            >
              {VISUALIZATION_STYLES.map((visualStyle) => (
                <option key={visualStyle.id} value={visualStyle.id}>
                  {visualStyle.label}
                </option>
              ))}
            </select>
          </div>
          
          {/* Color Scheme Select */}
          <div className="flex-1 min-w-[120px]">
            <label htmlFor="color-scheme-select" className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
              Color
            </label>
            <select
              id="color-scheme-select"
              value={colorScheme}
              onChange={(e) => setColorScheme(e.target.value)}
              className="w-full text-xs bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 py-1 px-2"
            >
              {COLOR_SCHEMES.map((scheme) => (
                <option key={scheme.id} value={scheme.id}>
                  {scheme.label}
                </option>
              ))}
            </select>
          </div>
          
          {/* Background Color Select */}
          <div className="flex-1 min-w-[120px]">
            <label htmlFor="bg-color-select" className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
              Background
            </label>
            <select
              id="bg-color-select"
              value={backgroundColor}
              onChange={(e) => setBackgroundColor(e.target.value)}
              className="w-full text-xs bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 py-1 px-2"
            >
              {BACKGROUND_COLORS.map((color) => (
                <option key={color.id} value={color.id}>
                  {color.label}
                </option>
              ))}
            </select>
          </div>
        </div>
        
        {/* Checkboxes */}
        <div className="flex flex-wrap gap-x-4 gap-y-2">
          {/* Show Labels Checkbox */}
          <div className="flex items-center">
            <input
              id="show-labels"
              type="checkbox"
              checked={showLabels}
              onChange={(e) => setShowLabels(e.target.checked)}
              className="h-3 w-3 rounded border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
            />
            <label htmlFor="show-labels" className="ml-1.5 text-xs text-gray-600 dark:text-gray-400">
              Labels
            </label>
          </div>
          
          {/* Spin Checkbox */}
          <div className="flex items-center">
            <input
              id="spin-toggle"
              type="checkbox"
              checked={spin}
              onChange={(e) => setSpin(e.target.checked)}
              className="h-3 w-3 rounded border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
            />
            <label htmlFor="spin-toggle" className="ml-1.5 text-xs text-gray-600 dark:text-gray-400">
              Spin
            </label>
          </div>
          
          {/* Download Button */}
          <div className="ml-auto flex space-x-1">
            <button
              onClick={downloadMolblock}
              disabled={!molData}
              className="py-1 px-2 rounded-md text-xs bg-green-50 text-green-700 dark:bg-green-900 dark:text-green-200 hover:bg-green-100 dark:hover:bg-green-800 transition-colors focus:outline-none focus:ring-1 focus:ring-green-500 disabled:opacity-50 disabled:cursor-not-allowed border border-green-200 dark:border-green-800"
              title="Download 3D molfile"
            >
              <span className="flex items-center">
                <HiOutlineDownload className="h-3.5 w-3.5 mr-1" />
                3D Mol
              </span>
            </button>
          </div>
        </div>
      </div>
      
      {/* Error Display */}
      {error && !loading && (
        <div className="p-3 bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border-t border-red-300 dark:border-red-700 text-xs">
          <div className="flex items-start">
            <HiOutlineExclamationCircle className="h-4 w-4 mr-2 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
            <span>{error}</span>
          </div>
        </div>
      )}
      
      {/* Instructions */}
      <div className="px-4 py-2 bg-gray-100 dark:bg-gray-700 border-t border-gray-200 dark:border-gray-600 text-xs text-gray-500 dark:text-gray-400">
        <div className="flex items-center">
          <span className="mr-2">Rotate: click & drag</span>
          <span className="mr-2">•</span>
          <span className="mr-2">Zoom: scroll</span>
          <span className="mr-2">•</span>
          <span>Pan: right-click drag</span>
        </div>
      </div>
    </div>
  );
};

export default MoleculeDepiction3D;