// Description: 3D visualization component using $3Dmol.js
// This component allows users to input SMILES strings, generate 3D structures, and visualize them.
import React, { useState, useRef, useEffect, useCallback } from "react";
// Ensure all used icons are imported
import {
  HiOutlineCube,
  HiOutlineCamera,
  HiOutlineViewGrid,
  HiOutlineRefresh,
  HiOutlineColorSwatch,
  HiOutlineClipboard,
  HiOutlineDownload,
  HiOutlineCheck,
  HiOutlineExclamationCircle,
} from "react-icons/hi";
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
// Assuming this service is configured correctly
import convertService from "../../services/convertService"; // Assuming this service exists

// Visualization styles configuration
const VISUALIZATION_STYLES = [
  { id: "stick", label: "Stick" },
  { id: "line", label: "Line" },
  { id: "sphere", label: "Sphere" },
];

// Background colors configuration (hex values used by $3Dmol)
const BACKGROUND_COLORS = [
  { id: "black", label: "Black", value: "#000000" },
  { id: "white", label: "White", value: "#ffffff" },
  { id: "gray", label: "Dark Gray", value: "#2d3748" }, // Adjusted label for clarity
  { id: "blue", label: "Dark Blue", value: "#1a202c" },
];

// Color schemes configuration (keys match $3Dmol options)
const COLOR_SCHEMES = [
  { id: "default", label: "Element (Default)" }, // Maps to undefined in $3Dmol style
  { id: "cyanCarbon", label: "Cyan Carbon" },
  { id: "greenCarbon", label: "Green Carbon" },
  { id: "magentaCarbon", label: "Magenta Carbon" },
  { id: "yellowCarbon", label: "Yellow Carbon" },
  { id: "whiteCarbon", label: "White Carbon" },
  { id: "spectrum", label: "Rainbow Spectrum" },
];

// Toolkit options configuration for 3D generation
const TOOLKIT_OPTIONS_3D = [
  { id: "openbabel", label: "Open Babel" },
  { id: "rdkit", label: "RDKit" },
  // Add other toolkits your service supports for 3D generation
];

// Example molecule SMILES
const EXAMPLE_MOLECULE = "CC1(C)OC2COC3(COS(N)(=O)=O)OC(C)(C)OC3C2O1"; // Example: (+)-Camphorsultam

const Depict3DView = ({ isActive = true }) => {
  // isActive prop might control if viewer should mount/render
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false); // For API calls
  const [error, setError] = useState(null);
  const [toolkit, setToolkit] = useState("openbabel"); // Default toolkit for 3D
  const [style, setStyle] = useState("stick");
  const [backgroundColor, setBackgroundColor] = useState("white"); // Default bg color ID
  const [colorScheme, setColorScheme] = useState("default");
  const [showLabels, setShowLabels] = useState(false);
  const [spin, setSpin] = useState(false);
  const [molData, setMolData] = useState(null); // Stores the 3D molblock string
  const [viewerInitialized, setViewerInitialized] = useState(false); // Tracks if $3Dmol viewer is ready
  const [viewerMounted, setViewerMounted] = useState(false); // Tracks if component/container is mounted
  const [copied, setCopied] = useState(false); // State for copy button feedback

  // Refs for the viewer instance and its container
  const viewerRef = useRef(null); // Stores the $3Dmol viewer instance
  const containerRef = useRef(null); // Ref to the div where the viewer will be mounted
  const modelRef = useRef(null); // Stores the $3Dmol model instance
  const mountTimerRef = useRef(null); // Timer ref for delayed mounting

  // Generate 3D structure from SMILES using the selected toolkit
  const generateStructure = useCallback(
    async (smilesStr) => {
      if (!smilesStr) return;
      if (
        !convertService ||
        typeof convertService.generate3DCoordinates !== "function"
      ) {
        console.error("convertService.generate3DCoordinates is not available.");
        setError("3D Conversion service is not configured correctly.");
        return;
      }
      setLoading(true);
      setError(null);
      setMolData(null);
      try {
        const molblockResult = await convertService.generate3DCoordinates(
          smilesStr,
          toolkit
        );
        if (
          typeof molblockResult === "string" &&
          molblockResult.trim() !== ""
        ) {
          setMolData(molblockResult);
        } else {
          throw new Error("Invalid or empty molecule data received from API.");
        }
      } catch (err) {
        console.error("Error generating 3D structure:", err);
        setError(
          `Error generating 3D structure: ${err.message || "Unknown error"}`
        );
        setMolData(null);
      } finally {
        setLoading(false);
      }
    },
    [toolkit]
  );

  // Handle form submission
  const handleSubmit = async (e) => {
    if (e) e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (trimmedSmiles) {
      await generateStructure(trimmedSmiles);
    } else {
      setError("Please enter a SMILES string.");
    }
  };

  // Destroy and cleanup the $3Dmol viewer instance
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
        setViewerInitialized(false);
      }
    }
  }, []);

  // Create a new $3Dmol viewer instance
  const createViewer = useCallback(() => {
    if (!containerRef.current || !window.$3Dmol) return null;
    const containerWidth = containerRef.current.clientWidth;
    const containerHeight = containerRef.current.clientHeight;
    if (containerWidth <= 0 || containerHeight <= 0) return null;
    try {
      const bgColorHex =
        BACKGROUND_COLORS.find((c) => c.id === backgroundColor)?.value ||
        "#000000";
      destroyViewer();
      const viewer = window.$3Dmol.createViewer(containerRef.current, {
        backgroundColor: bgColorHex,
        width: "100%",
        height: "100%",
        antialias: true,
        defaultcolors: window.$3Dmol.elementColors.rasmol,
      });
      if (!viewer)
        throw new Error("$3Dmol.createViewer returned null or undefined.");
      viewer.resize();
      return viewer;
    } catch (err) {
      console.error("Failed to create $3Dmol viewer:", err);
      setError(`Failed to initialize 3D viewer: ${err.message}`);
      return null;
    }
  }, [backgroundColor, destroyViewer]);

  // --- Apply visual style settings (LABEL LOGIC MODIFIED HERE) ---
  const applyStyle = useCallback(() => {
    // Ensure viewer and model references are valid
    if (!viewerRef.current || !modelRef.current) return false;
    const viewer = viewerRef.current;

    try {
      // Reset existing styles and labels before applying new ones
      viewer.setStyle({}, {}); // Clear previous styles for the model
      viewer.removeAllLabels(); // Remove all existing labels

      // Determine the color scheme for the molecule visualization
      // Use undefined for the default element-based coloring in $3Dmol
      const styleColorScheme =
        colorScheme === "default" ? undefined : colorScheme;
      // Base configuration for styles, including the selected color scheme
      const styleConfig = { colorscheme: styleColorScheme };

      // Apply the selected visualization style (stick, line, sphere)
      switch (style) {
        case "stick":
          viewer.setStyle({}, { stick: { radius: 0.15, ...styleConfig } });
          break;
        case "line":
          viewer.setStyle({}, { line: { linewidth: 2, ...styleConfig } });
          break;
        case "sphere":
          viewer.setStyle({}, { sphere: { scale: 0.3, ...styleConfig } });
          break;
        default: // Default to stick style if selection is somehow invalid
          viewer.setStyle({}, { stick: { ...styleConfig } });
      }

      // --- MODIFIED LABEL LOGIC ---
      // Add atom labels if the 'showLabels' option is enabled and a model exists
      if (showLabels && modelRef.current) {
        try {
          const model = modelRef.current;
          // Ensure the model has atoms data
          const atoms = model.atoms || [];

          // Detect if dark mode is active to choose the appropriate label color
          // This assumes a standard 'dark' class on the <html> element
          const isDarkMode =
            typeof window !== "undefined" &&
            document.documentElement.classList.contains("dark");

          // Set label font color: white for dark mode, a very dark gray for light mode for better contrast
          // Determine label font color based on background color for optimal contrast
          let labelFontColor;
          switch (backgroundColor) {
            case "white":
              labelFontColor = "#000000"; // Dark text on white background
              break;
            case "black":
            case "gray":
            case "blue":
              labelFontColor = "white"; // White text on dark backgrounds
              break;
            default:
              // Fallback to dark mode detection
              labelFontColor = isDarkMode ? "white" : "#1f2937";
          }

          // Define a small offset to position labels slightly away from the atom center
          const labelOffset = { x: 0.1, y: 0.1, z: 0 };

          // Iterate through each atom in the model
          for (let i = 0; i < atoms.length; i++) {
            const atom = atoms[i];
            // Add label only for non-Hydrogen atoms with valid coordinates
            if (
              atom &&
              atom.elem &&
              atom.elem !== "H" &&
              typeof atom.x === "number" &&
              typeof atom.y === "number" &&
              typeof atom.z === "number"
            ) {
              viewer.addLabel(atom.elem, {
                // Add the element symbol as the label text
                position: {
                  // Position the label near the atom with the offset
                  x: atom.x + labelOffset.x,
                  y: atom.y + labelOffset.y,
                  z: atom.z + labelOffset.z,
                },
                useScreen: false, // Position label in 3D space (relative to atom)
                fontColor: labelFontColor, // Apply the theme-aware font color
                fontSize: 14, // Set the font size
                showBackground: false,
                borderThickness: 1.0,
                borderColor: "cyan", // Set border color to cyan
                alignment: "center", // Center the label text horizontally
              });
            }
          }
        } catch (labelError) {
          // Log a warning if adding labels fails, but don't block rendering
          console.warn("Could not add atom labels:", labelError);
        }
      }
      // --- END MODIFIED LABEL LOGIC ---

      // Apply spin animation if enabled
      viewer.spin(spin);
      // Render the changes in the viewer
      viewer.render();
      return true; // Indicate success
    } catch (e) {
      // Catch and log errors during style application
      console.error("Error applying style:", e);
      setError(`Failed to apply style: ${e.message}`);
      return false; // Indicate failure
    }
  }, [style, colorScheme, showLabels, spin, backgroundColor]); // Dependencies for the useCallback hook

  // Render the molecule data in the viewer
  const renderMolecule = useCallback(() => {
    if (!viewerRef.current || !molData) return;
    const viewer = viewerRef.current;
    try {
      viewer.removeAllModels();
      modelRef.current = null;
      if (typeof molData !== "string" || !molData.includes("M  END"))
        throw new Error("Invalid Molblock data format.");
      const model = viewer.addModel(molData, "mol");
      if (!model) throw new Error("$3Dmol.addModel failed to return a model.");
      modelRef.current = model;
      const styleApplied = applyStyle();
      if (styleApplied) {
        viewer.zoomTo();
        viewer.render();
      } else {
        viewer.render();
      }
    } catch (e) {
      console.error("Error rendering molecule:", e);
      setError(
        `Failed to render molecule: ${e.message}. Try reinitializing or refreshing.`
      );
      destroyViewer();
    }
  }, [molData, applyStyle, destroyViewer]);

  // Initialize the viewer instance
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

  // Effect to mount/initialize the viewer
  useEffect(() => {
    if (isActive && !viewerMounted && containerRef.current) {
      mountTimerRef.current = setTimeout(() => {
        setViewerMounted(true);
        initializeViewer();
      }, 150);
    }
    return () => {
      if (mountTimerRef.current) clearTimeout(mountTimerRef.current);
    };
  }, [isActive, viewerMounted, initializeViewer]);

  // Effect to re-render when molData changes
  useEffect(() => {
    if (molData && viewerInitialized) renderMolecule();
    else if (!molData && viewerRef.current) {
      viewerRef.current.removeAllModels();
      modelRef.current = null;
      viewerRef.current.render();
    }
  }, [molData, viewerInitialized, renderMolecule]);

  // Effect to apply style changes
  useEffect(() => {
    if (viewerInitialized && viewerRef.current && modelRef.current)
      applyStyle();
  }, [style, colorScheme, showLabels, spin, viewerInitialized, applyStyle]);

  // Effect for background color changes
  useEffect(() => {
    if (viewerRef.current && viewerInitialized) {
      try {
        const bgColorHex =
          BACKGROUND_COLORS.find((c) => c.id === backgroundColor)?.value ||
          "#000000";
        viewerRef.current.setBackgroundColor(bgColorHex);
        viewerRef.current.render();
      } catch (e) {
        console.warn("Error changing background color:", e);
      }
    }
  }, [backgroundColor, viewerInitialized]);

  // Effect for component unmount cleanup
  useEffect(() => {
    return () => {
      destroyViewer();
    };
  }, [destroyViewer]);

  // Effect for window resize
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
    window.addEventListener("resize", handleResize);
    return () => window.removeEventListener("resize", handleResize);
  }, [viewerInitialized]);

  // Handler to use the example molecule
  const handleUseExampleMolecule = () => {
    setSmiles(EXAMPLE_MOLECULE);
    generateStructure(EXAMPLE_MOLECULE);
  };

  // Handler to take a screenshot
  const handleTakeScreenshot = () => {
    if (!viewerRef.current || !viewerInitialized) {
      setError("Viewer not ready.");
      return;
    }
    try {
      const pngDataUri = viewerRef.current.pngURI();
      if (!pngDataUri) throw new Error("Failed to generate PNG data URI.");
      const a = document.createElement("a");
      a.href = pngDataUri;
      a.download = "molecule-3d-screenshot.png";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    } catch (err) {
      console.error("Error taking screenshot:", err);
      setError("Failed to take screenshot.");
    }
  };

  // Handler to reset the camera view
  const handleResetView = () => {
    if (!viewerRef.current || !viewerInitialized) return;
    try {
      viewerRef.current.zoomTo();
      viewerRef.current.render();
    } catch (e) {
      console.warn("Error resetting view:", e);
    }
  };

  // Handler to attempt reinitialization if viewer failed
  const handleReinitialize = () => {
    setError(null);
    setViewerMounted(false);
    destroyViewer();
    setTimeout(() => {
      setViewerMounted(true);
    }, 100);
  };

  // Handle copying the molblock to clipboard
  const handleCopy = () => {
    if (!molData || !navigator.clipboard) return;
    navigator.clipboard
      .writeText(molData)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      })
      .catch((err) => {
        console.error("Failed to copy 3D molblock:", err);
        setError("Failed to copy 3D Molblock.");
      });
  };

  // Handle downloading the molblock as a .mol file
  const downloadMolblock = () => {
    if (!molData) return;
    try {
      const blob = new Blob([molData], {
        type: "chemical/x-mdl-molfile;charset=utf-8",
      });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      const filenameBase = smiles
        ? smiles.replace(/[^a-z0-9]/gi, "_").substring(0, 30)
        : "molecule";
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
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input and Options Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        {/* Title */}
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          3D Coordinate Generation & Visualization
        </h2>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          {/* SMILES Input */}
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            label="Input SMILES"
            required
          />

          {/* Options Grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 pt-2">
            {/* Toolkit Selection */}
            <div>
              <label
                htmlFor="toolkit-select-3d"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Toolkit (for 3D generation)
              </label>
              <select
                id="toolkit-select-3d"
                value={toolkit}
                onChange={(e) => setToolkit(e.target.value)}
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                {TOOLKIT_OPTIONS_3D.map((option) => (
                  <option key={option.id} value={option.id}>
                    {option.label}
                  </option>
                ))}
              </select>
            </div>
            {/* Visualization Style */}
            <div>
              <label
                htmlFor="style-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Visualization Style
              </label>
              <select
                id="style-select"
                value={style}
                onChange={(e) => setStyle(e.target.value)}
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                {VISUALIZATION_STYLES.map((visualStyle) => (
                  <option key={visualStyle.id} value={visualStyle.id}>
                    {visualStyle.label}
                  </option>
                ))}
              </select>
            </div>
            {/* Background Color */}
            <div>
              <label
                htmlFor="bg-color-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Background Color
              </label>
              <select
                id="bg-color-select"
                value={backgroundColor}
                onChange={(e) => setBackgroundColor(e.target.value)}
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                {BACKGROUND_COLORS.map((color) => (
                  <option key={color.id} value={color.id}>
                    {color.label}
                  </option>
                ))}
              </select>
            </div>
            {/* Color Scheme */}
            <div>
              <label
                htmlFor="color-scheme-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Color Scheme
              </label>
              <select
                id="color-scheme-select"
                value={colorScheme}
                onChange={(e) => setColorScheme(e.target.value)}
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                {COLOR_SCHEMES.map((scheme) => (
                  <option key={scheme.id} value={scheme.id}>
                    {scheme.label}
                  </option>
                ))}
              </select>
            </div>
          </div>

          {/* Checkbox Options */}
          <div className="flex flex-wrap gap-x-6 gap-y-2 pt-2">
            <div className="flex items-center">
              <input
                id="show-labels"
                type="checkbox"
                checked={showLabels}
                onChange={(e) => setShowLabels(e.target.checked)}
                className="h-4 w-4 rounded border-gray-300 dark:border-gray-600 text-blue-600 dark:text-blue-500 shadow-sm focus:ring-indigo-500 dark:focus:ring-blue-500 dark:focus:ring-offset-gray-800 bg-white dark:bg-gray-700"
              />
              <label
                htmlFor="show-labels"
                className="ml-2 text-sm text-gray-700 dark:text-gray-300"
              >
                Show atom labels (non-H)
              </label>
            </div>
          </div>

          {/* Action Buttons */}
          <div className="flex flex-wrap gap-3 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Generate Button */}
            <button
              type="submit"
              disabled={!smiles.trim() || loading}
              className={`px-5 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !smiles.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm"
              }`}
            >
              <HiOutlineCube className="mr-2 h-5 w-5" />
              {loading ? "Generating..." : "Generate 3D View"}
            </button>
            {/* Screenshot Button */}
            <button
              type="button"
              onClick={handleTakeScreenshot}
              disabled={loading || !molData || !viewerInitialized}
              className={`px-4 py-2 rounded-lg font-medium flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 ${
                loading || !molData || !viewerInitialized
                  ? "bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed"
                  : "bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600 shadow-sm"
              }`}
            >
              <HiOutlineCamera className="mr-2 h-5 w-5" />
              Screenshot
            </button>
            {/* Spin Toggle Button */}
            <button
              type="button"
              onClick={() => setSpin(!spin)}
              disabled={loading || !molData || !viewerInitialized}
              className={`px-4 py-2 rounded-lg font-medium flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 ${
                loading || !molData || !viewerInitialized
                  ? "bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed"
                  : "bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600 shadow-sm"
              }`}
            >
              <HiOutlineRefresh className="mr-2 h-5 w-5" />
              {spin ? "Stop Rotation" : "Rotate"}
            </button>
            {/* Reset View Button */}
            <button
              type="button"
              onClick={handleResetView}
              disabled={loading || !molData || !viewerInitialized}
              className={`px-4 py-2 rounded-lg font-medium flex items-center transition-colors duration-150 focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 ${
                loading || !molData || !viewerInitialized
                  ? "bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed"
                  : "bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-700 dark:text-gray-200 border border-gray-300 dark:border-gray-600 shadow-sm"
              }`}
            >
              <HiOutlineColorSwatch className="mr-2 h-5 w-5" />{" "}
              {/* Using ColorSwatch as placeholder for reset */}
              Reset View
            </button>
          </div>

          {/* Use Example Button */}
          {!smiles &&
            !molData && ( // Show only if input is empty and no data loaded
              <div className="mt-3 text-center">
                <button
                  type="button"
                  onClick={handleUseExampleMolecule}
                  className="text-sm font-medium text-blue-600 hover:text-blue-700 dark:text-blue-400 dark:hover:text-blue-300 focus:outline-none focus-visible:ring-2 focus-visible:ring-blue-500 rounded"
                >
                  Use example molecule?
                </button>
              </div>
            )}

          {/* Error Display with Retry */}
          {error && !loading && (
            <div
              className="p-4 bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 rounded-lg mt-4 shadow-sm flex justify-between items-start"
              role="alert"
            >
              <div className="flex items-start">
                <HiOutlineExclamationCircle
                  className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
                  aria-hidden="true"
                />
                <span>{error}</span>
              </div>
              {/* Show retry button only if viewer initialization failed */}
              {error.includes("initialize") && (
                <button
                  type="button"
                  onClick={handleReinitialize}
                  className="ml-4 px-2 py-1 bg-red-600 hover:bg-red-700 dark:bg-red-700 dark:hover:bg-red-600 rounded text-white text-xs font-medium focus:outline-none focus:ring-2 focus:ring-offset-1 dark:focus:ring-offset-red-900 focus:ring-red-500"
                >
                  Retry Init
                </button>
              )}
            </div>
          )}
        </form>
      </div>

      {/* Loading state for API call (distinct from viewer init) */}
      {loading && (
        <div className="flex justify-center items-center p-4 bg-blue-50 dark:bg-gray-700 rounded-lg border border-blue-200 dark:border-gray-600 shadow-sm">
          <div className="animate-spin h-5 w-5 border-2 border-blue-500 dark:border-blue-400 border-t-transparent rounded-full mr-3"></div>
          <span className="text-blue-700 dark:text-blue-300">
            Generating 3D structure...
          </span>
        </div>
      )}

      {/* NEW: Side-by-side container for 3D viewer and molblock */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* 3D Viewer Container - modified height */}
        <div
          className="relative rounded-lg shadow-lg overflow-hidden border border-gray-200 dark:border-gray-700"
          style={{
            backgroundColor:
              BACKGROUND_COLORS.find((c) => c.id === backgroundColor)?.value ||
              "#000000",
          }}
        >
          {/* Container for $3Dmol.js - reduced height */}
          <div
            className="w-full h-[500px] relative"
            ref={containerRef}
            data-testid="molecule-3d-container"
          >
            {/* $3Dmol viewer mounts here */}
          </div>

          {/* Loading/Initialization Overlay - unchanged */}
          {isActive && (!viewerInitialized || !viewerMounted) && (
            <div className="absolute inset-0 flex flex-col items-center justify-center bg-gray-200/80 dark:bg-black/80 backdrop-blur-sm z-10">
              {/* Themed overlay content box */}
              <div className="text-center p-6 rounded-lg bg-white/70 dark:bg-gray-800/80 shadow-lg border border-gray-300 dark:border-gray-700">
                {/* Spinner */}
                <div className="animate-spin mb-4 h-10 w-10 border-4 border-blue-500 dark:border-blue-400 border-t-transparent rounded-full mx-auto"></div>
                <p className="text-gray-700 dark:text-white text-lg font-medium">
                  Initializing 3D viewer...
                </p>
                <p className="text-gray-500 dark:text-gray-400 text-sm mt-1">
                  Please wait a moment
                </p>
              </div>
            </div>
          )}
        </div>

        {/* Molblock Display Section - moved to be side-by-side on larger screens */}
        {molData && !loading && !error ? (
          <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
            <div className="flex justify-between items-center mb-3 border-b border-gray-200 dark:border-gray-700 pb-2">
              <h3 className="text-lg font-semibold text-gray-900 dark:text-white">
                Generated 3D Molblock
              </h3>
              <div className="flex space-x-2">
                <button
                  onClick={handleCopy}
                  className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 ${
                    copied
                      ? "text-green-500 dark:text-green-500"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700"
                  }`}
                  title={copied ? "Copied!" : "Copy Molblock"}
                  aria-label={copied ? "Molblock Copied" : "Copy Molblock"}
                >
                  {copied ? (
                    <HiOutlineCheck className="h-5 w-5" />
                  ) : (
                    <HiOutlineClipboard className="h-5 w-5" />
                  )}
                </button>
                <button
                  onClick={downloadMolblock}
                  className="p-1.5 rounded-md text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700 transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500"
                  title="Download Molblock (.mol)"
                  aria-label="Download Molblock"
                >
                  <HiOutlineDownload className="h-5 w-5" />
                </button>
              </div>
            </div>
            <div className="mt-2 p-3 bg-gray-100 dark:bg-gray-900 rounded-lg font-mono text-xs overflow-auto max-h-[500px] border border-gray-200 dark:border-gray-700 shadow-sm">
              <pre className="whitespace-pre text-gray-700 dark:text-gray-300">
                {molData}
              </pre>
            </div>
          </div>
        ) : (
          // Empty placeholder div to maintain grid layout when molblock isn't shown
          <div className="hidden lg:block"></div>
        )}
      </div>
      {/* Method and Disclaimer Section */}
      {molData && !loading && (
        <div className="bg-yellow-50 dark:bg-yellow-900 dark:bg-opacity-20 border border-yellow-200 dark:border-yellow-800 rounded-lg p-4 text-sm shadow-sm">
          <div className="flex items-start space-x-3">
            <HiOutlineExclamationCircle
              className="h-5 w-5 text-yellow-600 dark:text-yellow-400 flex-shrink-0 mt-0.5"
              aria-hidden="true"
            />
            <div>
              <h4 className="font-semibold text-yellow-800 dark:text-yellow-300">
                Important Note on 3D Structure
              </h4>
              <p className="mt-1 text-gray-700 dark:text-gray-300">
                This 3D structure is computationally generated using{" "}
                {toolkit === "rdkit"
                  ? "RDKit (ETKDG/MMFF94)"
                  : "OpenBabel (MMFF94)"}
                and may not represent the actual molecular conformation found in
                nature or experimental conditions. For accurate structural
                information, refer to experimental data from X-ray
                crystallography, NMR studies, or high-level quantum-chemical
                calculations.
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Interaction Instructions */}
      {/* Themed instructions box */}
      <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 text-sm shadow">
        <div className="flex items-start space-x-3">
          <HiOutlineViewGrid
            className="h-5 w-5 text-blue-600 dark:text-blue-400 flex-shrink-0 mt-0.5"
            aria-hidden="true"
          />
          <div>
            <h4 className="font-semibold text-blue-800 dark:text-blue-300">
              Interacting with the 3D View
            </h4>
            <ul className="mt-2 space-y-1 list-disc list-inside text-gray-700 dark:text-gray-300">
              <li>
                <strong>Rotate:</strong> Click & drag
              </li>
              <li>
                <strong>Zoom:</strong> Scroll wheel / Pinch gesture
              </li>
              <li>
                <strong>Pan:</strong> Right-click drag / Two-finger drag
              </li>
              <li>
                <strong>Reset View:</strong> Double-click / Reset View button
              </li>
            </ul>
          </div>
        </div>
      </div>
    </div>
  );
};

export default Depict3DView;
