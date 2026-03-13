// Description: This component allows users to generate 2D depictions of multiple molecules from SMILES strings.
import React, { useState } from "react"; // Removed useEffect
import { motion, AnimatePresence, LayoutGroup } from "motion/react";
// Assuming these components are correctly implemented and styled for dark/light mode
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { EmptyState } from "@/components/feedback/EmptyState";
import { getErrorMessage } from "@/lib/error-messages";
// Assuming this service is configured correctly
import depictService from "../../services/depictService"; // Assuming this service exists
import { AlertCircle, ArrowLeftRight, Check, Clipboard, Download, Image } from "lucide-react";
import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { Input } from "@/components/ui/input";
import { Textarea } from "@/components/ui/textarea";

// Animation variants
const resultsContainerVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { staggerChildren: 0.05 } }, // Stagger children animation
};

const depictionCardVariant = {
  hidden: { opacity: 0, y: 20, scale: 0.95 },
  visible: {
    opacity: 1,
    y: 0,
    scale: 1,
    transition: {
      duration: 0.5,
      ease: [0.25, 0.46, 0.45, 0.94],
      scale: {
        type: "spring",
        stiffness: 300,
        damping: 25,
      },
    },
  },
};

const toggleSpring = {
  type: "spring",
  stiffness: 600,
  damping: 30,
};

// Toggle Switch Component - Enhanced with stunning visuals
const ToggleSwitch = ({ id, checked, onChange, label, disabled = false }) => {
  return (
    <div className="flex items-center space-x-3 py-1 px-0.5">
      <div
        className={`relative flex items-center w-12 h-6 rounded-full p-0.5 cursor-pointer transition-all duration-300 ease-in-out ${
          checked
            ? "bg-linear-to-r from-blue-600 via-blue-500 to-indigo-500 border-2 border-cyan-400 dark:border-cyan-300"
            : "bg-linear-to-r from-gray-300 via-gray-350 to-gray-400 dark:from-gray-600 dark:via-gray-650 dark:to-gray-700"
        } ${disabled ? "opacity-50 cursor-not-allowed" : "hover:scale-105"}`}
        onClick={() => !disabled && onChange(!checked)}
        role="switch"
        aria-checked={checked}
        aria-label={label}
      >
        <LayoutGroup>
          <motion.div
            className={`relative h-5 w-5 rounded-full ${
              checked
                ? "bg-linear-to-br from-white via-blue-50 to-indigo-50"
                : "bg-linear-to-br from-white via-gray-50 to-gray-100 dark:from-gray-200 dark:via-gray-300 dark:to-gray-400"
            }`}
            layout
            transition={toggleSpring}
            style={{
              left: checked ? "auto" : "2px",
              right: checked ? "2px" : "auto",
              position: "absolute",
            }}
          />
        </LayoutGroup>
      </div>
      <label
        htmlFor={id}
        className={`text-sm font-medium select-none transition-colors duration-200 ${
          disabled
            ? "text-gray-400 dark:text-gray-500"
            : checked
              ? "text-blue-700 dark:text-blue-400 cursor-pointer"
              : "text-gray-700 dark:text-gray-300 cursor-pointer hover:text-gray-900 dark:hover:text-gray-100"
        }`}
        onClick={() => !disabled && onChange(!checked)}
      >
        {label}
      </label>
    </div>
  );
};

const BatchDepictionView = () => {
  // Input state
  const [inputText, setInputText] = useState("");

  // Processing state
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [depictions, setDepictions] = useState([]); // Array of { smiles, title, imageUrl, id }

  // Toolkit selection
  const [toolkit, setToolkit] = useState("rdkit");

  // Depiction options - Basic (shared by both toolkits)
  const [width, setWidth] = useState(300);
  const [height, setHeight] = useState(200);
  const [showCIP, setShowCIP] = useState(false);
  const [useUnicolor, setUseUnicolor] = useState(false);
  const [highlight, setHighlight] = useState("");
  const [showAtomNumbers, setShowAtomNumbers] = useState(false);
  const [hydrogenDisplay, setHydrogenDisplay] = useState("Smart");

  // Depiction options - Enhanced features (CDK only)
  const [abbreviate, setAbbreviate] = useState("off");
  const [dative, setDative] = useState("metals");
  const [multicenter, setMulticenter] = useState("provided");
  const [style, setStyle] = useState("cow");
  const [annotate, setAnnotate] = useState("none");
  const [donuts, setDonuts] = useState(false);
  const [zoom, setZoom] = useState(1.0);
  const [ratio, setRatio] = useState(1.0);
  const [flip, setFlip] = useState(false);
  const [showtitle, setShowtitle] = useState(false);
  const [perceiveRadicals, setPerceiveRadicals] = useState(false);

  // UI state
  const [copiedSmiles, setCopiedSmiles] = useState(false); // For "Copy All SMILES" button
  const [copiedSingle, setCopiedSingle] = useState(null); // Tracks ID of molecule whose SMILES was copied
  const [downloadFormat, setDownloadFormat] = useState("svg"); // svg or png
  const [showToolsSection, setShowToolsSection] = useState(false); // Toggle visibility of options

  // Per-molecule rotation state { [id]: rotationValue }
  const [rotations, setRotations] = useState({});

  // Parse SMILES (and optional titles) from input text
  const parseSmiles = (text) => {
    return text
      .split(/[\n\r]+/) // Split by new lines
      .map((line) => line.trim()) // Trim whitespace
      .filter((line) => line.length > 0 && !line.startsWith("#")); // Remove empty lines and comments
  };

  // Regenerate depiction URLs when options change (or called manually)
  const regenerateDepictions = () => {
    if (depictions.length === 0) return;

    console.log(`Regenerating depictions with toolkit: ${toolkit}`);

    setDepictions((currentDepictions) =>
      currentDepictions.map((dep) => {
        const rotation = rotations[dep.id] || 0;

        // Use appropriate endpoint based on toolkit
        if (toolkit === "cdk") {
          // CDK uses enhanced endpoint with all features
          if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
            console.error("depictService.get2DDepictionUrlEnhanced is not available.");
            return dep;
          }

          const options = {
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            hydrogen_display: hydrogenDisplay,
            abbreviate,
            dative,
            multicenter,
            style,
            annotate,
            donuts,
            zoom,
            ratio,
            flip,
            showtitle,
            perceive_radicals: perceiveRadicals,
          };
          const updatedImageUrl = depictService.get2DDepictionUrlEnhanced(dep.smiles, options);
          return { ...dep, imageUrl: updatedImageUrl };
        } else {
          // RDKit uses basic endpoint
          if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
            console.error("depictService.get2DDepictionUrl is not available.");
            return dep;
          }

          const options = {
            toolkit: "rdkit",
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
          };
          const updatedImageUrl = depictService.get2DDepictionUrl(dep.smiles, options);
          return { ...dep, imageUrl: updatedImageUrl };
        }
      })
    );
  };

  // Handle toolkit change and regenerate
  const handleToolkitChange = (newToolkit) => {
    setToolkit(newToolkit);
    // Regenerate depictions if results already exist
    if (depictions.length > 0) {
      // We'll regenerate after state updates
      setTimeout(() => {
        regenerateDepictions();
      }, 0);
    }
  };

  // Switch between available toolkits
  const handleSwitchToolkit = () => {
    handleToolkitChange(toolkit === "rdkit" ? "cdk" : "rdkit");
  };

  // Handle text input change, toggle tools visibility
  const handleInputChange = (value) => {
    setInputText(value);
    // Show tools section automatically if there's input
    if (value.trim() && !showToolsSection) {
      setShowToolsSection(true);
    } else if (!value.trim() && showToolsSection) {
      // Keep tools section visible even if input is cleared for easier option changes
      // setShowToolsSection(false);
    }
  };

  // --- Removed useEffect hook for automatic regeneration on option changes ---

  // Handle rotation change for a specific molecule
  const handleRotationChange = (id, value) => {
    const rotation = Number(value);
    // Update rotation state for the specific molecule
    setRotations((prev) => ({ ...prev, [id]: rotation }));

    // Update the image URL immediately for this molecule
    setDepictions((currentDepictions) =>
      currentDepictions.map((dep) => {
        if (dep.id === id) {
          // Use appropriate endpoint based on toolkit
          if (toolkit === "cdk") {
            if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
              console.error("depictService.get2DDepictionUrlEnhanced is not available.");
              return dep;
            }

            const options = {
              width,
              height,
              rotate: rotation,
              CIP: showCIP,
              unicolor: useUnicolor,
              highlight: highlight || undefined,
              showAtomNumbers,
              hydrogen_display: hydrogenDisplay,
              abbreviate,
              dative,
              multicenter,
              style,
              annotate,
              donuts,
              zoom,
              ratio,
              flip,
              showtitle,
              perceive_radicals: perceiveRadicals,
            };
            const updatedImageUrl = depictService.get2DDepictionUrlEnhanced(dep.smiles, options);
            return { ...dep, imageUrl: updatedImageUrl };
          } else {
            if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
              console.error("depictService.get2DDepictionUrl is not available.");
              return dep;
            }

            const options = {
              toolkit: "rdkit",
              width,
              height,
              rotate: rotation,
              CIP: showCIP,
              unicolor: useUnicolor,
              highlight: highlight || undefined,
              showAtomNumbers,
            };
            const updatedImageUrl = depictService.get2DDepictionUrl(dep.smiles, options);
            return { ...dep, imageUrl: updatedImageUrl };
          }
        }
        return dep;
      })
    );
  };

  // Handle form submission to generate initial depictions
  const handleSubmit = async (e) => {
    e.preventDefault();

    const smilesList = parseSmiles(inputText);
    if (smilesList.length === 0) {
      setError("Please enter at least one valid SMILES string.");
      setDepictions([]);
      return;
    }

    setLoading(true);
    setError(null);
    setDepictions([]);
    setRotations({}); // Reset rotations

    try {
      const results = [];
      const maxToProcess = 50; // Limit processing
      const processingList = smilesList.slice(0, maxToProcess);

      if (smilesList.length > maxToProcess) {
        setError(`Processing limited to the first ${maxToProcess} molecules.`);
      }

      const initialRotations = {};
      for (let i = 0; i < processingList.length; i++) {
        const smilesLine = processingList[i];

        // Extract SMILES and Title with CXSMILES support
        let smiles;
        let title;

        // Check if this is CXSMILES format (contains |...|)
        const pipeIndex = smilesLine.indexOf("|");
        if (pipeIndex > -1) {
          // Find the closing pipe
          const closingPipeIndex = smilesLine.indexOf("|", pipeIndex + 1);
          if (closingPipeIndex > -1) {
            // CXSMILES format: everything up to and including the closing | is the SMILES
            smiles = smilesLine.substring(0, closingPipeIndex + 1).trim();
            // Title is everything after the closing | (if any)
            const remainingText = smilesLine.substring(closingPipeIndex + 1).trim();
            title = remainingText.length > 0 ? remainingText : `Molecule ${i + 1}`;
          } else {
            // Malformed CXSMILES (opening | but no closing |), treat as error or use fallback
            console.warn(`Malformed CXSMILES (missing closing |): ${smilesLine}`);
            const spaceIndex = smilesLine.search(/[\s\t]/);
            smiles = spaceIndex > 0 ? smilesLine.substring(0, spaceIndex) : smilesLine;
            title =
              spaceIndex > 0 ? smilesLine.substring(spaceIndex + 1).trim() : `Molecule ${i + 1}`;
          }
        } else {
          // Regular SMILES format: first space separates SMILES from title
          const spaceIndex = smilesLine.search(/[\s\t]/);
          smiles = spaceIndex > 0 ? smilesLine.substring(0, spaceIndex) : smilesLine;
          title =
            spaceIndex > 0 ? smilesLine.substring(spaceIndex + 1).trim() : `Molecule ${i + 1}`;
        }

        // Basic SMILES validation (can be improved)
        if (!smiles || smiles.length < 1) {
          console.warn(`Skipping invalid SMILES line: ${smilesLine}`);
          continue;
        }

        const id = `mol-${i}-${Date.now()}`; // Unique ID for key and state

        initialRotations[id] = 0; // Set initial rotation

        let imageUrl;

        // Use appropriate endpoint based on toolkit
        if (toolkit === "cdk") {
          if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
            console.error("depictService.get2DDepictionUrlEnhanced is not available.");
            setError("Depiction service is not configured correctly.");
            return;
          }

          const options = {
            width,
            height,
            rotate: 0,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            hydrogen_display: hydrogenDisplay,
            abbreviate,
            dative,
            multicenter,
            style,
            annotate,
            donuts,
            zoom,
            ratio,
            flip,
            showtitle,
            perceive_radicals: perceiveRadicals,
          };
          imageUrl = depictService.get2DDepictionUrlEnhanced(smiles, options);
        } else {
          if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
            console.error("depictService.get2DDepictionUrl is not available.");
            setError("Depiction service is not configured correctly.");
            return;
          }

          const options = {
            toolkit: "rdkit",
            width,
            height,
            rotate: 0,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
          };
          imageUrl = depictService.get2DDepictionUrl(smiles, options);
        }

        results.push({ smiles, title, imageUrl, id });
      }

      setRotations(initialRotations); // Set all initial rotations at once
      setDepictions(results); // Set all results at once

      if (results.length === 0 && !error) {
        setError("No valid SMILES found to process.");
      }
    } catch (err) {
      console.error("Error generating depictions:", err);
      setError(getErrorMessage("depict", err));
      setDepictions([]);
      setRotations({});
    } finally {
      setLoading(false);
    }
  };

  // Download all depictions as a ZIP file
  // NOTE: This requires the 'jszip' library to be installed (npm install jszip)
  const downloadAllDepictions = async () => {
    if (depictions.length === 0) return;

    setLoading(true); // Show loading indicator during zip creation
    setError(null);

    try {
      // Dynamically import JSZip only when needed
      const JSZip = (await import("jszip")).default;
      const zip = new JSZip();

      // Fetch all depiction blobs concurrently
      const fetchPromises = depictions.map(async (depiction) => {
        const safeTitle =
          depiction.title.replace(/[^a-z0-9]/gi, "_").toLowerCase() || `molecule_${depiction.id}`;
        const filename = `${safeTitle}.${downloadFormat}`; // Use selected format
        const rotation = rotations[depiction.id] || 0; // Get current rotation

        let downloadUrl;

        // Use appropriate endpoint based on toolkit
        if (toolkit === "cdk") {
          if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
            throw new Error("Depiction service (enhanced) is not configured correctly.");
          }

          const options = {
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            format: downloadFormat,
            hydrogen_display: hydrogenDisplay,
            abbreviate,
            dative,
            multicenter,
            style,
            annotate,
            donuts,
            zoom,
            ratio,
            flip,
            showtitle,
            perceive_radicals: perceiveRadicals,
          };
          downloadUrl = depictService.get2DDepictionUrlEnhanced(depiction.smiles, options);
        } else {
          if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
            throw new Error("Depiction service is not configured correctly.");
          }

          const options = {
            toolkit: "rdkit",
            width,
            height,
            rotate: rotation,
            CIP: showCIP,
            unicolor: useUnicolor,
            highlight: highlight || undefined,
            showAtomNumbers,
            format: downloadFormat,
          };
          downloadUrl = depictService.get2DDepictionUrl(depiction.smiles, options);
        }

        try {
          const response = await fetch(downloadUrl);
          if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);
          const blob = await response.blob();
          zip.file(filename, blob); // Add file to zip
        } catch (fetchErr) {
          console.error(`Failed to fetch or add ${filename}:`, fetchErr);
          // Add an error file to the zip instead
          zip.file(
            `${filename}.error.txt`,
            `Failed to download image for ${depiction.smiles}. Error: ${fetchErr.message}`
          );
        }
      });

      await Promise.all(fetchPromises); // Wait for all fetches

      // Generate zip file content
      const content = await zip.generateAsync({ type: "blob" });
      const downloadUrl = URL.createObjectURL(content);

      // Trigger download
      const a = document.createElement("a");
      a.href = downloadUrl;
      a.download = `depictions_${toolkit}_${downloadFormat}.zip`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(downloadUrl); // Clean up
    } catch (err) {
      console.error("Error creating zip file:", err);
      setError(getErrorMessage("depict", err));
    } finally {
      setLoading(false); // Hide loading indicator
    }
  };

  // Copy all SMILES strings to clipboard
  const copyAllSmiles = () => {
    if (depictions.length === 0 || !navigator.clipboard) return;
    const smilesText = depictions.map((d) => d.smiles).join("\n");
    navigator.clipboard
      .writeText(smilesText)
      .then(() => {
        setCopiedSmiles(true);
        setTimeout(() => setCopiedSmiles(false), 2000);
      })
      .catch((err) => {
        console.error("Failed to copy SMILES:", err);
        setError("Failed to copy SMILES to clipboard.");
      });
  };

  // Copy a single SMILES string
  const copySingleSmiles = (smiles, id) => {
    if (!navigator.clipboard) return;
    navigator.clipboard
      .writeText(smiles)
      .then(() => {
        setCopiedSingle(id);
        setTimeout(() => setCopiedSingle(null), 1500);
      })
      .catch((err) => {
        console.error("Failed to copy SMILES:", err);
        setError("Failed to copy SMILES.");
      });
  };

  // Download a single depiction image
  const downloadSingleDepiction = async (depiction) => {
    try {
      setLoading(true); // Show loading indicator while downloading
      const rotation = rotations[depiction.id] || 0;

      let url;

      // Use appropriate endpoint based on toolkit
      if (toolkit === "cdk") {
        if (!depictService || typeof depictService.get2DDepictionUrlEnhanced !== "function") {
          console.error("depictService.get2DDepictionUrlEnhanced is not available.");
          setError("Depiction service is not configured correctly.");
          return;
        }

        const options = {
          width,
          height,
          rotate: rotation,
          CIP: showCIP,
          unicolor: useUnicolor,
          highlight: highlight || undefined,
          showAtomNumbers,
          format: downloadFormat,
          hydrogen_display: hydrogenDisplay,
          abbreviate,
          dative,
          multicenter,
          style,
          annotate,
          donuts,
          zoom,
          ratio,
          flip,
          showtitle,
          perceive_radicals: perceiveRadicals,
        };
        url = depictService.get2DDepictionUrlEnhanced(depiction.smiles, options);
      } else {
        if (!depictService || typeof depictService.get2DDepictionUrl !== "function") {
          console.error("depictService.get2DDepictionUrl is not available.");
          setError("Depiction service is not configured correctly.");
          return;
        }

        const options = {
          toolkit: "rdkit",
          width,
          height,
          rotate: rotation,
          CIP: showCIP,
          unicolor: useUnicolor,
          highlight: highlight || undefined,
          showAtomNumbers,
          format: downloadFormat,
        };
        url = depictService.get2DDepictionUrl(depiction.smiles, options);
      }

      // Fetch the image data as a blob
      const response = await fetch(url);
      if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);

      const blob = await response.blob();
      const downloadUrl = URL.createObjectURL(blob);

      // Trigger download using an anchor tag
      const a = document.createElement("a");
      a.href = downloadUrl;
      a.download = `${depiction.title.replace(/[^a-z0-9]/gi, "_").toLowerCase() || "molecule"}.${downloadFormat}`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);

      // Clean up the blob URL
      URL.revokeObjectURL(downloadUrl);
    } catch (err) {
      console.error("Error downloading single depiction:", err);
      setError(getErrorMessage("depict", err));
    } finally {
      setLoading(false); // Hide loading indicator
    }
  };

  // Generate example molecules and populate input
  const generateExamples = () => {
    const examples = [
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C Caffeine",
      "CC(=O)OC1=CC=CC=C1C(=O)O Aspirin",
      "CC(=O)NC1=CC=C(C=C1)O Paracetamol",
      "Cl*.Cl*.c1ccccc1-c1ccccc1 |m:1:4.5.6.7.8.9,3:10.11.12.13.14.15| Dichlorobiphenyl",
      "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O Ibuprofen",
      "COc1cc2c(cc1OC)NC(=O)C3=C2C=CN3C Emetine",
      "Cl*.Cl*.c1ccccc1-c1ccccc1 |m:1:4.5.6.7.8.9,3:10.11.12.13.14.15|",
      "C1=CC=C(C=C1)C(=O)C(=O)O Phenylglyoxylic acid",
      "C1=CC=C(C=C1)C2=CC=C(C=C2)C(=O)O Biphenyl-4-carboxylic acid",
      "C1=CC=C2C(=C1)C=CC=C2 Naphthalene",
      "CC1=CC=C(C=C1)C Toluene",
    ];
    const exampleText = examples.join("\n");
    setInputText(exampleText);
    setDepictions([]); // Clear previous depictions
    setError(null);
    setShowToolsSection(true); // Show tools when examples are loaded
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Section Card */}
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.5, ease: [0.25, 0.46, 0.45, 0.94] }}
        className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-xl dark:shadow-2xl border border-gray-200 dark:border-gray-700"
      >
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Batch 2D Depiction
        </h2>
        <p className="text-sm text-gray-600 dark:text-gray-400 mb-4">
          Generate 2D molecular depictions with RDKit or CDK. Supports SMILES and CXSMILES formats.
          CDK offers enhanced features like abbreviations and advanced styling.
        </p>

        <form onSubmit={handleSubmit} className="space-y-4">
          {/* Input Text Area */}
          <div>
            <div className="flex justify-between items-center mb-1">
              <label
                htmlFor="smiles-input-batch"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300"
              >
                Enter SMILES (one per line, optional title after space/tab)
              </label>
              {/* Load Examples Button */}
              <Button
                variant="ghost"
                type="button"
                onClick={generateExamples}
                className="text-sm font-medium text-blue-600 hover:text-blue-700 dark:text-blue-400 dark:hover:text-blue-300 focus:outline-hidden focus-visible:ring-2 focus-visible:ring-blue-500 rounded-sm"
              >
                Load Examples
              </Button>
            </div>
            {/* Textarea Styling */}
            <Textarea
              id="smiles-input-batch"
              value={inputText}
              onChange={(e) => handleInputChange(e.target.value)}
              className="w-full h-40 font-mono text-xs sm:text-sm resize-y"
              placeholder="CN1C=NC2=C1C(=O)N(C(=O)N2C)C Caffeine&#10;CCO Ethanol&#10;Cl*.Cl*.c1ccccc1-c1ccccc1 |m:1:4.5.6.7.8.9,3:10.11.12.13.14.15| Dichlorobiphenyl&#10;# Lines starting with # are ignored"
              required
              aria-required="true"
            />
            <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
              Format: SMILES/CXSMILES [space/tab] Title (optional). CXSMILES with |...| notation
              fully supported. Max 50 molecules processed.
            </p>
          </div>

          {/* Tools Section - Animated visibility */}
          <AnimatePresence>
            {showToolsSection && (
              <motion.div
                key="tools-section"
                initial={{ opacity: 0, height: 0, scale: 0.98 }}
                animate={{
                  opacity: 1,
                  height: "auto",
                  scale: 1,
                  transition: {
                    duration: 0.4,
                    ease: [0.25, 0.46, 0.45, 0.94],
                    opacity: { duration: 0.3 },
                    height: { duration: 0.4 },
                    scale: { duration: 0.3, delay: 0.1 },
                  },
                }}
                exit={{
                  opacity: 0,
                  height: 0,
                  scale: 0.98,
                  transition: {
                    duration: 0.3,
                    ease: "easeIn",
                  },
                }}
                className="space-y-4 overflow-hidden border-t border-gray-200 dark:border-gray-700 pt-4"
              >
                <h3 className="text-md font-medium text-gray-700 dark:text-gray-300">
                  Depiction Options
                </h3>
                {/* Tools Grid */}
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4 px-0.5">
                  {/* Toolkit Select */}
                  <div>
                    <label
                      htmlFor="toolkit-select-batch"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Toolkit
                    </label>
                    <div className="flex items-center">
                      <Select value={toolkit} onValueChange={handleToolkitChange}>
                        <SelectTrigger id="toolkit-select-batch" className="w-full">
                          <SelectValue placeholder="Select toolkit" />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="rdkit">RDKit</SelectItem>
                          <SelectItem value="cdk">CDK (Enhanced)</SelectItem>
                        </SelectContent>
                      </Select>
                      <Button
                        variant="ghost"
                        type="button"
                        onClick={handleSwitchToolkit}
                        className="ml-2 p-2 bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 text-gray-600 dark:text-gray-300 rounded-md border border-gray-300 dark:border-gray-600 shadow-xs"
                        title="Switch toolkit"
                      >
                        <ArrowLeftRight className="h-5 w-5" />
                      </Button>
                    </div>
                  </div>

                  {/* Width Input */}
                  <div>
                    <label
                      htmlFor="width-input"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Width (px)
                    </label>
                    <div className="flex items-center space-x-2">
                      <motion.div whileHover={{ scale: 1.05 }} whileTap={{ scale: 0.95 }}>
                        <Button
                          type="button"
                          variant="outline"
                          size="sm"
                          onClick={() => setWidth(Math.max(50, width - 10))}
                          className="px-3 py-2 font-bold"
                          aria-label="Decrease width"
                        >
                          −
                        </Button>
                      </motion.div>
                      <Input
                        id="width-input"
                        type="number"
                        value={width}
                        onChange={(e) =>
                          setWidth(Math.max(50, Math.min(2000, Number(e.target.value))))
                        }
                        min={50}
                        max={2000}
                        step={10}
                        className="flex-1 text-center [appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none"
                      />
                      <motion.div whileHover={{ scale: 1.05 }} whileTap={{ scale: 0.95 }}>
                        <Button
                          type="button"
                          variant="outline"
                          size="sm"
                          onClick={() => setWidth(Math.min(2000, width + 10))}
                          className="px-3 py-2 font-bold"
                          aria-label="Increase width"
                        >
                          +
                        </Button>
                      </motion.div>
                    </div>
                  </div>

                  {/* Height Input */}
                  <div>
                    <label
                      htmlFor="height-input"
                      className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                    >
                      Height (px)
                    </label>
                    <div className="flex items-center space-x-2">
                      <motion.div whileHover={{ scale: 1.05 }} whileTap={{ scale: 0.95 }}>
                        <Button
                          type="button"
                          variant="outline"
                          size="sm"
                          onClick={() => setHeight(Math.max(50, height - 10))}
                          className="px-3 py-2 font-bold"
                          aria-label="Decrease height"
                        >
                          −
                        </Button>
                      </motion.div>
                      <Input
                        id="height-input"
                        type="number"
                        value={height}
                        onChange={(e) =>
                          setHeight(Math.max(50, Math.min(2000, Number(e.target.value))))
                        }
                        min={50}
                        max={2000}
                        step={10}
                        className="flex-1 text-center [appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none"
                      />
                      <motion.div whileHover={{ scale: 1.05 }} whileTap={{ scale: 0.95 }}>
                        <Button
                          type="button"
                          variant="outline"
                          size="sm"
                          onClick={() => setHeight(Math.min(2000, height + 10))}
                          className="px-3 py-2 font-bold"
                          aria-label="Increase height"
                        >
                          +
                        </Button>
                      </motion.div>
                    </div>
                  </div>

                  {/* Hydrogen Display (CDK only) */}
                  {toolkit === "cdk" && (
                    <div>
                      <label
                        htmlFor="hydrogen-display-select"
                        className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                      >
                        Hydrogen Display
                      </label>
                      <Select value={hydrogenDisplay} onValueChange={setHydrogenDisplay}>
                        <SelectTrigger id="hydrogen-display-select" className="w-full">
                          <SelectValue placeholder="Select display" />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="Smart">Smart</SelectItem>
                          <SelectItem value="Provided">Provided</SelectItem>
                          <SelectItem value="Minimal">Minimal</SelectItem>
                          <SelectItem value="Explicit">Explicit</SelectItem>
                          <SelectItem value="Stereo">Stereo</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                  )}
                </div>

                {/* CDK Enhanced Options - Animated */}
                <AnimatePresence>
                  {toolkit === "cdk" && (
                    <motion.div
                      initial={{ opacity: 0, height: 0, y: -10 }}
                      animate={{
                        opacity: 1,
                        height: "auto",
                        y: 0,
                        transition: {
                          duration: 0.5,
                          ease: [0.25, 0.46, 0.45, 0.94],
                          opacity: { duration: 0.4, delay: 0.1 },
                          height: { duration: 0.5 },
                          y: { duration: 0.4, delay: 0.1 },
                        },
                      }}
                      exit={{
                        opacity: 0,
                        height: 0,
                        y: -10,
                        transition: {
                          duration: 0.3,
                          ease: "easeIn",
                        },
                      }}
                      className="space-y-4 overflow-hidden"
                    >
                      <div className="border-t border-blue-200 dark:border-blue-800 pt-4 mt-2">
                        <h4 className="text-sm font-medium text-blue-700 dark:text-blue-400 mb-3 flex items-center">
                          <span className="mr-2">✨</span>
                          Enhanced CDK Features
                        </h4>

                        <div className="grid grid-cols-1 md:grid-cols-3 gap-4 px-0.5">
                          {/* Style Preset */}
                          <div>
                            <label
                              htmlFor="style-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Color Scheme
                            </label>
                            <Select value={style} onValueChange={setStyle}>
                              <SelectTrigger id="style-select" className="w-full">
                                <SelectValue placeholder="Select style" />
                              </SelectTrigger>
                              <SelectContent>
                                <SelectItem value="cow">Color on White</SelectItem>
                                <SelectItem value="cob">Color on Black</SelectItem>
                                <SelectItem value="cot">Color on Transparent</SelectItem>
                                <SelectItem value="bow">Black on White</SelectItem>
                                <SelectItem value="bot">Black on Transparent</SelectItem>
                                <SelectItem value="wob">White on Black</SelectItem>
                                <SelectItem value="nob">Neon on Black</SelectItem>
                              </SelectContent>
                            </Select>
                          </div>

                          {/* Abbreviations */}
                          <div>
                            <label
                              htmlFor="abbreviate-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Abbreviations
                            </label>
                            <Select value={abbreviate} onValueChange={setAbbreviate}>
                              <SelectTrigger id="abbreviate-select" className="w-full">
                                <SelectValue placeholder="Select abbreviation" />
                              </SelectTrigger>
                              <SelectContent>
                                <SelectItem value="off">None</SelectItem>
                                <SelectItem value="groups">
                                  Functional Groups (Ph, Me, Et)
                                </SelectItem>
                                <SelectItem value="reagents">Reagents (THF, DMF)</SelectItem>
                                <SelectItem value="on">Both Groups & Reagents</SelectItem>
                              </SelectContent>
                            </Select>
                          </div>

                          {/* Annotations */}
                          <div>
                            <label
                              htmlFor="annotate-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Annotations
                            </label>
                            <Select value={annotate} onValueChange={setAnnotate}>
                              <SelectTrigger id="annotate-select" className="w-full">
                                <SelectValue placeholder="Select annotation" />
                              </SelectTrigger>
                              <SelectContent>
                                <SelectItem value="none">None</SelectItem>
                                <SelectItem value="number">Atom Numbers</SelectItem>
                                <SelectItem value="cip">CIP Labels (R/S, E/Z)</SelectItem>
                                <SelectItem value="mapidx">Atom Mapping</SelectItem>
                                <SelectItem value="colmap">Color-Coded Mapping</SelectItem>
                                <SelectItem value="bondnumber">Bond Numbers</SelectItem>
                                <SelectItem value="atomvalue">Atom Values</SelectItem>
                              </SelectContent>
                            </Select>
                          </div>

                          {/* Dative Bonds */}
                          <div>
                            <label
                              htmlFor="dative-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Dative Bonds
                            </label>
                            <Select value={dative} onValueChange={setDative}>
                              <SelectTrigger id="dative-select" className="w-full">
                                <SelectValue placeholder="Select dative bonds" />
                              </SelectTrigger>
                              <SelectContent>
                                <SelectItem value="never">Never</SelectItem>
                                <SelectItem value="metals">Metal Complexes</SelectItem>
                                <SelectItem value="always">Always</SelectItem>
                              </SelectContent>
                            </Select>
                          </div>

                          {/* Multicenter Bonds */}
                          <div>
                            <label
                              htmlFor="multicenter-select"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Multicenter Bonds
                            </label>
                            <Select value={multicenter} onValueChange={setMulticenter}>
                              <SelectTrigger id="multicenter-select" className="w-full">
                                <SelectValue placeholder="Select multicenter" />
                              </SelectTrigger>
                              <SelectContent>
                                <SelectItem value="provided">As Provided</SelectItem>
                                <SelectItem value="dative">Dative Arrows</SelectItem>
                                <SelectItem value="dashed">Dashed Lines</SelectItem>
                                <SelectItem value="dashed_neutral">Dashed (Neutral)</SelectItem>
                                <SelectItem value="hidden">Hidden</SelectItem>
                                <SelectItem value="hidden_neutral">Hidden (Neutral)</SelectItem>
                              </SelectContent>
                            </Select>
                          </div>

                          {/* Zoom */}
                          <div>
                            <label
                              htmlFor="zoom-input"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Zoom: {zoom.toFixed(1)}x
                            </label>
                            <input
                              id="zoom-input"
                              type="range"
                              value={zoom}
                              onChange={(e) => setZoom(Number(e.target.value))}
                              min="0.5"
                              max="3.0"
                              step="0.1"
                              className="w-full h-2 bg-gray-200 dark:bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500 dark:accent-blue-400"
                            />
                          </div>

                          {/* Ratio */}
                          <div>
                            <label
                              htmlFor="ratio-input"
                              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                            >
                              Bond Thickness: {ratio.toFixed(1)}x
                            </label>
                            <input
                              id="ratio-input"
                              type="range"
                              value={ratio}
                              onChange={(e) => setRatio(Number(e.target.value))}
                              min="0.5"
                              max="2.0"
                              step="0.1"
                              className="w-full h-2 bg-gray-200 dark:bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500 dark:accent-blue-400"
                            />
                          </div>
                        </div>

                        {/* CDK Advanced Toggles */}
                        <div className="flex flex-wrap gap-6 pt-4 border-t border-gray-200 dark:border-gray-700 mt-4">
                          <ToggleSwitch
                            id="flip"
                            checked={flip}
                            onChange={setFlip}
                            label="Flip structure horizontally"
                          />

                          <ToggleSwitch
                            id="showtitle"
                            checked={showtitle}
                            onChange={setShowtitle}
                            label="Show titles in depiction"
                          />

                          <ToggleSwitch
                            id="perceive-radicals"
                            checked={perceiveRadicals}
                            onChange={setPerceiveRadicals}
                            label="Perceive radicals"
                          />
                        </div>
                      </div>
                    </motion.div>
                  )}
                </AnimatePresence>

                {/* Highlight Input */}
                <div>
                  <label
                    htmlFor="highlight-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
                  >
                    Highlight Substructure (SMARTS - optional)
                  </label>
                  <Input
                    id="highlight-input"
                    type="text"
                    value={highlight}
                    onChange={(e) => setHighlight(e.target.value)}
                    placeholder="e.g., c1ccccc1"
                    className="w-full"
                  />
                </div>

                {/* Checkboxes */}
                <div className="flex flex-wrap gap-6 pt-3">
                  <ToggleSwitch
                    id="cip"
                    checked={showCIP}
                    onChange={setShowCIP}
                    label="Show CIP (R/S, E/Z)"
                  />
                  <ToggleSwitch
                    id="unicolor"
                    checked={useUnicolor}
                    onChange={setUseUnicolor}
                    label="Black & white"
                  />
                  <ToggleSwitch
                    id="atomNumbers"
                    checked={showAtomNumbers}
                    onChange={setShowAtomNumbers}
                    label="Show atom numbers"
                  />
                  {/* Aromatic Donuts Toggle (CDK only) */}
                  {toolkit === "cdk" && (
                    <motion.div
                      initial={{ opacity: 0, scale: 0.9 }}
                      animate={{ opacity: 1, scale: 1 }}
                      exit={{ opacity: 0, scale: 0.9 }}
                      transition={{ duration: 0.2 }}
                    >
                      <ToggleSwitch
                        id="donuts"
                        checked={donuts}
                        onChange={setDonuts}
                        label="Aromatic rings (circles)"
                      />
                    </motion.div>
                  )}
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Action Buttons Row */}
          <div className="flex flex-wrap items-center gap-3 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Generate Button */}
            <motion.div
              whileHover={!inputText.trim() || loading ? {} : { scale: 1.02, y: -1 }}
              whileTap={!inputText.trim() || loading ? {} : { scale: 0.98 }}
              transition={{ type: "spring", stiffness: 400, damping: 17 }}
            >
              <Button
                type="submit"
                disabled={!inputText.trim() || loading}
                className={cn(
                  "flex items-center justify-center",
                  !inputText.trim() || loading ? "" : "shadow-lg hover:shadow-xl"
                )}
              >
                <Image className="mr-2 h-5 w-5" />
                {loading ? "Generating..." : "Generate Depictions"}
              </Button>
            </motion.div>

            {/* Regenerate Button - only shown when depictions exist and not loading */}
            {depictions.length > 0 && !loading && (
              <motion.div
                initial={{ opacity: 0, scale: 0.9 }}
                animate={{ opacity: 1, scale: 1 }}
                exit={{ opacity: 0, scale: 0.9 }}
                whileHover={{ scale: 1.02, y: -1 }}
                whileTap={{ scale: 0.98 }}
                transition={{ type: "spring", stiffness: 400, damping: 17 }}
              >
                <Button
                  variant="secondary"
                  onClick={() => regenerateDepictions()}
                  className="flex items-center justify-center shadow-md hover:shadow-lg"
                  title="Apply current options to all depictions"
                >
                  <Image className="mr-2 h-4 w-4" />
                  Update All
                </Button>
              </motion.div>
            )}

            {/* Download and Copy buttons - only shown when depictions exist */}
            {depictions.length > 0 && !loading && (
              <>
                {/* Download Format Selector and Button */}
                <div className="flex items-center shadow-md rounded-lg overflow-hidden px-0.5">
                  {/* Select Styling */}
                  <Select value={downloadFormat} onValueChange={setDownloadFormat}>
                    <SelectTrigger
                      className="h-full border-r text-sm rounded-r-none"
                      aria-label="Download format"
                    >
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="svg">SVG</SelectItem>
                      <SelectItem value="png">PNG</SelectItem>
                    </SelectContent>
                  </Select>
                  {/* Download Button with Animation */}
                  <motion.div
                    whileHover={{ scale: 1.02 }}
                    whileTap={{ scale: 0.98 }}
                    transition={{ type: "spring", stiffness: 400, damping: 17 }}
                  >
                    <Button
                      type="button"
                      onClick={downloadAllDepictions}
                      className="bg-green-600 hover:bg-green-700 dark:bg-green-500 dark:hover:bg-green-600 text-white flex items-center text-sm"
                      title="Download all depictions as a ZIP file"
                    >
                      <Download className="mr-1.5 h-4 w-4" />
                      Download All (.zip)
                    </Button>
                  </motion.div>
                </div>

                {/* Copy All SMILES Button */}
                <motion.div
                  whileHover={{ scale: 1.02 }}
                  whileTap={{ scale: 0.98 }}
                  transition={{ type: "spring", stiffness: 400, damping: 17 }}
                >
                  <Button
                    type="button"
                    variant={copiedSmiles ? "secondary" : "outline"}
                    onClick={copyAllSmiles}
                    className={cn(
                      "flex items-center text-sm shadow-md hover:shadow-lg",
                      copiedSmiles &&
                        "bg-green-100 dark:bg-green-700 text-green-700 dark:text-green-200"
                    )}
                    title="Copy all input SMILES to clipboard"
                  >
                    <motion.div
                      animate={{ rotate: copiedSmiles ? [0, -10, 10, -10, 0] : 0 }}
                      transition={{ duration: 0.5 }}
                    >
                      {copiedSmiles ? (
                        <Check className="mr-1.5 h-5 w-5" />
                      ) : (
                        <Clipboard className="mr-1.5 h-5 w-5" />
                      )}
                    </motion.div>
                    {copiedSmiles ? "Copied!" : "Copy All SMILES"}
                  </Button>
                </motion.div>
              </>
            )}
          </div>
        </form>
      </motion.div>

      {/* Loading State */}
      {loading && depictions.length === 0 && <ToolSkeleton variant="molecule" />}

      {/* Error Display */}
      {error && !loading && (
        <GlassErrorCard
          message={error}
          onRetry={() => {
            setError(null);
          }}
        />
      )}

      {/* Results Grid */}
      {depictions.length > 0 && !loading && (
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, ease: [0.25, 0.46, 0.45, 0.94] }}
          className="space-y-4"
        >
          {/* Results Header */}
          <motion.h3
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ duration: 0.4, delay: 0.1 }}
            className="text-lg font-semibold text-gray-800 dark:text-gray-200 flex items-center"
          >
            <span className="mr-2">🧪</span>
            Generated Depictions ({depictions.length})
          </motion.h3>

          {/* Grid Container with Animation */}
          <motion.div
            className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-4" // Responsive grid
            variants={resultsContainerVariants}
            initial="hidden"
            animate="visible"
          >
            {depictions.map((depiction) => (
              // Individual Card with Animation
              <motion.div
                key={depiction.id}
                className="bg-white dark:bg-gray-800 rounded-xl overflow-hidden shadow-lg dark:shadow-2xl border border-gray-200 dark:border-gray-700 flex flex-col hover:shadow-2xl dark:hover:shadow-3xl transition-shadow duration-300"
                variants={depictionCardVariant}
                layout
                whileHover={{ y: -4, transition: { duration: 0.2 } }}
              >
                {/* Image container - Kept white background for consistency */}
                <div className="p-2 grow">
                  {/* Ensure depiction image background is consistent */}
                  <div className="bg-white rounded-md p-2 flex items-center justify-center h-48">
                    <img
                      src={depiction.imageUrl}
                      alt={depiction.title}
                      className="max-w-full max-h-full object-contain"
                      loading="lazy" // Lazy load images
                      onError={(e) => {
                        console.error(
                          `Error loading image for ${depiction.title}: ${depiction.imageUrl}`
                        );
                        e.target.onerror = null;
                        // Use neutral fallback SVG
                        e.target.src =
                          "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj48dGV4dCB4PSI1MCIgeT0iNTAiIGRvbWluYW50LWJhc2VsaW5lPSJtaWRkbGUiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZvbnQtc2l6ZT0iMTAiIGZpbGw9IiM4ODg4ODg4Ij5FcnJvciBsb2FkaW5nPC90ZXh0Pjwvc3ZnPg==";
                        e.target.classList.add("p-4"); // Add padding to text
                      }}
                    />
                  </div>
                </div>

                {/* Info Section */}
                <div className="px-3 py-2 border-t border-gray-200 dark:border-gray-700">
                  <h4
                    className="font-medium text-gray-800 dark:text-gray-100 text-sm truncate"
                    title={depiction.title}
                  >
                    {depiction.title}
                  </h4>
                  <p
                    className="text-gray-500 dark:text-gray-400 text-xs truncate font-mono"
                    title={depiction.smiles}
                  >
                    {depiction.smiles}
                  </p>
                </div>

                {/* Rotation Slider */}
                <div className="px-3 py-2 bg-gray-50 dark:bg-gray-900 border-t border-gray-200 dark:border-gray-700">
                  <label
                    htmlFor={`rotate-${depiction.id}`}
                    className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1"
                  >
                    Rotate: {rotations[depiction.id] || 0}°
                  </label>
                  {/* Slider Styling */}
                  <input
                    id={`rotate-${depiction.id}`}
                    type="range"
                    value={rotations[depiction.id] || 0}
                    onChange={(e) => handleRotationChange(depiction.id, e.target.value)}
                    min="0"
                    max="359"
                    step="1"
                    // Enhanced slider styling
                    className="w-full h-2 bg-gray-200 dark:bg-gray-700 rounded-lg appearance-none cursor-pointer range-sm accent-blue-500 dark:accent-blue-400 focus:outline-hidden focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-900 focus:ring-blue-500"
                    aria-label={`Rotation for ${depiction.title}`}
                  />
                  {/* Slider Labels */}
                  <div className="flex justify-between text-xs text-gray-400 dark:text-gray-500 mt-1">
                    <span>0°</span>
                    <span>359°</span>
                  </div>
                </div>

                {/* Action Buttons */}
                <div className="px-3 py-2 bg-linear-to-r from-gray-50 to-gray-100 dark:from-gray-700 dark:to-gray-800 border-t border-gray-200 dark:border-gray-600 flex justify-end space-x-2">
                  {/* Download Button */}
                  <motion.div
                    whileHover={{ scale: 1.1, y: -2 }}
                    whileTap={{ scale: 0.95 }}
                    transition={{ type: "spring", stiffness: 400, damping: 17 }}
                  >
                    <Button
                      variant="ghost"
                      size="icon"
                      onClick={() => downloadSingleDepiction(depiction)}
                      className="rounded-lg text-gray-600 dark:text-gray-400 hover:text-blue-600 dark:hover:text-blue-400 hover:bg-blue-50 dark:hover:bg-blue-900/30 shadow-xs hover:shadow-md"
                      title={`Download ${downloadFormat.toUpperCase()}`}
                      aria-label={`Download ${depiction.title} as ${downloadFormat.toUpperCase()}`}
                      disabled={loading}
                    >
                      <Download className="h-5 w-5" />
                    </Button>
                  </motion.div>

                  {/* Copy SMILES Button */}
                  <motion.div
                    whileHover={{ scale: 1.1, y: -2 }}
                    whileTap={{ scale: 0.95 }}
                    transition={{ type: "spring", stiffness: 400, damping: 17 }}
                  >
                    <Button
                      variant="ghost"
                      size="icon"
                      onClick={() => copySingleSmiles(depiction.smiles, depiction.id)}
                      className="rounded-lg text-gray-600 dark:text-gray-400 hover:text-green-600 dark:hover:text-green-400 hover:bg-green-50 dark:hover:bg-green-900/30 relative shadow-xs hover:shadow-md"
                      title="Copy SMILES"
                      aria-label={`Copy SMILES for ${depiction.title}`}
                    >
                      {/* Animated Check/Clipboard Icon */}
                      <AnimatePresence mode="wait">
                        {copiedSingle === depiction.id ? (
                          <motion.div
                            key="check"
                            initial={{ scale: 0, rotate: -180 }}
                            animate={{ scale: 1, rotate: 0 }}
                            exit={{ scale: 0, rotate: 180 }}
                            transition={{
                              type: "spring",
                              stiffness: 500,
                              damping: 25,
                              duration: 0.3,
                            }}
                            className="flex items-center justify-center"
                          >
                            <Check className="h-5 w-5 text-green-500 dark:text-green-400" />
                          </motion.div>
                        ) : (
                          <motion.div
                            key="clipboard"
                            initial={{ scale: 1, rotate: 0 }}
                            exit={{ scale: 0, rotate: -180 }}
                            transition={{ duration: 0.2 }}
                            className="flex items-center justify-center"
                          >
                            <Clipboard className="h-5 w-5" />
                          </motion.div>
                        )}
                      </AnimatePresence>
                    </Button>
                  </motion.div>
                </div>
              </motion.div>
            ))}
          </motion.div>
        </motion.div>
      )}

      {/* Initial Placeholder */}
      {!depictions.length && !loading && !error && (
        // Initial state card styling
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg flex flex-col items-center justify-center text-center min-h-[300px] border border-gray-200 dark:border-gray-700">
          <Image className="h-16 w-16 text-gray-400 dark:text-gray-600 mb-4" />
          <p className="text-gray-500 dark:text-gray-400">
            Enter SMILES or CXSMILES strings (one per line, optionally with titles) and click
            "Generate Depictions".
          </p>
          <Button
            variant="ghost"
            type="button"
            onClick={generateExamples}
            className="mt-4 text-sm font-medium text-blue-600 hover:text-blue-700 dark:text-blue-400 dark:hover:text-blue-300 focus:outline-hidden focus-visible:ring-2 focus-visible:ring-blue-500 rounded-sm"
          >
            Or load example molecules
          </Button>
        </div>
      )}
    </div>
  );
};
export default BatchDepictionView;
