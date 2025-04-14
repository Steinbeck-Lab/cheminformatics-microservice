import React, { useState, useEffect, useRef } from "react";
import {
  HiOutlineClipboardCopy,
  HiOutlineCheck,
  HiOutlineExclamationCircle,
  HiOutlineInformationCircle,
  HiOutlineRefresh,
  HiOutlinePencil,
  HiOutlineX,
  HiOutlineCode,
  HiOutlineDocumentText,
  HiOutlineUpload,
  HiOutlineQuestionMarkCircle,
  HiOutlineArrowsExpand,
  HiOutlineChevronDown,
  HiOutlineChevronRight,
  HiOutlineSwitchHorizontal,
} from "react-icons/hi";

// Import utility functions
import {
  RINCHI_VERSIONS,
  loadRinchiModule,
  convertRxnfileToRinchi,
  generateRinchiKey,
  convertRinchiToFileText,
} from "../../utils/rinchiUtils.js";

// Tooltip component for RInChI options
const OptionTooltip = ({ content }) => (
  <div className="group relative inline-block">
    <HiOutlineQuestionMarkCircle className="h-4 w-4 ml-1 text-gray-500 dark:text-gray-400 inline-block align-text-bottom cursor-help" />
    <div className="absolute z-10 opacity-0 invisible group-hover:visible group-hover:opacity-100 transition-opacity duration-300 w-64 bg-white dark:bg-gray-800 p-2 rounded-md shadow-lg border border-gray-200 dark:border-gray-700 text-xs text-gray-700 dark:text-gray-300 -mt-1 left-6">
      {content}
    </div>
  </div>
);

// RInChI options component
const RInChIOptions = ({
  onChange,
  rinchiVersion,
  setRinchiVersion,
  forceEquilibrium,
  setForceEquilibrium,
  outputFormat,
  setOutputFormat,
}) => {
  return (
    <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
      <div className="flex items-center mb-4">
        <div className="bg-indigo-100 dark:bg-indigo-900/50 p-2 rounded-lg mr-3">
          <HiOutlineCode className="h-5 w-5 text-indigo-700 dark:text-indigo-400" />
        </div>
        <h2 className="text-lg font-bold text-gray-800 dark:text-white">
          RInChI Options
        </h2>
      </div>

      <div className="mb-4">
        <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
          RInChI Version
        </label>
        <select
          className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500"
          value={rinchiVersion}
          onChange={(e) => setRinchiVersion(e.target.value)}
        >
          {Object.entries(RINCHI_VERSIONS).map(([value, version]) => (
            <option key={value} value={value}>
              {version.label}
            </option>
          ))}
        </select>
      </div>

      <div className="space-y-4">
        <div className="flex items-center">
          <input
            id="force-equilibrium"
            type="checkbox"
            checked={forceEquilibrium}
            onChange={(e) => setForceEquilibrium(e.target.checked)}
            className="h-4 w-4 text-indigo-600 focus:ring-indigo-500 border-gray-300 rounded"
          />
          <label
            htmlFor="force-equilibrium"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Force Equilibrium
          </label>
          <OptionTooltip content="When checked, the reaction is treated as an equilibrium reaction, regardless of the arrow type used in the drawing." />
        </div>

        <div className="border-t border-gray-200 dark:border-gray-700 pt-3">
          <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2 flex items-center">
            File Format Options
            <OptionTooltip content="Configure the output format for file conversion" />
          </p>

          <div className="flex flex-col space-y-2 ml-2">
            <label className="inline-flex items-center">
              <input
                type="radio"
                name="outputFormat"
                value="RXN"
                checked={outputFormat === "RXN"}
                onChange={() => setOutputFormat("RXN")}
                className="h-4 w-4 text-indigo-600 focus:ring-indigo-500 border-gray-300"
              />
              <span className="ml-2 text-sm text-gray-700 dark:text-gray-300">
                RXN File Format
              </span>
              <OptionTooltip content="Standard MDL RXN file format for chemical reactions" />
            </label>

            <label className="inline-flex items-center">
              <input
                type="radio"
                name="outputFormat"
                value="RD"
                checked={outputFormat === "RD"}
                onChange={() => setOutputFormat("RD")}
                className="h-4 w-4 text-indigo-600 focus:ring-indigo-500 border-gray-300"
              />
              <span className="ml-2 text-sm text-gray-700 dark:text-gray-300">
                RD File Format
              </span>
              <OptionTooltip content="MDL RD file format for reaction data with properties" />
            </label>
          </div>
        </div>
      </div>
    </div>
  );
};

// Result block component for displaying RInChI, RInChIKeys, and RAuxInfo
const ResultBlock = ({
  title,
  value,
  onCopy,
  copyState,
  icon,
  collapsible = false,
}) => {
  const [isCollapsed, setIsCollapsed] = useState(false);

  return (
    <div className="mt-4 bg-white dark:bg-gray-800 p-4 rounded-lg border border-gray-200 dark:border-gray-700 shadow">
      <div className="flex justify-between items-center mb-2">
        <h3
          className={`text-sm font-medium text-gray-700 dark:text-gray-300 flex items-center ${
            collapsible
              ? "cursor-pointer hover:text-indigo-600 dark:hover:text-indigo-400"
              : ""
          }`}
          onClick={collapsible ? () => setIsCollapsed(!isCollapsed) : undefined}
        >
          {collapsible &&
            (isCollapsed ? (
              <HiOutlineChevronRight className="h-4 w-4 mr-1" />
            ) : (
              <HiOutlineChevronDown className="h-4 w-4 mr-1" />
            ))}
          {icon && <span className="mr-2">{icon}</span>}
          {title}
        </h3>
        {value && (
          <button
            onClick={onCopy}
            className="inline-flex items-center px-2 py-1 text-xs text-indigo-700 dark:text-indigo-300 bg-indigo-100 dark:bg-indigo-900/30 rounded hover:bg-indigo-200 dark:hover:bg-indigo-800/50"
          >
            {copyState ? (
              <>
                <HiOutlineCheck className="w-3.5 h-3.5 mr-1" />
                Copied!
              </>
            ) : (
              <>
                <HiOutlineClipboardCopy className="w-3.5 h-3.5 mr-1" />
                Copy
              </>
            )}
          </button>
        )}
      </div>
      {(!isCollapsed || !collapsible) && (
        <div className="font-mono text-sm overflow-x-auto bg-gray-50 dark:bg-gray-900 p-3 rounded border border-gray-200 dark:border-gray-700 break-all max-h-24 overflow-y-auto">
          {value || (
            <span className="text-gray-400 dark:text-gray-500">
              No data generated yet
            </span>
          )}
        </div>
      )}
    </div>
  );
};

const RInChIView = () => {
  // State variables
  const [rinchi, setRinchi] = useState("");
  const [rauxInfo, setRauxInfo] = useState("");
  const [longRinchiKey, setLongRinchiKey] = useState("");
  const [shortRinchiKey, setShortRinchiKey] = useState("");
  const [webRinchiKey, setWebRinchiKey] = useState("");
  const [inputRinchi, setInputRinchi] = useState("");
  const [inputRauxInfo, setInputRauxInfo] = useState("");
  const [rxnfileContent, setRxnfileContent] = useState("");
  const [outputRxnFile, setOutputRxnFile] = useState("");
  const [logMessage, setLogMessage] = useState("");
  const [forceEquilibrium, setForceEquilibrium] = useState(false);
  const [outputFormat, setOutputFormat] = useState("RXN");
  const [rinchiVersion, setRinchiVersion] = useState(
    Object.keys(RINCHI_VERSIONS).find((key) => RINCHI_VERSIONS[key].default) ||
      "1.1"
  );
  const [activeInputType, setActiveInputType] = useState("reaction"); // 'reaction', 'rxnfile', 'rinchi', 'rinchi-to-file'

  const [isLoading, setIsLoading] = useState(false);
  const [isEditorReady, setIsEditorReady] = useState(false);
  const [error, setError] = useState(null);
  const [copySuccess, setCopySuccess] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [copyModalText, setCopyModalText] = useState("");
  const [rinchiModuleLoaded, setRinchiModuleLoaded] = useState(false);
  const [retryAttempt, setRetryAttempt] = useState(0);

  // Refs
  const ketcherFrame = useRef(null);
  const copyTextRef = useRef(null);
  const fileInputRef = useRef(null);
  const messageHandlers = useRef({});
  const messageId = useRef(0);

  // Reset error when needed
  useEffect(() => {
    if (error) setError(null);
  }, [inputRinchi, rinchiVersion, activeInputType, error]);

  // Auto-select text in copy modal when it appears
  useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      copyTextRef.current.select();
    }
  }, [showCopyModal]);

  // Set up message communication with the Ketcher iframe
  useEffect(() => {
    const handleMessage = (event) => {
      // Only accept messages from our iframe
      if (
        ketcherFrame.current &&
        event.source === ketcherFrame.current.contentWindow
      ) {
        const { id, type, status, data, error } = event.data;

        // Handle response to our message
        if (id && messageHandlers.current[id]) {
          const handler = messageHandlers.current[id];

          if (status === "success") {
            handler.resolve(data);
          } else if (status === "error") {
            handler.reject(new Error(error || "Unknown error"));
          }

          // Remove the handler after it's been used
          delete messageHandlers.current[id];
        }

        // Handle initialization message
        if (type === "ketcher-ready") {
          console.log("Ketcher ready message received");
          setIsEditorReady(true);
        }
      }
    };

    // Add listener for messages
    window.addEventListener("message", handleMessage);

    // Clean up on unmount
    return () => {
      window.removeEventListener("message", handleMessage);
      messageHandlers.current = {};
    };
  }, []);

  // Function to send messages to the iframe and wait for response
  const sendMessage = (type, payload = {}) => {
    return new Promise((resolve, reject) => {
      if (!ketcherFrame.current || !ketcherFrame.current.contentWindow) {
        reject(new Error("Ketcher frame not available"));
        return;
      }

      // Generate unique ID for this message
      const id = messageId.current++;

      // Store promise handlers
      messageHandlers.current[id] = { resolve, reject };

      // Send message to iframe
      const message = { id, type, payload };

      try {
        ketcherFrame.current.contentWindow.postMessage(message, "*");
      } catch (err) {
        delete messageHandlers.current[id];
        reject(new Error(`Failed to send message: ${err.message}`));
      }

      // Set timeout to reject if no response
      setTimeout(() => {
        if (messageHandlers.current[id]) {
          delete messageHandlers.current[id];
          reject(new Error("Timeout waiting for response"));
        }
      }, 10000); // 10 second timeout
    });
  };

  // Function to communicate with ketcher through both direct access and messaging
  const executeKetcherCommand = async (command, args = []) => {
    // First try direct access method
    try {
      if (ketcherFrame.current && ketcherFrame.current.contentWindow.ketcher) {
        const ketcher = ketcherFrame.current.contentWindow.ketcher;
        if (typeof ketcher[command] === "function") {
          return await ketcher[command](...args);
        }
      }
    } catch (directError) {
      console.debug(`Direct Ketcher access failed for ${command}`, directError);
      // Fall through to postMessage method
    }

    // Fall back to message passing
    try {
      return await sendMessage("ketcher-command", { command, args });
    } catch (msgError) {
      console.error(`Message passing failed for ${command}`, msgError);
      throw new Error(`Failed to execute ${command}: ${msgError.message}`);
    }
  };

  // Initialize Ketcher iframe communication
  const initializeKetcher = () => {
    // Function to inject communication script into iframe
    const injectCommunicationScript = () => {
      try {
        const iframeWindow = ketcherFrame.current.contentWindow;
        const iframeDocument = iframeWindow.document;

        // Create a script element
        const script = iframeDocument.createElement("script");
        script.textContent = `
          // Set up message handler in the Ketcher iframe
          window.addEventListener('message', async function(event) {
            // Check source
            if (event.source !== window.parent) return;
            
            const { id, type, payload } = event.data;
            
            // Handle command execution
            if (type === 'ketcher-command') {
              try {
                const { command, args } = payload;
                
                if (!window.ketcher) {
                  throw new Error('Ketcher not initialized');
                }
                
                if (typeof window.ketcher[command] !== 'function') {
                  throw new Error(\`Command \${command} not available\`);
                }
                
                const result = await window.ketcher[command](...(args || []));
                event.source.postMessage({ id, status: 'success', data: result }, event.origin);
              } catch (error) {
                event.source.postMessage({ id, status: 'error', error: error.message }, event.origin);
              }
            }
          });
          
          // Wait for Ketcher to initialize and then notify parent
          const checkKetcher = () => {
            if (window.ketcher) {
              // Notify parent that Ketcher is ready
              window.parent.postMessage({ type: 'ketcher-ready' }, '*');
              console.log('Ketcher ready, notified parent');
            } else {
              // Check again in 100ms
              setTimeout(checkKetcher, 100);
            }
          };
          
          // Start checking for ketcher
          checkKetcher();
        `;

        // Add script to iframe
        iframeDocument.head.appendChild(script);
        console.log("Communication script injected into iframe");
      } catch (error) {
        console.error("Failed to inject communication script:", error);
      }
    };

    // Wait for iframe to load, then inject script
    if (ketcherFrame.current) {
      // Clear any previous onload
      ketcherFrame.current.onload = () => {
        console.log("Ketcher iframe loaded, injecting script");
        // Let the iframe load completely before injecting
        setTimeout(injectCommunicationScript, 500);
      };
    }
  };

  // Re-initialize when retry attempt changes
  useEffect(() => {
    initializeKetcher();
  }, [retryAttempt]);

  // Enhanced copyToClipboard function with multiple fallback methods
  const copyToClipboard = async (text = null, type = "rinchi") => {
    let textToCopy;

    switch (type) {
      case "rinchi":
        textToCopy = text || rinchi;
        break;
      case "longkey":
        textToCopy = text || longRinchiKey;
        break;
      case "shortkey":
        textToCopy = text || shortRinchiKey;
        break;
      case "webkey":
        textToCopy = text || webRinchiKey;
        break;
      case "rauxinfo":
        textToCopy = text || rauxInfo;
        break;
      case "rxnfile":
        textToCopy = text || outputRxnFile;
        break;
      default:
        textToCopy = text || rinchi;
    }

    if (!textToCopy) {
      setError(
        `No ${type
          .replace("key", " RInChIKey")
          .toUpperCase()} to copy. Generate data first.`
      );
      return;
    }

    // Try multiple clipboard copy methods in sequence
    try {
      // Method 1: Use the Clipboard API (modern browsers)
      if (navigator.clipboard && navigator.clipboard.writeText) {
        await navigator.clipboard.writeText(textToCopy);
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      }

      // Method 2: Use execCommand (older browsers)
      const textArea = document.createElement("textarea");
      textArea.value = textToCopy;

      // Make the textarea out of viewport
      textArea.style.position = "fixed";
      textArea.style.left = "-999999px";
      textArea.style.top = "-999999px";
      document.body.appendChild(textArea);

      // Select and copy
      textArea.focus();
      textArea.select();

      const successful = document.execCommand("copy");
      document.body.removeChild(textArea);

      if (successful) {
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      } else {
        throw new Error("execCommand copy failed");
      }
    } catch (err) {
      console.error("Failed to copy text:", err);

      // Method 3: Show a modal with text to copy manually
      setCopyModalText(textToCopy);
      setShowCopyModal(true);
    }
  };

  // Check RInChI module status
  useEffect(() => {
    const checkModuleStatus = async () => {
      try {
        setRinchiModuleLoaded(false);

        // Try to load the module for the current version
        await loadRinchiModule(rinchiVersion)
          .then(() => {
            setRinchiModuleLoaded(true);
            console.log(`RInChI module ${rinchiVersion} loaded successfully`);
          })
          .catch((err) => {
            console.error(
              `Failed to load RInChI module ${rinchiVersion}:`,
              err
            );
            setError(
              `Failed to load RInChI module ${rinchiVersion}. Please check your network connection.`
            );
          });
      } catch (err) {
        console.error("Error checking RInChI module status:", err);
      }
    };

    checkModuleStatus();
  }, [rinchiVersion]);

  // Check if Ketcher has a reaction
  const checkForReaction = async () => {
    try {
      // First try using the containsReaction method if available
      try {
        const ketcher = ketcherFrame.current.contentWindow.ketcher;
        if (ketcher.containsReaction) {
          return ketcher.containsReaction();
        }
      } catch (e) {
        // If direct access fails, try another approach
      }

      // Alternative: try to get the RXN file and check if it has reaction components
      const rxnFile = await executeKetcherCommand("getRxn");
      return (
        rxnFile && (rxnFile.includes("$RXN") || rxnFile.includes("$RDFILE"))
      );
    } catch (err) {
      console.error("Error checking for reaction:", err);
      return false;
    }
  };

  // Check if reaction has equilibrium arrow
  const hasEquilibriumArrow = async () => {
    try {
      // This is a simplistic check - in a real implementation you would need more
      // complex logic to detect the arrow type
      const rxnFile = await executeKetcherCommand("getRxn");
      return rxnFile && rxnFile.includes("EQUILIBRIUM");
    } catch (err) {
      console.error("Error checking for equilibrium arrow:", err);
      return false;
    }
  };

  // Generate RInChI from Ketcher
  const generateRInChI = async () => {
    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!rinchiModuleLoaded) {
      setError("RInChI module not loaded. Please wait or refresh the page.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Generating RInChI...");
    setRinchi("");
    setRauxInfo("");
    setLongRinchiKey("");
    setShortRinchiKey("");
    setWebRinchiKey("");

    try {
      // Check if there's a reaction
      const hasReaction = await checkForReaction();
      if (!hasReaction) {
        setError(
          "No reaction found. Please draw a reaction with reactants and products."
        );
        setIsLoading(false);
        return;
      }

      // Get RXN file from Ketcher
      const rxnFile = await executeKetcherCommand("getRxn");

      if (!rxnFile || rxnFile.trim() === "") {
        setError(
          "Failed to get reaction data. Please ensure you've drawn a valid reaction."
        );
        setIsLoading(false);
        return;
      }

      // Save the rxnFile for reference
      setRxnfileContent(rxnFile);

      // Check if reaction has equilibrium arrow (if not forcing equilibrium)
      let isEquilibrium = forceEquilibrium;
      if (!isEquilibrium) {
        isEquilibrium = await hasEquilibriumArrow();
      }

      // Generate RInChI from the RXN file
      await convertReactionToRinchiAndWriteResults(rxnFile, isEquilibrium);
    } catch (err) {
      console.error("Failed to generate RInChI:", err);
      setError(`Failed to generate RInChI: ${err.message}`);
      setIsLoading(false);
    }
  };

  // Convert reaction in Ketcher to RInChI
  const convertReactionToRinchiAndWriteResults = async (
    rxnFile,
    isEquilibrium = false
  ) => {
    try {
      if (!rxnFile || rxnFile.trim() === "") {
        throw new Error("No reaction data provided");
      }

      if (!rinchiModuleLoaded) {
        throw new Error(
          "RInChI module not loaded. Please wait or refresh the page."
        );
      }

      // Convert RXN to RInChI
      const result = await convertRxnfileToRinchi(
        rxnFile,
        isEquilibrium,
        rinchiVersion
      );

      if (!result) {
        throw new Error("No result returned from conversion function");
      }

      if (result.return_code !== 0) {
        throw new Error(result.error || "Failed to generate RInChI");
      }

      // Update state with results
      setRinchi(result.rinchi || "");
      setRauxInfo(result.rauxinfo || "");

      // Update log
      setLogMessage(
        `RInChI generated successfully${
          isEquilibrium ? " (equilibrium forced)" : ""
        }\n${result.error || ""}`
      );

      // Generate RInChIKeys
      if (result.rinchi) {
        await Promise.all([
          generateRinchiKeyAndUpdateState(result.rinchi, "Long"),
          generateRinchiKeyAndUpdateState(result.rinchi, "Short"),
          generateRinchiKeyAndUpdateState(result.rinchi, "Web"),
        ]);
      }

      return result;
    } catch (err) {
      console.error("Failed to generate RInChI from reaction:", err);
      setError(`Failed to generate RInChI: ${err.message}`);
      throw err;
    } finally {
      setIsLoading(false);
    }
  };

  // Handle file upload for RXN/RD file
  const handleRxnFileUpload = (event) => {
    const file = event.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      setRxnfileContent(e.target.result);
    };
    reader.readAsText(file);
  };

  // Process RXN/RD file
  const processRxnFile = async () => {
    if (!rxnfileContent.trim()) {
      setError("Please enter or upload a RXN/RD file");
      return;
    }

    if (!rinchiModuleLoaded) {
      setError("RInChI module not loaded. Please wait or refresh the page.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Processing RXN/RD file...");
    setRinchi("");
    setRauxInfo("");
    setLongRinchiKey("");
    setShortRinchiKey("");
    setWebRinchiKey("");

    try {
      // Validate the file format
      if (
        !rxnfileContent.includes("$RXN") &&
        !rxnfileContent.includes("$RDFILE")
      ) {
        setError("Invalid RXN/RD file format. Please check your input.");
        setIsLoading(false);
        return;
      }

      // Convert to RInChI
      await convertReactionToRinchiAndWriteResults(
        rxnfileContent,
        forceEquilibrium
      );

      // If editor is ready, also load the reaction in Ketcher
      if (isEditorReady) {
        try {
          await executeKetcherCommand("setMolecule", [rxnfileContent]);
          setLogMessage((prevLog) => `Reaction loaded in editor.\n${prevLog}`);
        } catch (ketcherErr) {
          console.error("Failed to load reaction in Ketcher:", ketcherErr);
          // Don't fail the whole process if just the visualization fails
        }
      }
    } catch (err) {
      console.error("Failed to process RXN/RD file:", err);
      setError(`Failed to process RXN/RD file: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Generate RInChIKey of specific type and update state
  const generateRinchiKeyAndUpdateState = async (rinchiString, keyType) => {
    if (!rinchiString) {
      console.warn(`Cannot generate ${keyType}-RInChIKey: No RInChI provided`);
      return;
    }

    if (!rinchiString.startsWith("RInChI=")) {
      console.warn(
        `Invalid RInChI format for ${keyType}-RInChIKey generation:`,
        rinchiString
      );
      return;
    }

    try {
      const result = await generateRinchiKey(
        rinchiString,
        keyType,
        rinchiVersion
      );

      if (!result) {
        console.error(
          `No result returned from ${keyType}-RInChIKey generation`
        );
        return;
      }

      if (result.return_code !== 0) {
        console.error(`Failed to generate ${keyType}-RInChIKey:`, result.error);
        return;
      }

      // Update the appropriate state based on key type
      switch (keyType) {
        case "Long":
          setLongRinchiKey(result.rinchikey);
          break;
        case "Short":
          setShortRinchiKey(result.rinchikey);
          break;
        case "Web":
          setWebRinchiKey(result.rinchikey);
          break;
        default:
          console.warn(`Unknown RInChIKey type: ${keyType}`);
      }

      return result.rinchikey;
    } catch (err) {
      console.error(`Error generating ${keyType}-RInChIKey:`, err);
      return null;
    }
  };

  // Load RInChI into Ketcher
  const loadRInChI = async () => {
    if (!inputRinchi.trim()) {
      setError("Please enter a RInChI string");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!rinchiModuleLoaded) {
      setError("RInChI module not loaded. Please wait or refresh the page.");
      return;
    }

    if (!inputRinchi.startsWith("RInChI=")) {
      setError('Invalid RInChI string. RInChI should start with "RInChI="');
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Converting RInChI to reaction...");

    try {
      // Call convertRinchiToFileText to generate a RXN file
      const result = await convertRinchiToFileText(
        inputRinchi,
        inputRauxInfo || "",
        "RXN",
        rinchiVersion
      );

      if (result.return_code !== 0) {
        throw new Error(result.error || "Failed to convert RInChI to reaction");
      }

      if (!result.fileText) {
        throw new Error("No reaction data generated from RInChI");
      }

      // Load the RXN file into Ketcher
      await executeKetcherCommand("setMolecule", [result.fileText]);

      // Update state
      setRinchi(inputRinchi);
      setRauxInfo(inputRauxInfo || "");
      setRxnfileContent(result.fileText);
      setLogMessage(
        `Reaction loaded from RInChI successfully${
          result.error ? `\n${result.error}` : ""
        }`
      );

      // Generate RInChIKeys
      if (inputRinchi) {
        await Promise.all([
          generateRinchiKeyAndUpdateState(inputRinchi, "Long"),
          generateRinchiKeyAndUpdateState(inputRinchi, "Short"),
          generateRinchiKeyAndUpdateState(inputRinchi, "Web"),
        ]);
      }
    } catch (err) {
      console.error("Failed to load RInChI:", err);
      setError(`Failed to convert RInChI to reaction: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Convert RInChI to RXN/RD file
  const convertRInChIToFile = async () => {
    if (!inputRinchi.trim()) {
      setError("Please enter a RInChI string");
      return;
    }

    if (!rinchiModuleLoaded) {
      setError("RInChI module not loaded. Please wait or refresh the page.");
      return;
    }

    if (!inputRinchi.startsWith("RInChI=")) {
      setError('Invalid RInChI string. RInChI should start with "RInChI="');
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage(`Converting RInChI to ${outputFormat} file...`);
    setOutputRxnFile("");

    try {
      // Call convertRinchiToFileText to generate a file
      const result = await convertRinchiToFileText(
        inputRinchi,
        inputRauxInfo || "",
        outputFormat,
        rinchiVersion
      );

      if (result.return_code !== 0) {
        throw new Error(
          result.error || `Failed to convert RInChI to ${outputFormat} file`
        );
      }

      if (!result.fileText) {
        throw new Error(`No ${outputFormat} file data generated from RInChI`);
      }

      // Update state
      setOutputRxnFile(result.fileText);
      setLogMessage(
        `${outputFormat} file generated successfully${
          result.error ? `\n${result.error}` : ""
        }`
      );
    } catch (err) {
      console.error("Failed to convert RInChI to file:", err);
      setError(
        `Failed to convert RInChI to ${outputFormat} file: ${err.message}`
      );
    } finally {
      setIsLoading(false);
    }
  };

  // Clear the editor
  const clearEditor = async () => {
    if (!isEditorReady) {
      console.warn("Editor not ready for clearing");
      return;
    }

    try {
      await executeKetcherCommand("setMolecule", [""]);
      setRinchi("");
      setRauxInfo("");
      setLongRinchiKey("");
      setShortRinchiKey("");
      setWebRinchiKey("");
      setRxnfileContent("");
      setLogMessage("Editor cleared");
      console.log("Editor cleared successfully");
    } catch (err) {
      console.error("Failed to clear editor:", err);
      setError("Failed to clear the editor");
    }
  };

  // Function to retry initialization
  const handleRetryInit = () => {
    setIsEditorReady(false);
    setError(null);
    setRetryAttempt((prev) => prev + 1);
  };

  // Examples of common reactions with their RInChI
  const examples = [
    {
      name: "Esterification",
      rinchi:
        "RInChI=1.00.1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)!C2H6O/c1-2-3/h3H,2H2,1H3<>C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3!H2O/h1H2/d+",
      description: "Acetic acid + ethanol → ethyl acetate + water",
    },
    {
      name: "Diels-Alder",
      rinchi:
        "RInChI=1.00.1S/C4H6/c1-3-4-2/h3-4H,1-2H2!C2H2/c1-2/h1-2H<>C6H8/c1-2-4-6-5-3-1/h1-6H/d-",
      description: "1,3-butadiene + ethyne → cyclohexadiene",
    },
    {
      name: "Equilibrium",
      rinchi:
        "RInChI=1.00.1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1<>C6H12O6/c7-1-3-4(9)5(10)6(11)2-8/h3-11H,1-2H2/t3-,4+,5-,6-/m1/s1/d=",
      description: "Glucose ⇌ fructose",
    },
  ];

  return (
    <div className="flex flex-col gap-6 p-4 md:p-6 bg-gradient-to-br from-indigo-50 to-blue-100 dark:from-gray-900 dark:to-slate-900 min-h-screen">
      {/* Copy Modal - For fallback copying */}
      {showCopyModal && (
        <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4 animate-fadeIn">
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-2xl max-w-lg w-full">
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-bold text-gray-900 dark:text-white">
                Copy Text
              </h3>
              <button
                onClick={() => setShowCopyModal(false)}
                className="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200"
              >
                <HiOutlineX className="w-5 h-5" />
              </button>
            </div>
            <p className="mb-4 text-gray-700 dark:text-gray-300">
              Automatic copying failed. Please select and copy this text
              manually:
            </p>
            <div className="mb-4">
              <input
                type="text"
                ref={copyTextRef}
                value={copyModalText}
                readOnly
                className="w-full p-2 border border-gray-300 dark:border-gray-600 rounded font-mono text-sm bg-gray-50 dark:bg-gray-900"
                onClick={(e) => e.target.select()}
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

      {/* Main content area */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
        {/* Left sidebar with controls */}
        <div className="lg:col-span-3 flex flex-col gap-4">
          {/* RInChI Options */}
          <RInChIOptions
            onChange={() => {}} // RInChI doesn't have many options that affect generation
            rinchiVersion={rinchiVersion}
            setRinchiVersion={setRinchiVersion}
            forceEquilibrium={forceEquilibrium}
            setForceEquilibrium={setForceEquilibrium}
            outputFormat={outputFormat}
            setOutputFormat={setOutputFormat}
          />

          {/* Global Editor Controls */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
            <div className="flex items-center mb-4">
              <div className="bg-blue-100 dark:bg-blue-900/50 p-2 rounded-lg mr-3">
                <HiOutlineArrowsExpand className="h-5 w-5 text-blue-700 dark:text-blue-400" />
              </div>
              <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                Editor Controls
              </h2>
            </div>

            <div className="space-y-4">
              <button
                onClick={clearEditor}
                disabled={isLoading || !isEditorReady}
                className={`w-full px-4 py-2 rounded-lg text-sm font-medium transition-colors flex items-center justify-center ${
                  isLoading || !isEditorReady
                    ? "bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed"
                    : "bg-white dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-700 border border-gray-300 dark:border-gray-600"
                }`}
              >
                <HiOutlineRefresh className="h-4 w-4 mr-2" />
                Clear Editor
              </button>

              <button
                onClick={generateRInChI}
                disabled={isLoading || !isEditorReady || !rinchiModuleLoaded}
                className={`w-full px-4 py-2.5 rounded-lg text-sm font-medium transition-all duration-300 flex items-center justify-center ${
                  isLoading || !isEditorReady || !rinchiModuleLoaded
                    ? "bg-gray-400 dark:bg-gray-600 text-white cursor-not-allowed"
                    : "bg-gradient-to-r from-indigo-600 to-blue-700 hover:from-indigo-700 hover:to-blue-800 text-white"
                }`}
              >
                <HiOutlineDocumentText className="h-4 w-4 mr-2" />
                {isLoading ? "Generating..." : "Generate RInChI"}
              </button>
            </div>
          </div>

          {/* Tab-specific controls */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
            <div className="flex items-center mb-4">
              <div className="bg-indigo-100 dark:bg-indigo-900/50 p-2 rounded-lg mr-3">
                <HiOutlinePencil className="h-5 w-5 text-indigo-700 dark:text-indigo-400" />
              </div>
              <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                Reaction Input
              </h2>
            </div>

            {/* Input Type Tabs */}
            <div className="mb-4">
              <div className="flex flex-wrap border-b border-gray-200 dark:border-gray-700">
                <button
                  className={`px-4 py-2 text-sm font-medium ${
                    activeInputType === "reaction"
                      ? "text-indigo-600 dark:text-indigo-400 border-b-2 border-indigo-600 dark:border-indigo-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("reaction")}
                >
                  Draw
                </button>
                <button
                  className={`px-4 py-2 text-sm font-medium ${
                    activeInputType === "rxnfile"
                      ? "text-indigo-600 dark:text-indigo-400 border-b-2 border-indigo-600 dark:border-indigo-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("rxnfile")}
                >
                  RXN/RD File
                </button>
                <button
                  className={`px-4 py-2 text-sm font-medium ${
                    activeInputType === "rinchi"
                      ? "text-indigo-600 dark:text-indigo-400 border-b-2 border-indigo-600 dark:border-indigo-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("rinchi")}
                >
                  RInChI
                </button>
                <button
                  className={`px-4 py-2 text-sm font-medium ${
                    activeInputType === "rinchi-to-file"
                      ? "text-indigo-600 dark:text-indigo-400 border-b-2 border-indigo-600 dark:border-indigo-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("rinchi-to-file")}
                >
                  RInChI → File
                </button>
              </div>
            </div>

            {/* RXN/RD File Input Option */}
            {activeInputType === "rxnfile" && (
              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="rxnfile-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter RXN/RD File Content
                  </label>
                  <textarea
                    id="rxnfile-input"
                    value={rxnfileContent}
                    onChange={(e) => setRxnfileContent(e.target.value)}
                    placeholder="Paste RXN/RD file content here..."
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 transition-all duration-200"
                    rows={6}
                  />
                </div>

                <div className="flex flex-col space-y-3">
                  <div className="relative">
                    <input
                      type="file"
                      id="rxnfile-upload"
                      accept=".rxn,.rd"
                      ref={fileInputRef}
                      onChange={handleRxnFileUpload}
                      className="hidden"
                    />
                    <label
                      htmlFor="rxnfile-upload"
                      className="cursor-pointer flex items-center justify-center w-full px-4 py-2.5 bg-gray-100 dark:bg-gray-700 text-gray-800 dark:text-gray-200 rounded-lg hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors duration-300"
                    >
                      <HiOutlineUpload className="w-5 h-5 mr-2" />
                      Upload RXN/RD File
                    </label>
                  </div>

                  <button
                    onClick={processRxnFile}
                    disabled={
                      isLoading || !rxnfileContent.trim() || !rinchiModuleLoaded
                    }
                    className={`relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 ${
                      isLoading || !rxnfileContent.trim() || !rinchiModuleLoaded
                        ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                        : "bg-gradient-to-r from-indigo-600 to-blue-700 hover:from-indigo-700 hover:to-blue-800"
                    }`}
                  >
                    <HiOutlineDocumentText className="mr-2 h-5 w-5" />
                    {isLoading ? "Processing..." : "Process RXN/RD File"}
                  </button>
                </div>
              </div>
            )}

            {/* RInChI Input Option for conversion to structure */}
            {activeInputType === "rinchi" && (
              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="rinchi-input-structure"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter RInChI to Convert
                  </label>
                  <textarea
                    id="rinchi-input-structure"
                    value={inputRinchi}
                    onChange={(e) => setInputRinchi(e.target.value)}
                    placeholder="e.g., RInChI=1.00.1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)!C2H6O/c1-2-3/h3H,2H2,1H3<>C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3!H2O/h1H2/d+"
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 transition-all duration-200"
                    rows={4}
                  />
                </div>

                <div>
                  <label
                    htmlFor="rauxinfo-input-structure"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    RAuxInfo (Optional)
                  </label>
                  <textarea
                    id="rauxinfo-input-structure"
                    value={inputRauxInfo}
                    onChange={(e) => setInputRauxInfo(e.target.value)}
                    placeholder="e.g., RAuxInfo=1.00.1/1/N:1,2,3,4/E:(3,4)/rA:4nCCOO/rB:s1;s2;d2;/rC:-3.8549,-.5552,0;..."
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 transition-all duration-200"
                    rows={3}
                  />
                </div>

                <button
                  onClick={loadRInChI}
                  disabled={
                    isLoading ||
                    !isEditorReady ||
                    !inputRinchi.trim() ||
                    !rinchiModuleLoaded
                  }
                  className={`w-full relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 ${
                    isLoading ||
                    !isEditorReady ||
                    !inputRinchi.trim() ||
                    !rinchiModuleLoaded
                      ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                      : "bg-gradient-to-r from-indigo-600 to-blue-700 hover:from-indigo-700 hover:to-blue-800"
                  }`}
                >
                  {isLoading
                    ? "Converting..."
                    : !isEditorReady
                    ? "Initializing..."
                    : !rinchiModuleLoaded
                    ? "Loading RInChI Module..."
                    : "Convert to Reaction"}
                </button>

                {/* Quick Examples */}
                <div className="mt-2">
                  <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                    Quick Examples:
                  </p>
                  <div className="flex flex-wrap gap-2">
                    {examples.map((example, index) => (
                      <button
                        key={index}
                        onClick={() => setInputRinchi(example.rinchi)}
                        className="inline-flex items-center px-3 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded-full text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-indigo-500 dark:focus:ring-offset-gray-800 transition-all duration-200"
                        title={`${example.name}: ${example.description}`}
                      >
                        {example.name}
                      </button>
                    ))}
                  </div>
                </div>
              </div>
            )}

            {/* RInChI Input Option for conversion to RXN/RD file */}
            {activeInputType === "rinchi-to-file" && (
              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="rinchi-input-file"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter RInChI to Convert
                  </label>
                  <textarea
                    id="rinchi-input-file"
                    value={inputRinchi}
                    onChange={(e) => setInputRinchi(e.target.value)}
                    placeholder="e.g., RInChI=1.00.1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)!C2H6O/c1-2-3/h3H,2H2,1H3<>C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3!H2O/h1H2/d+"
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 transition-all duration-200"
                    rows={4}
                  />
                </div>

                <div>
                  <label
                    htmlFor="rauxinfo-input-file"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    RAuxInfo (Optional)
                  </label>
                  <textarea
                    id="rauxinfo-input-file"
                    value={inputRauxInfo}
                    onChange={(e) => setInputRauxInfo(e.target.value)}
                    placeholder="e.g., RAuxInfo=1.00.1/1/N:1,2,3,4/E:(3,4)/rA:4nCCOO/rB:s1;s2;d2;/rC:-3.8549,-.5552,0;..."
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 transition-all duration-200"
                    rows={3}
                  />
                </div>

                <button
                  onClick={convertRInChIToFile}
                  disabled={
                    isLoading || !inputRinchi.trim() || !rinchiModuleLoaded
                  }
                  className={`w-full relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-indigo-500 ${
                    isLoading || !inputRinchi.trim() || !rinchiModuleLoaded
                      ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                      : "bg-gradient-to-r from-indigo-600 to-blue-700 hover:from-indigo-700 hover:to-blue-800"
                  }`}
                >
                  <HiOutlineSwitchHorizontal className="mr-2 h-5 w-5" />
                  {isLoading
                    ? "Converting..."
                    : !rinchiModuleLoaded
                    ? "Loading RInChI Module..."
                    : `Convert to ${outputFormat} File`}
                </button>

                {/* Quick Examples */}
                <div className="mt-2">
                  <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                    Quick Examples:
                  </p>
                  <div className="flex flex-wrap gap-2">
                    {examples.map((example, index) => (
                      <button
                        key={index}
                        onClick={() => setInputRinchi(example.rinchi)}
                        className="inline-flex items-center px-3 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded-full text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-indigo-500 dark:focus:ring-offset-gray-800 transition-all duration-200"
                        title={`${example.name}: ${example.description}`}
                      >
                        {example.name}
                      </button>
                    ))}
                  </div>
                </div>
              </div>
            )}

            {/* Drawing Instructions */}
            {activeInputType === "reaction" && (
              <div className="text-sm text-gray-700 dark:text-gray-300">
                <p>
                  Use the Ketcher editor to draw your chemical reaction, then
                  click "Generate RInChI" to convert it.
                </p>
                <ul className="list-disc list-inside mt-2 space-y-1 text-gray-600 dark:text-gray-400">
                  <li>Draw reactants on the left side</li>
                  <li>Add reaction arrow from the toolbar</li>
                  <li>Draw products on the right side</li>
                  <li>Check "Force Equilibrium" for equilibrium reactions</li>
                </ul>
              </div>
            )}
          </div>

          {/* Status and Information */}
          <div
            className={`px-6 py-4 rounded-xl flex items-center justify-between ${
              isEditorReady && rinchiModuleLoaded
                ? "bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800/40"
                : "bg-yellow-50 dark:bg-yellow-900/20 border border-yellow-200 dark:border-yellow-800/40"
            }`}
          >
            <div className="flex flex-col space-y-2">
              <div className="flex items-center">
                <div
                  className={`h-3 w-3 rounded-full mr-3 ${
                    isEditorReady
                      ? "bg-green-500 animate-pulse"
                      : "bg-yellow-500 animate-pulse"
                  }`}
                ></div>
                <span
                  className={`text-sm font-medium ${
                    isEditorReady
                      ? "text-green-800 dark:text-green-300"
                      : "text-yellow-800 dark:text-yellow-300"
                  }`}
                >
                  {isEditorReady ? "Editor Ready" : "Initializing Editor..."}
                </span>
              </div>

              <div className="flex items-center">
                <div
                  className={`h-3 w-3 rounded-full mr-3 ${
                    rinchiModuleLoaded
                      ? "bg-green-500 animate-pulse"
                      : "bg-yellow-500 animate-pulse"
                  }`}
                ></div>
                <span
                  className={`text-sm font-medium ${
                    rinchiModuleLoaded
                      ? "text-green-800 dark:text-green-300"
                      : "text-yellow-800 dark:text-yellow-300"
                  }`}
                >
                  {rinchiModuleLoaded
                    ? "RInChI Module Ready"
                    : "Loading RInChI Module..."}
                </span>
              </div>
            </div>

            <div className="flex gap-3">
              {(!isEditorReady || !rinchiModuleLoaded) && (
                <button
                  onClick={handleRetryInit}
                  className="inline-flex items-center px-3 py-1.5 rounded-lg text-sm font-medium transition-colors bg-yellow-100 dark:bg-yellow-900/30 text-yellow-700 dark:text-yellow-300 hover:bg-yellow-200 dark:hover:bg-yellow-900/50 border border-yellow-300 dark:border-yellow-700/50"
                >
                  <HiOutlineRefresh className="h-4 w-4 mr-1.5" />
                  Retry
                </button>
              )}
            </div>
          </div>

          {/* Information Box */}
          <div className="bg-gradient-to-br from-blue-50 to-indigo-50 dark:from-blue-900/30 dark:to-indigo-900/30 border border-blue-200 dark:border-blue-800/50 rounded-xl p-5 text-sm shadow-lg">
            <h4 className="font-bold text-blue-800 dark:text-blue-300 mb-3 flex items-center">
              <HiOutlineInformationCircle className="h-5 w-5 mr-2 text-blue-500 dark:text-blue-400" />
              About RInChI
            </h4>
            <div className="space-y-3 text-gray-700 dark:text-gray-300">
              <p>
                The IUPAC RInChI (Reaction InChI) is a textual identifier for
                chemical reactions, extending the InChI concept to represent
                chemical transformations.
              </p>
              <div>
                <h5 className="font-medium mb-1 text-gray-800 dark:text-gray-200">
                  Key Features:
                </h5>
                <ul className="list-disc list-inside space-y-1 pl-1 text-gray-600 dark:text-gray-400">
                  <li>
                    Represents chemical reactions and their directionality
                  </li>
                  <li>Combines InChI strings of reactants and products</li>
                  <li>Supports multiple reactants, agents, and products</li>
                  <li>
                    Includes directionality indicators (forward, reverse,
                    equilibrium)
                  </li>
                  <li>RInChIKeys allow for easy reaction searching</li>
                </ul>
              </div>
            </div>
          </div>
        </div>

        {/* Main Editor Area */}
        <div className="lg:col-span-9 flex flex-col gap-4">
          {/* Error Display */}
          {error && (
            <div
              className="p-4 rounded-xl bg-red-50 dark:bg-red-900/30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700/50 flex items-start shadow-lg animate-fadeIn"
              role="alert"
            >
              <HiOutlineExclamationCircle
                className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
                aria-hidden="true"
              />
              <p>{error}</p>
            </div>
          )}

          {/* Ketcher Editor Container (only show when needed) */}
          {activeInputType !== "rinchi-to-file" && (
            <div
              className="bg-white dark:bg-gray-800 rounded-xl shadow-xl border border-gray-200 dark:border-gray-700 overflow-hidden"
              style={{ height: "900px" }}
            >
              <iframe
                ref={ketcherFrame}
                src={`${process.env.PUBLIC_URL}/standalone/index.html`}
                title="Ketcher Editor"
                className="w-full h-full"
                frameBorder="0"
                onLoad={(e) => {
                  // Don't prevent default - just let communication script be injected
                  console.log("Ketcher iframe loaded");
                }}
              />
            </div>
          )}

          {/* Results */}
          <div className="bg-white dark:bg-gray-800 p-5 rounded-xl border border-gray-200 dark:border-gray-700 shadow-lg">
            <h3 className="text-lg font-bold text-gray-800 dark:text-white mb-3 flex items-center">
              <HiOutlineDocumentText className="h-5 w-5 mr-2 text-indigo-500 dark:text-indigo-400" />
              Results
            </h3>

            {/* RInChI */}
            <ResultBlock
              title="RInChI"
              value={rinchi}
              onCopy={() => copyToClipboard(rinchi, "rinchi")}
              copyState={copySuccess}
              icon={
                <HiOutlineCode className="h-4 w-4 text-indigo-600 dark:text-indigo-500" />
              }
            />

            {/* RAuxInfo */}
            <ResultBlock
              title="RAuxInfo"
              value={rauxInfo}
              onCopy={() => copyToClipboard(rauxInfo, "rauxinfo")}
              copyState={copySuccess}
              icon={
                <HiOutlineDocumentText className="h-4 w-4 text-blue-600 dark:text-blue-500" />
              }
            />

            {/* RInChIKeys */}
            <div className="mt-4">
              <h4 className="text-sm font-bold text-gray-700 dark:text-gray-300 mb-2 flex items-center">
                <HiOutlineInformationCircle className="h-4 w-4 mr-2 text-indigo-600 dark:text-indigo-500" />
                RInChIKeys
              </h4>

              <div className="space-y-2">
                {/* Long-RInChIKey */}
                <ResultBlock
                  title="Long-RInChIKey"
                  value={longRinchiKey}
                  onCopy={() => copyToClipboard(longRinchiKey, "longkey")}
                  copyState={copySuccess}
                  collapsible={true}
                />

                {/* Short-RInChIKey */}
                <ResultBlock
                  title="Short-RInChIKey"
                  value={shortRinchiKey}
                  onCopy={() => copyToClipboard(shortRinchiKey, "shortkey")}
                  copyState={copySuccess}
                  collapsible={true}
                />

                {/* Web-RInChIKey */}
                <ResultBlock
                  title="Web-RInChIKey"
                  value={webRinchiKey}
                  onCopy={() => copyToClipboard(webRinchiKey, "webkey")}
                  copyState={copySuccess}
                  collapsible={true}
                />
              </div>
            </div>

            {/* Output RXN/RD File */}
            {activeInputType === "rinchi-to-file" && outputRxnFile && (
              <ResultBlock
                title={`${outputFormat} File Output`}
                value={outputRxnFile}
                onCopy={() => copyToClipboard(outputRxnFile, "rxnfile")}
                copyState={copySuccess}
                icon={
                  <HiOutlineDocumentText className="h-4 w-4 text-blue-600 dark:text-blue-500" />
                }
              />
            )}

            {/* Log messages */}
            <div className="mt-4 p-3 bg-gray-50 dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700">
              <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                <HiOutlineInformationCircle className="h-4 w-4 mr-1 text-gray-600 dark:text-gray-400" />
                Log
              </h4>
              <pre className="text-xs text-gray-600 dark:text-gray-400 font-mono whitespace-pre-wrap max-h-48 overflow-auto">
                {logMessage || "No log messages yet."}
              </pre>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default RInChIView;
