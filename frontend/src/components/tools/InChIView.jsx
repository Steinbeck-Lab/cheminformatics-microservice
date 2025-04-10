import React, { useState, useEffect, useRef } from "react";
import {
  HiOutlineClipboardCopy,
  HiOutlineCheck,
  HiOutlineExclamationCircle,
  HiOutlineBeaker,
  HiOutlineInformationCircle,
  HiOutlineRefresh,
  HiOutlinePencil,
  HiOutlineX,
  HiOutlineCode,
  HiOutlineDocumentText,
  HiOutlineDocumentDownload,
} from "react-icons/hi";

// Import utility functions
import {
  INCHI_VERSIONS,
  useMockImplementation,
  convertMolfileToInchi,
  generateInchiKey,
  convertInchiToMolfile,
  mockConvertMolfileToInchi,
  mockGenerateInchiKey,
  mockConvertInchiToMolfile,
} from "../../utils/inchiUtils.js";

// InChI options component
const InChIOptions = ({ onChange, inchiVersion, setInchiVersion }) => {
  const [mobileH, setMobileH] = useState(true);
  const [includeStereo, setIncludeStereo] = useState(true);
  const [stereoMode, setStereoMode] = useState("SAbs");
  const [recMet, setRecMet] = useState(false);
  const [ket, setKet] = useState(false);
  const [t15, setT15] = useState(false);
  const [polymers, setPolymers] = useState(false);
  const [NPZz, setNPZz] = useState(false);
  const [orgMet, setOrgMet] = useState(false);

  // Generate and propagate options string when settings change
  useEffect(() => {
    const options = [];

    // Add options based on current settings
    if (!mobileH) options.push("FixedH");
    if (!includeStereo) options.push("SNon");

    if (includeStereo) {
      if (stereoMode === "SRel") options.push("SRel");
      if (stereoMode === "SRac") options.push("SRac");
      if (stereoMode === "SUCF") options.push("SUCF");
    }

    if (recMet) options.push("RecMet");
    if (ket) options.push("KET");
    if (t15) options.push("15T");
    if (polymers) options.push("Polymers");
    if (NPZz) options.push("NPZz");

    // orgMet is only available in the 1.07.3-orgmet version
    if (orgMet && inchiVersion === "1.07.3-orgmet") options.push("OrgMet");

    // Format the options string as required by the API
    const optionsString = options.map((opt) => `-${opt}`).join(" ");
    onChange(optionsString);
  }, [
    mobileH,
    includeStereo,
    stereoMode,
    recMet,
    ket,
    t15,
    polymers,
    NPZz,
    orgMet,
    inchiVersion,
    onChange,
  ]);

  // Update OrgMet option when InChI version changes
  useEffect(() => {
    if (inchiVersion !== "1.07.3-orgmet" && orgMet) {
      setOrgMet(false);
    }
  }, [inchiVersion, orgMet]);

  return (
    <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
      <div className="flex items-center mb-4">
        <div className="bg-purple-100 dark:bg-purple-900/50 p-2 rounded-lg mr-3">
          <HiOutlineCode className="h-5 w-5 text-purple-700 dark:text-purple-400" />
        </div>
        <h2 className="text-lg font-bold text-gray-800 dark:text-white">
          InChI Options
        </h2>
      </div>

      <div className="mb-4">
        <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
          InChI Version
        </label>
        <select
          className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-green-500 focus:border-green-500"
          value={inchiVersion}
          onChange={(e) => setInchiVersion(e.target.value)}
        >
          {Object.entries(INCHI_VERSIONS).map(([value, version]) => (
            <option key={value} value={value}>
              {version.label}
            </option>
          ))}
        </select>
      </div>

      <div className="space-y-4">
        <div className="flex items-center">
          <input
            id="mobileH"
            type="checkbox"
            checked={mobileH}
            onChange={(e) => setMobileH(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="mobileH"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Mobile H Perception
          </label>
        </div>

        <div>
          <div className="flex items-center">
            <input
              id="includeStereo"
              type="checkbox"
              checked={includeStereo}
              onChange={(e) => setIncludeStereo(e.target.checked)}
              className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
            />
            <label
              htmlFor="includeStereo"
              className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
            >
              Include Stereo
            </label>
          </div>

          {includeStereo && (
            <div className="ml-6 mt-2 grid grid-cols-2 gap-2">
              <div className="flex items-center">
                <input
                  id="stereo-abs"
                  type="radio"
                  name="stereoMode"
                  value="SAbs"
                  checked={stereoMode === "SAbs"}
                  onChange={() => setStereoMode("SAbs")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-abs"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Absolute
                </label>
              </div>

              <div className="flex items-center">
                <input
                  id="stereo-rel"
                  type="radio"
                  name="stereoMode"
                  value="SRel"
                  checked={stereoMode === "SRel"}
                  onChange={() => setStereoMode("SRel")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-rel"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Relative
                </label>
              </div>

              <div className="flex items-center">
                <input
                  id="stereo-rac"
                  type="radio"
                  name="stereoMode"
                  value="SRac"
                  checked={stereoMode === "SRac"}
                  onChange={() => setStereoMode("SRac")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-rac"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Racemic
                </label>
              </div>

              <div className="flex items-center">
                <input
                  id="stereo-ucf"
                  type="radio"
                  name="stereoMode"
                  value="SUCF"
                  checked={stereoMode === "SUCF"}
                  onChange={() => setStereoMode("SUCF")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-ucf"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  From chiral flag
                </label>
              </div>
            </div>
          )}
        </div>

        <div className="flex items-center">
          <input
            id="recmet"
            type="checkbox"
            checked={recMet}
            onChange={(e) => setRecMet(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="recmet"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Include Bonds to Metal
          </label>
        </div>

        <div className="flex items-center">
          <input
            id="ket"
            type="checkbox"
            checked={ket}
            onChange={(e) => setKet(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="ket"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Keto-enol Tautomerism
          </label>
        </div>

        <div className="flex items-center">
          <input
            id="15t"
            type="checkbox"
            checked={t15}
            onChange={(e) => setT15(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="15t"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            1,5-tautomerism
          </label>
        </div>

        <div className="flex items-center">
          <input
            id="polymers"
            type="checkbox"
            checked={polymers}
            onChange={(e) => setPolymers(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="polymers"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Treat Polymers
          </label>
        </div>

        <div className="flex items-center">
          <input
            id="npzz"
            type="checkbox"
            checked={NPZz}
            onChange={(e) => setNPZz(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="npzz"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Allow non-polymer Zz pseudoatoms
          </label>
        </div>

        {inchiVersion === "1.07.3-orgmet" && (
          <div className="flex items-center">
            <input
              id="orgmet"
              type="checkbox"
              checked={orgMet}
              onChange={(e) => setOrgMet(e.target.checked)}
              className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
            />
            <label
              htmlFor="orgmet"
              className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
            >
              Molecular inorganics
            </label>
          </div>
        )}
      </div>
    </div>
  );
};

// Result block component for displaying InChI, InChIKey, and AuxInfo
const ResultBlock = ({ title, value, onCopy, copyState }) => {
  return (
    <div className="mt-4 bg-white dark:bg-gray-800 p-4 rounded-lg border border-gray-200 dark:border-gray-700 shadow">
      <div className="flex justify-between items-center mb-2">
        <h3 className="text-sm font-medium text-gray-700 dark:text-gray-300">
          {title}
        </h3>
        {value && (
          <button
            onClick={onCopy}
            className="inline-flex items-center px-2 py-1 text-xs text-green-700 dark:text-green-300 bg-green-100 dark:bg-green-900/30 rounded hover:bg-green-200 dark:hover:bg-green-800/50"
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
      <div className="font-mono text-sm overflow-x-auto bg-gray-50 dark:bg-gray-900 p-3 rounded border border-gray-200 dark:border-gray-700 break-all max-h-24 overflow-y-auto">
        {value || (
          <span className="text-gray-400 dark:text-gray-500">
            No data generated yet
          </span>
        )}
      </div>
    </div>
  );
};

const InChIView = () => {
  // State variables
  const [inchi, setInchi] = useState("");
  const [auxInfo, setAuxInfo] = useState("");
  const [inchiKey, setInchiKey] = useState("");
  const [inputInchi, setInputInchi] = useState("");
  const [molfileContent, setMolfileContent] = useState("");
  const [logMessage, setLogMessage] = useState("");
  const [options, setOptions] = useState("");
  const [inchiVersion, setInchiVersion] = useState(
    Object.keys(INCHI_VERSIONS).find((key) => INCHI_VERSIONS[key].default) ||
      "1.07.3"
  );
  const [activeTab, setActiveTab] = useState("draw");

  const [isLoading, setIsLoading] = useState(false);
  const [isEditorReady, setIsEditorReady] = useState(false);
  const [error, setError] = useState(null);
  const [copySuccess, setCopySuccess] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [copyModalText, setCopyModalText] = useState("");
  const [inchiModuleLoaded, setInchiModuleLoaded] = useState(false);
  const [retryAttempt, setRetryAttempt] = useState(0);

  // Refs
  const ketcherFrame = useRef(null);
  const copyTextRef = useRef(null);
  const messageHandlers = useRef({});
  const messageId = useRef(0);

  // Reset error when needed
  useEffect(() => {
    if (error) setError(null);
  }, [inputInchi, inchiVersion, options, activeTab, error]); // Added error to dependency array

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
  const copyToClipboard = async (text = null, type = "inchi") => {
    let textToCopy;

    switch (type) {
      case "inchi":
        textToCopy = text || inchi;
        break;
      case "key":
        textToCopy = text || inchiKey;
        break;
      case "auxinfo":
        textToCopy = text || auxInfo;
        break;
      default:
        textToCopy = text || inchi;
    }

    if (!textToCopy) {
      setError(
        `No ${type.toUpperCase()} to copy. Generate ${type.toUpperCase()} first.`
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

  // Check InChI module status
  useEffect(() => {
    const checkModuleStatus = async () => {
      try {
        // This just tries to load the module but doesn't use it yet
        // We set inchiModuleLoaded to true if useMockImplementation is true
        // or if we successfully load the real module
        setInchiModuleLoaded(useMockImplementation ? true : false);

        if (!useMockImplementation) {
          // Try to load the real module
          // This will set a "loading" state but won't actually use the results
          // Our utility functions will handle the actual module interaction
          await import("../../utils/inchiUtils")
            .then(() => {
              setInchiModuleLoaded(true);
            })
            .catch((err) => {
              console.error("Failed to load InChI utilities:", err);
              setError(
                "Failed to load InChI utilities. Mock implementation will be used."
              );
              setInchiModuleLoaded(true); // We'll fall back to mocks
            });
        }
      } catch (err) {
        console.error("Error checking InChI module status:", err);
        setInchiModuleLoaded(useMockImplementation);
      }
    };

    checkModuleStatus();
  }, [inchiVersion]);

  // Load InChI into Ketcher
  const loadInChI = async () => {
    if (!inputInchi.trim()) {
      setError("Please enter an InChI string");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!inputInchi.startsWith("InChI=")) {
      setError('Invalid InChI string. InChI should start with "InChI="');
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("");

    try {
      // Call the appropriate function based on whether we're using mocks
      const result = useMockImplementation
        ? await mockConvertInchiToMolfile(inputInchi)
        : await convertInchiToMolfile(inputInchi, "", inchiVersion);

      if (result.return_code === -1) {
        throw new Error(
          result.message || "Failed to convert InChI to structure"
        );
      }

      if (result.molfile) {
        // Load the molfile into Ketcher
        await executeKetcherCommand("setMolecule", [result.molfile]);
        console.log("InChI loaded successfully");

        // Update state
        setInchi(inputInchi);
        setLogMessage(result.message || "Structure loaded successfully");

        // Generate InChIKey
        generateInChIKeyFromInChI(inputInchi);
      } else {
        throw new Error("No structure data generated from InChI");
      }
    } catch (err) {
      console.error("Failed to load InChI:", err);
      setError(`Failed to convert InChI to structure: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Generate InChI from Ketcher
  const generateInChI = async () => {
    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("");
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      // Get molfile from Ketcher
      const molfile = await executeKetcherCommand("getMolfile");

      if (!molfile || molfile.trim() === "") {
        setError("No structure drawn. Please draw a molecule first.");
        setIsLoading(false);
        return;
      }

      // Check if it's a reaction (not supported by InChI)
      if (molfile.includes("$RXN")) {
        setError("Reactions are not supported. Please draw a single molecule.");
        setIsLoading(false);
        return;
      }

      // Save the molfile for reference
      setMolfileContent(molfile);

      // Convert molfile to InChI using appropriate function
      const result = useMockImplementation
        ? await mockConvertMolfileToInchi(molfile, options)
        : await convertMolfileToInchi(molfile, options, inchiVersion);

      if (result.return_code === -1) {
        throw new Error(result.message || "Failed to generate InChI");
      }

      // Update state with results
      setInchi(result.inchi);
      setAuxInfo(result.auxinfo);
      setLogMessage(result.message || "InChI generated successfully");

      // Generate InChIKey
      if (result.inchi) {
        generateInChIKeyFromInChI(result.inchi);
      }
    } catch (err) {
      console.error("Failed to generate InChI:", err);
      setError(`Failed to generate InChI: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  // Generate InChIKey from InChI
  const generateInChIKeyFromInChI = async (inchiString) => {
    if (!inchiString) return;

    try {
      const result = useMockImplementation
        ? await mockGenerateInchiKey(inchiString)
        : await generateInchiKey(inchiString, inchiVersion);

      if (result.return_code === -1) {
        console.error("Failed to generate InChIKey:", result.message);
        return;
      }

      setInchiKey(result.inchikey);
    } catch (err) {
      console.error("Error generating InChIKey:", err);
    }
  };

  // Convert molfile string to InChI
  const convertMolfileToInChI = async () => {
    if (!molfileContent.trim()) {
      setError("Please enter or upload a Molfile");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("");
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      // Convert molfile to InChI using appropriate function
      const result = useMockImplementation
        ? await mockConvertMolfileToInchi(molfileContent, options)
        : await convertMolfileToInchi(molfileContent, options, inchiVersion);

      if (result.return_code === -1) {
        throw new Error(
          result.message || "Failed to generate InChI from molfile"
        );
      }

      // Update state with results
      setInchi(result.inchi);
      setAuxInfo(result.auxinfo);
      setLogMessage(result.message || "InChI generated successfully");

      // Generate InChIKey
      if (result.inchi) {
        generateInChIKeyFromInChI(result.inchi);
      }
    } catch (err) {
      console.error("Failed to convert molfile to InChI:", err);
      setError(`Failed to generate InChI: ${err.message}`);
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
      setInchi("");
      setAuxInfo("");
      setInchiKey("");
      setMolfileContent("");
      setLogMessage("");
      console.log("Editor cleared successfully");
    } catch (err) {
      console.error("Failed to clear editor:", err);
      setError("Failed to clear the editor");
    }
  };

  // Handle file upload for molfile
  const handleMolfileUpload = (event) => {
    const file = event.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      setMolfileContent(e.target.result);
    };
    reader.readAsText(file);
  };

  // Function to retry initialization
  const handleRetryInit = () => {
    setIsEditorReady(false);
    setError(null);
    setRetryAttempt((prev) => prev + 1);
  };

  // Examples of common molecules with their InChI
  const examples = [
    {
      name: "Ethanol",
      value: "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
      description: "Alcohol",
    },
    {
      name: "Aspirin",
      value:
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
      description: "Pain reliever",
    },
    {
      name: "Caffeine",
      value:
        "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
      description: "Stimulant",
    },
    {
      name: "Ibuprofen",
      value:
        "InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)",
      description: "Anti-inflammatory",
    },
  ];

  return (
    <div className="flex flex-col gap-6 p-4 md:p-6 bg-gradient-to-br from-green-50 to-emerald-100 dark:from-gray-900 dark:to-slate-900 min-h-screen">
      {/* Header with animated background */}
      <div className="bg-gradient-to-r from-green-600 to-emerald-700 dark:from-green-700 dark:to-emerald-900 rounded-xl shadow-xl overflow-hidden relative">
        {/* Animated background elements */}
        <div className="absolute inset-0 overflow-hidden opacity-10">
          <div className="absolute left-0 top-0 w-40 h-40 rounded-full bg-white transform -translate-x-1/2 -translate-y-1/2"></div>
          <div className="absolute right-20 bottom-0 w-64 h-64 rounded-full bg-white transform translate-x-1/2 translate-y-1/2"></div>
          <div className="absolute left-1/3 top-1/2 w-32 h-32 rounded-full bg-white transform -translate-y-1/2"></div>
        </div>

        <div className="relative p-8 md:p-10 flex flex-col md:flex-row items-center justify-between gap-6">
          {/* Title and description */}
          <div className="text-white space-y-2 text-center md:text-left">
            <h1 className="text-3xl md:text-4xl font-bold">
              InChI Structure Editor
            </h1>
            <p className="text-green-100 dark:text-green-200 opacity-90 max-w-xl">
              Draw, edit, and convert chemical structures to InChI notation
            </p>
          </div>

          {/* Decorative icon */}
          <div className="hidden md:flex items-center justify-center bg-white/10 backdrop-blur-sm p-6 rounded-full w-24 h-24 border border-white/20 shadow-lg">
            <HiOutlineBeaker className="w-12 h-12 text-white" />
          </div>
        </div>
      </div>

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

      {/* Navigation tabs */}
      <div className="bg-white dark:bg-gray-800 rounded-xl shadow-lg overflow-hidden">
        <div className="flex border-b border-gray-200 dark:border-gray-700">
          <button
            className={`px-4 py-3 text-sm font-medium flex-1 ${
              activeTab === "draw"
                ? "text-green-600 dark:text-green-400 border-b-2 border-green-600 dark:border-green-400"
                : "text-gray-600 dark:text-gray-400 hover:text-gray-800 hover:bg-gray-50 dark:hover:text-gray-300 dark:hover:bg-gray-700"
            }`}
            onClick={() => setActiveTab("draw")}
          >
            Draw Structure
          </button>
          <button
            className={`px-4 py-3 text-sm font-medium flex-1 ${
              activeTab === "molfile"
                ? "text-green-600 dark:text-green-400 border-b-2 border-green-600 dark:border-green-400"
                : "text-gray-600 dark:text-gray-400 hover:text-gray-800 hover:bg-gray-50 dark:hover:text-gray-300 dark:hover:bg-gray-700"
            }`}
            onClick={() => setActiveTab("molfile")}
          >
            Convert Molfile
          </button>
          <button
            className={`px-4 py-3 text-sm font-medium flex-1 ${
              activeTab === "inchi"
                ? "text-green-600 dark:text-green-400 border-b-2 border-green-600 dark:border-green-400"
                : "text-gray-600 dark:text-gray-400 hover:text-gray-800 hover:bg-gray-50 dark:hover:text-gray-300 dark:hover:bg-gray-700"
            }`}
            onClick={() => setActiveTab("inchi")}
          >
            Convert InChI
          </button>
        </div>
      </div>

      {/* Main content area */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
        {/* Left sidebar with controls */}
        <div className="lg:col-span-3 flex flex-col gap-4">
          {/* InChI Options */}
          <InChIOptions
            onChange={setOptions}
            inchiVersion={inchiVersion}
            setInchiVersion={setInchiVersion}
          />

          {/* Tab-specific controls */}
          {activeTab === "inchi" && (
            <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
              <div className="flex items-center mb-4">
                <div className="bg-indigo-100 dark:bg-indigo-900/50 p-2 rounded-lg mr-3">
                  <HiOutlinePencil className="h-5 w-5 text-indigo-700 dark:text-indigo-400" />
                </div>
                <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                  InChI Input
                </h2>
              </div>

              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="inchi-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter InChI to Convert
                  </label>
                  <textarea
                    id="inchi-input"
                    value={inputInchi}
                    onChange={(e) => setInputInchi(e.target.value)}
                    placeholder="e.g., InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-green-500 focus:border-green-500 transition-all duration-200"
                    rows={4}
                  />
                </div>

                <button
                  onClick={loadInChI}
                  disabled={
                    isLoading ||
                    !isEditorReady ||
                    !inputInchi.trim() ||
                    !inchiModuleLoaded
                  }
                  className={`w-full relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-green-500 ${
                    isLoading ||
                    !isEditorReady ||
                    !inputInchi.trim() ||
                    !inchiModuleLoaded
                      ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                      : "bg-gradient-to-r from-green-600 to-emerald-700 hover:from-green-700 hover:to-emerald-800"
                  }`}
                >
                  {isLoading
                    ? "Loading..."
                    : !isEditorReady
                    ? "Initializing..."
                    : !inchiModuleLoaded
                    ? "Loading InChI Module..."
                    : "Convert to Structure"}
                </button>
              </div>

              {/* Quick Examples */}
              <div className="mt-6">
                <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                  Quick Examples:
                </p>
                <div className="flex flex-wrap gap-2">
                  {examples.map((example, index) => (
                    <button
                      key={index}
                      onClick={() => setInputInchi(example.value)}
                      className="inline-flex items-center px-3 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded-full text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-green-500 dark:focus:ring-offset-gray-800 transition-all duration-200"
                      title={`${example.name}: ${example.description}`}
                    >
                      {example.name}
                    </button>
                  ))}
                </div>
              </div>
            </div>
          )}

          {activeTab === "molfile" && (
            <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
              <div className="flex items-center mb-4">
                <div className="bg-blue-100 dark:bg-blue-900/50 p-2 rounded-lg mr-3">
                  <HiOutlineDocumentText className="h-5 w-5 text-blue-700 dark:text-blue-400" />
                </div>
                <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                  Molfile Input
                </h2>
              </div>

              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="molfile-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter Molfile Content
                  </label>
                  <textarea
                    id="molfile-input"
                    value={molfileContent}
                    onChange={(e) => setMolfileContent(e.target.value)}
                    placeholder="Paste Molfile content here..."
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-green-500 focus:border-green-500 transition-all duration-200"
                    rows={6}
                  />
                </div>

                <div className="flex flex-col space-y-3">
                  <div className="relative">
                    <input
                      type="file"
                      id="molfile-upload"
                      accept=".mol,.sdf,.mdl"
                      onChange={handleMolfileUpload}
                      className="hidden"
                    />
                    <label
                      htmlFor="molfile-upload"
                      className="cursor-pointer flex items-center justify-center w-full px-4 py-2.5 bg-gray-100 dark:bg-gray-700 text-gray-800 dark:text-gray-200 rounded-lg hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors duration-300"
                    >
                      <HiOutlineDocumentDownload className="w-5 h-5 mr-2" />
                      Upload Molfile
                    </label>
                  </div>

                  <button
                    onClick={convertMolfileToInChI}
                    disabled={
                      isLoading || !molfileContent.trim() || !inchiModuleLoaded
                    }
                    className={`relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-green-500 ${
                      isLoading || !molfileContent.trim() || !inchiModuleLoaded
                        ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                        : "bg-gradient-to-r from-green-600 to-emerald-700 hover:from-green-700 hover:to-emerald-800"
                    }`}
                  >
                    {isLoading
                      ? "Converting..."
                      : !inchiModuleLoaded
                      ? "Loading InChI Module..."
                      : "Convert to InChI"}
                  </button>
                </div>
              </div>
            </div>
          )}

          {/* Status and Information */}
          <div
            className={`px-6 py-4 rounded-xl flex items-center justify-between ${
              isEditorReady && inchiModuleLoaded
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
                    inchiModuleLoaded
                      ? "bg-green-500 animate-pulse"
                      : "bg-yellow-500 animate-pulse"
                  }`}
                ></div>
                <span
                  className={`text-sm font-medium ${
                    inchiModuleLoaded
                      ? "text-green-800 dark:text-green-300"
                      : "text-yellow-800 dark:text-yellow-300"
                  }`}
                >
                  {inchiModuleLoaded
                    ? "InChI Module Ready"
                    : "Loading InChI Module..."}
                </span>
              </div>
            </div>

            <div className="flex gap-3">
              {(!isEditorReady || !inchiModuleLoaded) && (
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
              About This Tool
            </h4>
            <div className="space-y-3 text-gray-700 dark:text-gray-300">
              <p>
                This tool converts between chemical structures and InChI
                notation using the official InChI algorithm.
              </p>
              <div>
                <h5 className="font-medium mb-1 text-gray-800 dark:text-gray-200">
                  Features:
                </h5>
                <ul className="list-disc list-inside space-y-1 pl-1 text-gray-600 dark:text-gray-400">
                  <li>Draw chemical structures with the Ketcher editor</li>
                  <li>Convert structures to InChI, InChIKey, and AuxInfo</li>
                  <li>Load structures from InChI strings</li>
                  <li>Convert molfiles to InChI notation</li>
                  <li>Multiple InChI versions and options</li>
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

          {/* Ketcher Editor Container */}
          {(activeTab === "draw" || activeTab === "inchi") && (
            <div
              className="bg-white dark:bg-gray-800 rounded-xl shadow-xl border border-gray-200 dark:border-gray-700 overflow-hidden flex-grow"
              style={{ minHeight: "500px" }}
            >
              <iframe
                ref={ketcherFrame}
                src={`${process.env.PUBLIC_URL}/standalone/index.html`}
                title="Ketcher Editor"
                className="w-full h-full"
                style={{ minHeight: "500px" }}
                frameBorder="0"
              />
            </div>
          )}

          {/* Controls for Draw tab */}
          {activeTab === "draw" && (
            <div className="flex flex-wrap gap-3">
              <button
                onClick={clearEditor}
                disabled={isLoading || !isEditorReady}
                className={`px-4 py-2 rounded-lg text-sm font-medium transition-colors ${
                  isLoading || !isEditorReady
                    ? "bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed"
                    : "bg-white dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-700 border border-gray-300 dark:border-gray-600"
                }`}
              >
                <HiOutlineRefresh className="h-4 w-4 inline mr-1.5" />
                Clear Editor
              </button>

              <button
                onClick={generateInChI}
                disabled={isLoading || !isEditorReady || !inchiModuleLoaded}
                className={`px-4 py-2 rounded-lg text-sm font-medium flex-1 transition-all duration-300 ${
                  isLoading || !isEditorReady || !inchiModuleLoaded
                    ? "bg-gray-400 dark:bg-gray-600 text-white cursor-not-allowed"
                    : "bg-gradient-to-r from-green-600 to-emerald-700 hover:from-green-700 hover:to-emerald-800 text-white"
                }`}
              >
                {isLoading ? "Generating..." : "Generate InChI"}
              </button>
            </div>
          )}

          {/* Results */}
          {(inchi || auxInfo || inchiKey || logMessage) && (
            <div className="bg-white dark:bg-gray-800 p-5 rounded-xl border border-gray-200 dark:border-gray-700 shadow-lg">
              <h3 className="text-lg font-bold text-gray-800 dark:text-white mb-3 flex items-center">
                <HiOutlineDocumentText className="h-5 w-5 mr-2 text-green-500 dark:text-green-400" />
                Results
              </h3>

              {/* InChI */}
              <ResultBlock
                title="InChI"
                value={inchi}
                onCopy={() => copyToClipboard(inchi, "inchi")}
                copyState={copySuccess}
              />

              {/* InChIKey */}
              <ResultBlock
                title="InChIKey"
                value={inchiKey}
                onCopy={() => copyToClipboard(inchiKey, "key")}
                copyState={copySuccess}
              />

              {/* AuxInfo */}
              <ResultBlock
                title="AuxInfo"
                value={auxInfo}
                onCopy={() => copyToClipboard(auxInfo, "auxinfo")}
                copyState={copySuccess}
              />

              {/* Log messages */}
              {logMessage && (
                <div className="mt-4 p-3 bg-gray-50 dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700">
                  <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                    Log
                  </h4>
                  <p className="text-xs text-gray-600 dark:text-gray-400 font-mono whitespace-pre-wrap">
                    {logMessage}
                  </p>
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default InChIView;
