import React, { useState, useEffect, useRef } from 'react';
import {
  HiOutlineClipboardCopy, 
  HiOutlineCheck, 
  HiOutlineExclamationCircle,
  HiOutlineInformationCircle,
  HiOutlineRefresh,
  HiOutlinePencil,
  HiOutlineX
} from 'react-icons/hi';

// Add custom styles for animations
const styles = `
@keyframes fadeIn {
  from { opacity: 0; transform: translateY(10px); }
  to { opacity: 1; transform: translateY(0); }
}
.animate-fadeIn {
  animation: fadeIn 0.5s ease-out forwards;
}
`;

// Add the styles to the document
if (typeof document !== 'undefined') {
  const styleSheet = document.createElement('style');
  styleSheet.type = 'text/css';
  styleSheet.innerText = styles;
  document.head.appendChild(styleSheet);
}

const StructureDrawView = () => {
  const [smiles, setSmiles] = useState('');
  const [inputSmiles, setInputSmiles] = useState('');
  const [copySuccess, setCopySuccess] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [copyModalText, setCopyModalText] = useState('');
  const [error, setError] = useState(null);
  const [isLoading, setIsLoading] = useState(false);
  const [isEditorReady, setIsEditorReady] = useState(false);
  const [retryAttempt, setRetryAttempt] = useState(0);
  const ketcherFrame = useRef(null);
  const copyTextRef = useRef(null);
  const messageHandlers = useRef({});
  const messageId = useRef(0);

  // Examples of common molecules
  const examples = [
    { name: 'Ethanol', value: 'CCO', description: 'Alcohol' },
    { name: 'Aspirin', value: 'CC(=O)OC1=CC=CC=C1C(=O)O', description: 'Pain reliever' },
    { name: 'Caffeine', value: 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C', description: 'Stimulant' },
    { name: 'Ibuprofen', value: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', description: 'Anti-inflammatory' },
  ];

  // Reset error state when input changes
  useEffect(() => {
    if (error) setError(null);
  }, [inputSmiles, error]);

  // Auto-select text in copy modal when it appears
  useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      copyTextRef.current.select();
    }
  }, [showCopyModal]);

  // Set up message communication with the iframe
  useEffect(() => {
    // Message handler for communications from the iframe
    const handleMessage = (event) => {
      // Only accept messages from our iframe
      if (ketcherFrame.current && event.source === ketcherFrame.current.contentWindow) {
        const { id, type, status, data, error } = event.data;
        
        // Handle response to our message
        if (id && messageHandlers.current[id]) {
          const handler = messageHandlers.current[id];
          
          if (status === 'success') {
            handler.resolve(data);
          } else if (status === 'error') {
            handler.reject(new Error(error || 'Unknown error'));
          }
          
          // Remove the handler after it's been used
          delete messageHandlers.current[id];
        }
        
        // Handle initialization message
        if (type === 'ketcher-ready') {
          console.log('Ketcher ready message received');
          setIsEditorReady(true);
        }
      }
    };

    // Add listener for messages
    window.addEventListener('message', handleMessage);

    // Clean up on unmount
    return () => {
      window.removeEventListener('message', handleMessage);
      messageHandlers.current = {};
    };
  }, []);
  
  // Function to send messages to the iframe and wait for response
  const sendMessage = (type, payload = {}) => {
    return new Promise((resolve, reject) => {
      if (!ketcherFrame.current || !ketcherFrame.current.contentWindow) {
        reject(new Error('Ketcher frame not available'));
        return;
      }
      
      // Generate unique ID for this message
      const id = messageId.current++;
      
      // Store promise handlers
      messageHandlers.current[id] = { resolve, reject };
      
      // Send message to iframe
      const message = { id, type, payload };
      
      try {
        ketcherFrame.current.contentWindow.postMessage(message, '*');
      } catch (err) {
        delete messageHandlers.current[id];
        reject(new Error(`Failed to send message: ${err.message}`));
      }
      
      // Set timeout to reject if no response
      setTimeout(() => {
        if (messageHandlers.current[id]) {
          delete messageHandlers.current[id];
          reject(new Error('Timeout waiting for response'));
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
        if (typeof ketcher[command] === 'function') {
          return await ketcher[command](...args);
        }
      }
    } catch (directError) {
      console.debug(`Direct Ketcher access failed for ${command}`, directError);
      // Fall through to postMessage method
    }
    
    // Fall back to message passing
    try {
      return await sendMessage('ketcher-command', { command, args });
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
        const script = iframeDocument.createElement('script');
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
        console.log('Communication script injected into iframe');
      } catch (error) {
        console.error('Failed to inject communication script:', error);
      }
    };
    
    // Wait for iframe to load, then inject script
    if (ketcherFrame.current) {
      // Clear any previous onload
      ketcherFrame.current.onload = () => {
        console.log('Ketcher iframe loaded, injecting script');
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
  const copyToClipboard = async (text = null) => {
    const textToCopy = text || smiles;
    
    if (!textToCopy) {
      setError('No SMILES to copy. Generate SMILES first.');
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
      const textArea = document.createElement('textarea');
      textArea.value = textToCopy;
      
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
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      } else {
        throw new Error('execCommand copy failed');
      }
    } catch (err) {
      console.error('Failed to copy text:', err);
      
      // Method 3: Show a modal with text to copy manually
      setCopyModalText(textToCopy);
      setShowCopyModal(true);
    }
  };

  // Load SMILES into Ketcher
  const loadSmiles = async () => {
    if (!inputSmiles.trim()) {
      setError('Please enter a SMILES string');
      return;
    }
    
    if (!isEditorReady) {
      setError('Editor not ready. Please try again in a moment.');
      return;
    }
    
    setIsLoading(true);
    setError(null);
    
    try {
      // Use the command execution function
      await executeKetcherCommand('setMolecule', [inputSmiles]);
      console.log('SMILES loaded successfully');
      setSmiles(inputSmiles);
    } catch (err) {
      console.error('Failed to load SMILES:', err);
      setError('Invalid SMILES string or error loading structure. Please check your input.');
    } finally {
      setIsLoading(false);
    }
  };

  // Get SMILES from Ketcher
  const getSmiles = async () => {
    if (!isEditorReady) {
      setError('Editor not ready. Please try again in a moment.');
      return '';
    }
    
    setIsLoading(true);
    setError(null);
    
    try {
      // Use the command execution function
      const newSmiles = await executeKetcherCommand('getSmiles');
      
      if (!newSmiles || newSmiles === '') {
        setError('No structure drawn. Please draw a molecule first.');
        setIsLoading(false);
        return '';
      }
      
      setSmiles(newSmiles);
      console.log('Generated SMILES:', newSmiles);
      return newSmiles;
    } catch (err) {
      console.error('Failed to generate SMILES:', err);
      setError('Could not generate SMILES from the current structure');
    } finally {
      setIsLoading(false);
    }
    return '';
  };

  // Clear the editor
  const clearEditor = async () => {
    if (!isEditorReady) {
      console.warn('Editor not ready for clearing');
      return;
    }
    
    try {
      // Use the command execution function
      await executeKetcherCommand('setMolecule', ['']);
      setSmiles('');
      console.log('Editor cleared successfully');
    } catch (err) {
      console.error('Failed to clear editor:', err);
      setError('Failed to clear the editor');
    }
  };

  const handleUseExample = (exampleValue) => {
    setInputSmiles(exampleValue);
  };

  // Function to retry initialization
  const handleRetryInit = () => {
    setIsEditorReady(false);
    setError(null);
    setRetryAttempt(prev => prev + 1);
  };

  return (
    <div className="flex flex-col gap-6 p-4 md:p-6 bg-white dark:bg-gray-900 min-h-screen">
      {/* Header with animated background */}
      <div className="bg-gradient-to-r from-blue-600 to-purple-700 dark:from-blue-700 dark:to-purple-900 rounded-xl shadow-xl overflow-hidden relative">
        {/* Animated background elements */}
        <div className="absolute inset-0 overflow-hidden opacity-10">
          <div className="absolute left-0 top-0 w-40 h-40 rounded-full bg-white transform -translate-x-1/2 -translate-y-1/2"></div>
          <div className="absolute right-20 bottom-0 w-64 h-64 rounded-full bg-white transform translate-x-1/2 translate-y-1/2"></div>
          <div className="absolute left-1/3 top-1/2 w-32 h-32 rounded-full bg-white transform -translate-y-1/2"></div>
        </div>

        
      </div>

      {/* Copy Modal - For fallback copying */}
      {showCopyModal && (
        <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4 animate-fadeIn">
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
                value={copyModalText}
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

      {/* Main content area */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
        {/* Left sidebar with controls */}
        <div className="lg:col-span-3 flex flex-col gap-4">
          {/* SMILES Input Section */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-xl border border-gray-100 dark:border-gray-700 backdrop-blur-sm">
            <div className="flex items-center mb-4">
              <div className="bg-indigo-100 dark:bg-indigo-900/50 p-2 rounded-lg mr-3">
                <HiOutlinePencil className="h-5 w-5 text-indigo-700 dark:text-indigo-400" />
              </div>
              <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                SMILES Input
              </h2>
            </div>
            
            <div className="space-y-4">
              <div>
                <label htmlFor="smiles-input" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                  Enter SMILES to Edit
                </label>
                <input
                  id="smiles-input"
                  type="text"
                  value={inputSmiles}
                  onChange={(e) => setInputSmiles(e.target.value)}
                  placeholder="e.g., CCO for ethanol"
                  className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 transition-all duration-200"
                />
              </div>
              
              <button
                onClick={loadSmiles}
                disabled={isLoading || !isEditorReady || !inputSmiles.trim()}
                className={`w-full relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                  isLoading || !isEditorReady || !inputSmiles.trim()
                    ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                    : 'bg-gradient-to-r from-blue-600 to-indigo-700 hover:from-blue-700 hover:to-indigo-800'
                }`}
              >
                {isLoading ? 'Loading...' : !isEditorReady ? 'Initializing...' : 'Load Structure'}
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
                    onClick={() => handleUseExample(example.value)}
                    className="inline-flex items-center px-3 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded-full text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-indigo-500 dark:focus:ring-offset-gray-800 transition-all duration-200"
                    title={`${example.name}: ${example.description}`}
                  >
                    {example.name}
                  </button>
                ))}
              </div>
            </div>
          </div>
          
          {/* Generated SMILES Section */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-xl border border-gray-100 dark:border-gray-700 backdrop-blur-sm">
            <div className="flex items-center mb-4">
              <div className="bg-green-100 dark:bg-green-900/50 p-2 rounded-lg mr-3">
                <HiOutlineClipboardCopy className="h-5 w-5 text-green-700 dark:text-green-400" />
              </div>
              <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                Generated SMILES
              </h2>
            </div>
            
            <div className="space-y-4">
              <div className="bg-gray-50 dark:bg-gray-900/50 p-4 rounded-lg border border-gray-200 dark:border-gray-700">
                <p className="text-sm text-gray-500 dark:text-gray-400 mb-2">Structure as SMILES:</p>
                <div className="font-mono text-sm bg-white dark:bg-gray-800 p-3 rounded border border-gray-200 dark:border-gray-700 break-all max-h-32 overflow-y-auto">
                  {smiles || <span className="text-gray-400 dark:text-gray-500">No structure generated yet</span>}
                </div>
              </div>
              
              <div className="flex flex-col sm:flex-row gap-3">
                <button
                  onClick={getSmiles}
                  disabled={isLoading || !isEditorReady}
                  className={`flex-1 relative overflow-hidden px-4 py-2.5 rounded-lg font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500 ${
                    isLoading || !isEditorReady
                      ? 'bg-gray-400 dark:bg-gray-600 text-white cursor-not-allowed'
                      : 'bg-blue-600 text-white hover:bg-blue-700'
                  }`}
                >
                  Generate SMILES
                </button>
                
                <button
                  onClick={() => copyToClipboard()}
                  disabled={!smiles}
                  className={`flex-1 relative overflow-hidden px-4 py-2.5 rounded-lg font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-indigo-500 ${
                    !smiles
                      ? 'bg-gray-400 dark:bg-gray-600 text-white cursor-not-allowed'
                      : 'bg-indigo-600 text-white hover:bg-indigo-700'
                  }`}
                >
                  {copySuccess ? (
                    <>
                      <HiOutlineCheck className="w-5 h-5 mr-1" />
                      Copied!
                    </>
                  ) : (
                    <>
                      <HiOutlineClipboardCopy className="w-5 h-5 mr-1" />
                      Copy
                    </>
                  )}
                </button>
              </div>
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
                This structure editor allows you to draw chemical structures and generate their SMILES notation.
              </p>
              <div>
                <h5 className="font-medium mb-1 text-gray-800 dark:text-gray-200">Features:</h5>
                <ul className="list-disc list-inside space-y-1 pl-1 text-gray-600 dark:text-gray-400">
                  <li>Draw chemical structures with an intuitive interface</li>
                  <li>Load existing structures from SMILES strings</li>
                  <li>Generate SMILES notation from drawn structures</li>
                  <li>Copy results to clipboard for use in other applications</li>
                </ul>
              </div>
            </div>
          </div>
        </div>
        
        {/* Main Editor Area */}
        <div className="lg:col-span-9 flex flex-col gap-4">
          {/* Editor Status */}
          <div className={`px-6 py-3 rounded-xl flex items-center justify-between ${
            isEditorReady 
              ? 'bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800/40' 
              : 'bg-yellow-50 dark:bg-yellow-900/20 border border-yellow-200 dark:border-yellow-800/40'
          }`}>
            <div className="flex items-center">
              <div className={`h-3 w-3 rounded-full mr-3 ${
                isEditorReady 
                  ? 'bg-green-500 animate-pulse' 
                  : 'bg-yellow-500 animate-pulse'
              }`}></div>
              <span className={`font-medium ${
                isEditorReady 
                  ? 'text-green-800 dark:text-green-300' 
                  : 'text-yellow-800 dark:text-yellow-300'
              }`}>
                {isEditorReady ? 'Editor Ready' : 'Initializing Editor...'}
              </span>
            </div>
            
            <div className="flex gap-3">
              {!isEditorReady && (
                <button
                  onClick={handleRetryInit}
                  className="inline-flex items-center px-3 py-1.5 rounded-lg text-sm font-medium transition-colors bg-yellow-100 dark:bg-yellow-900/30 text-yellow-700 dark:text-yellow-300 hover:bg-yellow-200 dark:hover:bg-yellow-900/50 border border-yellow-300 dark:border-yellow-700/50"
                >
                  <HiOutlineRefresh className="h-4 w-4 mr-1.5" />
                  Retry Init
                </button>
              )}
            
              <button
                onClick={clearEditor}
                disabled={isLoading || !isEditorReady}
                className={`inline-flex items-center px-3 py-1.5 rounded-lg text-sm font-medium transition-colors ${
                  isLoading || !isEditorReady
                    ? 'bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed'
                    : 'bg-white dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-700 border border-gray-300 dark:border-gray-600'
                }`}
              >
                <HiOutlineRefresh className="h-4 w-4 mr-1.5" />
                Clear Editor
              </button>
            </div>
          </div>
          
          {/* Error Display */}
          {error && (
            <div className="p-4 rounded-xl bg-red-50 dark:bg-red-900/30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700/50 flex items-start shadow-lg animate-fadeIn" role="alert">
              <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
              <p>{error}</p>
            </div>
          )}
          
          {/* Ketcher Editor Container */}
          <div className="bg-white dark:bg-gray-800 rounded-xl shadow-xl border border-gray-200 dark:border-gray-700 overflow-hidden flex-grow" style={{ minHeight: "600px" }}>
            <iframe
              ref={ketcherFrame}
              src={`${process.env.PUBLIC_URL}/standalone/index.html`}
              title="Ketcher Editor"
              className="w-full h-full"
              style={{ minHeight: "600px" }}
              frameBorder="0"
            />
          </div>
          
          {/* Instructions */}
          <div className="bg-white dark:bg-gray-800 p-5 rounded-xl border border-gray-200 dark:border-gray-700 shadow-lg">
            <h3 className="text-lg font-bold text-gray-800 dark:text-white mb-3 flex items-center">
              <HiOutlineInformationCircle className="h-5 w-5 mr-2 text-blue-500 dark:text-blue-400" />
              How to Use
            </h3>
            <ol className="list-decimal list-inside space-y-2 text-gray-700 dark:text-gray-300">
              <li>Draw a molecular structure using the editor tools</li>
              <li>Click "Generate SMILES" to convert your structure to SMILES notation</li>
              <li>Use "Copy" to save the SMILES to your clipboard</li>
              <li>Alternatively, paste a SMILES string in the input field and click "Load Structure"</li>
            </ol>
          </div>
        </div>
      </div>
    </div>
  );
};

export default StructureDrawView;