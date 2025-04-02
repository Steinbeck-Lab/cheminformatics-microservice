// Description: This component handles the Optical Chemical Structure Recognition (OCSR) functionality.
import React, { useState, useCallback, useEffect } from 'react'; // Added useEffect
import { useDropzone } from 'react-dropzone';
// Ensure all used icons are imported
import {
  HiOutlineUpload, // Added for dropzone hint
  HiOutlineLink, // Added for URL input section
  HiOutlineExclamationCircle // Added for error display
} from 'react-icons/hi';
// Assuming these components are correctly implemented and styled for dark/light mode
import MoleculeCard from '../common/MoleculeCard';
// import LoadingScreen from '../common/LoadingScreen'; // Removed unused import
// Assuming this service is configured correctly
import ocsrService from '../../services/ocsrService'; // Assuming this service exists

const OCRView = () => {
  const [files, setFiles] = useState([]); // Stores the uploaded file object
  const [imageUrl, setImageUrl] = useState(''); // Stores the image URL input
  const [results, setResults] = useState([]); // Stores detected SMILES strings
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [reference, setReference] = useState(''); // Optional reference input
  const [showUrlInput, setShowUrlInput] = useState(false); // Toggle between upload and URL input

  // Callback for react-dropzone
  const onDrop = useCallback((acceptedFiles, rejectedFiles) => { // Added rejectedFiles
    // Clear previous state on new drop attempt
    setFiles([]);
    setImageUrl('');
    setError(null);
    setResults([]);

    if (acceptedFiles.length > 0) {
      setFiles(acceptedFiles);
      console.log("Accepted file:", acceptedFiles[0].name);
    } else if (rejectedFiles.length > 0) {
      // Handle rejected files (e.g., wrong type, too many files)
      console.error("Rejected file:", rejectedFiles[0].errors);
      const errorMessages = rejectedFiles[0].errors.map(e => e.message).join(', ');
      setError(`File rejected: ${errorMessages}. Please upload a single image (PNG, JPG, GIF).`);
    }
  }, []);

  // Setup react-dropzone
  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: { // Define accepted MIME types and extensions
      'image/png': ['.png'],
      'image/jpeg': ['.jpeg', '.jpg'],
      'image/gif': ['.gif'],
    },
    maxFiles: 1, // Allow only one file
    multiple: false
  });

  // Handle image processing (either uploaded file or URL)
  const handleProcessImage = async () => {
    const hasFile = files.length > 0;
    const hasUrl = imageUrl.trim() !== '';

    if (!hasFile && !hasUrl) {
      setError('Please upload an image or provide an image URL.');
      return;
    }

    setLoading(true);
    setError(null);
    setResults([]); // Clear previous results

    try {
      let response;
      // Prioritize file upload if both are present
      if (hasFile) {
        // Ensure service function exists
        if (!ocsrService || typeof ocsrService.processImage !== 'function') {
          throw new Error("Image processing service (file) is not available.");
        }
        response = await ocsrService.processImage(files[0], reference || undefined); // Pass undefined if empty
      } else {
        // Ensure service function exists
        if (!ocsrService || typeof ocsrService.processImageUrl !== 'function') {
          throw new Error("Image processing service (URL) is not available.");
        }
        response = await ocsrService.processImageUrl(imageUrl.trim(), reference || undefined); // Pass undefined if empty
      }

      // Process the response
      // Assuming response structure is { smiles: "...", ... } or { smiles: ["...", "..."], ... }
      if (response && response.smiles) {
        const smilesArray = Array.isArray(response.smiles) ? response.smiles : [response.smiles];
        // Filter out potentially empty strings if API returns them
        const validSmiles = smilesArray.filter(s => typeof s === 'string' && s.trim() !== '');

        if (validSmiles.length > 0) {
          setResults(validSmiles);
        } else {
          setError('OCSR process completed, but no chemical structures were detected in the image.');
        }
      } else {
        // Handle cases where response might be missing 'smiles' key
        setError('OCSR process completed, but no chemical structures were detected.');
      }
    } catch (err) {
      console.error("Error processing image:", err);
      setError(`Error processing image: ${err.response?.data?.detail || err.message || 'Unknown error'}`);
      setResults([]); // Ensure results are cleared on error
    } finally {
      setLoading(false);
    }
  };

  // Clear all inputs and results
  const clearAll = () => {
    setFiles([]);
    setImageUrl('');
    setResults([]);
    setError(null);
    setReference('');
    // Optionally reset showUrlInput state
    // setShowUrlInput(false);
  };

  // Get preview URL for uploaded file
  const filePreviewUrl = files.length > 0 ? URL.createObjectURL(files[0]) : null;

  // --- Clean up object URL when component unmounts or file changes ---
  useEffect(() => {
    // This is the cleanup function
    return () => {
      if (filePreviewUrl) {
        console.log("Revoking Object URL:", filePreviewUrl); // Debug log
        URL.revokeObjectURL(filePreviewUrl);
      }
    };
  }, [filePreviewUrl]); // Dependency: run cleanup when the URL changes


  return (
    // Main container with theme-aware background
    <div className="space-y-6 p-4 md:p-6 bg-gray-50 dark:bg-gray-900 min-h-screen">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        {/* Title */}
        <h2 className="text-2xl font-bold text-gray-800 dark:text-blue-400 mb-4">
          Optical Chemical Structure Recognition (OCSR)
        </h2>

        {/* Input Method Tabs */}
        <div className="space-y-4">
          <div className="flex border-b border-gray-200 dark:border-gray-700 mb-4">
            {/* Upload Tab Button */}
            <button
              onClick={() => setShowUrlInput(false)}
              className={`px-4 py-2 text-sm font-medium focus:outline-none transition-colors duration-150 ${!showUrlInput
                  ? 'border-b-2 border-blue-500 text-blue-600 dark:border-blue-400 dark:text-blue-400'
                  : 'border-b-2 border-transparent text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-200 hover:border-gray-300 dark:hover:border-gray-600'
                }`}
              aria-current={!showUrlInput ? 'page' : undefined}
            >
              Upload Image
            </button>
            {/* URL Tab Button */}
            <button
              onClick={() => setShowUrlInput(true)}
              className={`px-4 py-2 text-sm font-medium focus:outline-none transition-colors duration-150 ${showUrlInput
                  ? 'border-b-2 border-blue-500 text-blue-600 dark:border-blue-400 dark:text-blue-400'
                  : 'border-b-2 border-transparent text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-200 hover:border-gray-300 dark:hover:border-gray-600'
                }`}
              aria-current={showUrlInput ? 'page' : undefined}
            >
              Use Image URL
            </button>
          </div>

          {/* Conditional Input Area */}
          {!showUrlInput ? (
            // Dropzone Styling
            <div
              {...getRootProps()}
              className={`border-2 border-dashed rounded-lg p-6 text-center cursor-pointer transition-colors duration-200 ease-in-out ${isDragActive
                  ? 'border-blue-500 dark:border-blue-400 bg-blue-50 dark:bg-blue-900/20' // Active drag state
                  : 'border-gray-300 dark:border-gray-600 hover:border-blue-400 dark:hover:border-blue-500 bg-white dark:bg-gray-800 hover:bg-gray-50 dark:hover:bg-gray-700/50' // Default state
                }`}
            >
              <input {...getInputProps()} />
              {filePreviewUrl ? (
                // File Preview
                <div className="flex flex-col items-center space-y-2">
                  <img
                    src={filePreviewUrl}
                    alt="Preview"
                    className="max-h-48 rounded-md border border-gray-200 dark:border-gray-700"
                  />
                  <p className="text-sm text-gray-700 dark:text-gray-300">{files[0]?.name || 'Uploaded Image'}</p>
                </div>
              ) : (
                // Dropzone Hint Text
                <div className="text-gray-500 dark:text-gray-400 flex flex-col items-center">
                  <HiOutlineUpload className="h-10 w-10 mb-2 text-gray-400 dark:text-gray-500" aria-hidden="true" />
                  <p className="text-base mb-1">
                    {isDragActive
                      ? 'Drop the image here...'
                      : 'Drag & drop image, or click to select'}
                  </p>
                  <p className="text-xs">
                    Supports PNG, JPG, GIF
                  </p>
                </div>
              )}
            </div>
          ) : (
            // URL Input Styling
            <div className="space-y-2">
              <label htmlFor="imageUrl" className="block text-sm font-medium text-gray-700 dark:text-gray-300">
                Image URL
              </label>
              <div className="relative">
                <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                  <HiOutlineLink className="h-5 w-5 text-gray-400 dark:text-gray-500" aria-hidden="true" />
                </div>
                <input
                  type="url"
                  id="imageUrl"
                  value={imageUrl}
                  onChange={(e) => { setImageUrl(e.target.value); setFiles([]); setError(null); setResults([]); }} // Clear file/results on URL change
                  placeholder="https://example.com/structure.png"
                  className="w-full pl-10 pr-4 py-2 rounded-md bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 text-gray-900 dark:text-white focus:outline-none focus:border-indigo-500 dark:focus:border-blue-500 focus:ring-1 focus:ring-indigo-500 dark:focus:ring-blue-500 shadow-sm"
                />
              </div>
            </div>
          )}

          {/* Optional Reference Input */}
          <div className="space-y-2 pt-2">
            <label htmlFor="reference" className="block text-sm font-medium text-gray-700 dark:text-gray-300">
              Reference Name (Optional)
            </label>
            <input
              type="text"
              id="reference"
              value={reference}
              onChange={(e) => setReference(e.target.value)}
              placeholder="e.g., Figure 1A, Compound X"
              className="w-full px-4 py-2 rounded-md bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 text-gray-900 dark:text-white focus:outline-none focus:border-indigo-500 dark:focus:border-blue-500 focus:ring-1 focus:ring-indigo-500 dark:focus:ring-blue-500 shadow-sm"
            />
          </div>

          {/* Action Buttons */}
          <div className="flex flex-wrap gap-3 pt-4 border-t border-gray-200 dark:border-gray-700">
            {/* Process Button */}
            <button
              type="button" // Changed from submit to prevent default form action
              onClick={handleProcessImage}
              disabled={loading || (!files.length && !imageUrl.trim())}
              // Button Styling
              className={`px-5 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${loading || (!files.length && !imageUrl.trim())
                  ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed' // Disabled
                  : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm' // Enabled
                }`}
            >
              {loading ? (
                <>
                  <svg className="animate-spin -ml-1 mr-3 h-5 w-5 text-white" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                  </svg>
                  Processing...
                </>
              ) : 'Process Image'}
            </button>

            {/* Clear Button */}
            {(files.length > 0 || imageUrl.trim() || results.length > 0 || error) && (
              <button
                type="button"
                onClick={clearAll}
                // Button Styling (Red color)
                className="px-5 py-2 rounded-lg font-medium bg-red-600 text-white hover:bg-red-700 dark:bg-red-500 dark:hover:bg-red-600 shadow-sm transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-red-500"
              >
                Clear
              </button>
            )}
          </div>
        </div>
      </div>

      {/* Loading State Full Screen (Optional, could use inline indicator instead) */}
      {/* {loading && <LoadingScreen text="Processing chemical structure..." />} */}

      {/* Error Display */}
      {error && !loading && (
        // Error message styling
        <div className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow" role="alert">
          <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
          <span>{error}</span>
        </div>
      )}

      {/* Results Display Section */}
      {results.length > 0 && !loading && (
        <div className="space-y-4">
          {/* Results Header */}
          <h3 className="text-xl font-semibold text-gray-800 dark:text-blue-300">
            Detected Structures ({results.length})
          </h3>
          {/* Results Grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
            {results.map((smilesResult, index) => (
              // Assuming MoleculeCard is theme-aware
              <MoleculeCard
                key={index}
                smiles={smilesResult}
                title={`Detected Structure ${index + 1}`}
                showActions={true} // Show copy/download actions on result cards
                size="md" // Adjust size as needed
              />
            ))}
          </div>
        </div>
      )}
    </div>
  );
};

export default OCRView;
