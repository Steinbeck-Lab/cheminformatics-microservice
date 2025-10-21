// Description: Enhanced OCRView component with production-ready UX/UI
import React, { useState, useCallback, useEffect } from "react";
import { useDropzone } from "react-dropzone";
import { motion, AnimatePresence } from "framer-motion";
import {
  HiOutlineUpload,
  HiOutlineLink,
  HiOutlineExclamationCircle,
  HiOutlineCheckCircle,
  HiOutlineXCircle,
  HiOutlinePhotograph,
  HiOutlineRefresh,
} from "react-icons/hi";
import MoleculeCard from "../common/MoleculeCard";
import ocsrService from "../../services/ocsrService";

const OCRView = () => {
  // State management
  const [files, setFiles] = useState([]);
  const [imageUrl, setImageUrl] = useState("");
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [reference, setReference] = useState("");
  const [inputMethod, setInputMethod] = useState("upload"); // 'upload' | 'url'
  const [handDrawn, setHandDrawn] = useState(false);
  const [processingStage, setProcessingStage] = useState("");

  // Dropzone callback
  const onDrop = useCallback((acceptedFiles, rejectedFiles) => {
    setFiles([]);
    setImageUrl("");
    setError(null);
    setResults([]);

    if (acceptedFiles.length > 0) {
      setFiles(acceptedFiles);
    } else if (rejectedFiles.length > 0) {
      const errorMessages = rejectedFiles[0].errors.map((e) => e.message).join(", ");
      setError(
        `File rejected: ${errorMessages}. Please upload a valid image file (PNG, JPG, or GIF).`
      );
    }
  }, []);

  // Configure dropzone
  const { getRootProps, getInputProps, isDragActive, isDragReject } = useDropzone({
    onDrop,
    accept: {
      "image/png": [".png"],
      "image/jpeg": [".jpeg", ".jpg"],
      "image/gif": [".gif"],
    },
    maxFiles: 1,
    multiple: false,
    disabled: loading,
  });

  // Image processing handler
  const handleProcessImage = async () => {
    const hasFile = files.length > 0;
    const hasUrl = imageUrl.trim() !== "";

    if (!hasFile && !hasUrl) {
      setError("Please upload an image or provide an image URL.");
      return;
    }

    setLoading(true);
    setError(null);
    setResults([]);
    setProcessingStage("Uploading image...");

    try {
      let response;

      if (hasFile) {
        setProcessingStage("Analyzing structure...");
        if (!ocsrService || typeof ocsrService.processImage !== "function") {
          throw new Error("Image processing service is not available.");
        }
        response = await ocsrService.processImage(files[0], reference || undefined, handDrawn);
      } else {
        setProcessingStage("Fetching image...");
        if (!ocsrService || typeof ocsrService.processImageUrl !== "function") {
          throw new Error("URL processing service is not available.");
        }
        response = await ocsrService.processImageUrl(
          imageUrl.trim(),
          reference || undefined,
          handDrawn
        );
      }

      setProcessingStage("Extracting structures...");

      if (response && response.smiles) {
        const smilesArray = Array.isArray(response.smiles) ? response.smiles : [response.smiles];
        const validSmiles = smilesArray.filter((s) => typeof s === "string" && s.trim() !== "");

        if (validSmiles.length > 0) {
          setResults(validSmiles);
          setProcessingStage("");
        } else {
          setError(
            "No chemical structures were detected in the image. Please try a different image or adjust the settings."
          );
        }
      } else {
        setError("Processing completed, but no chemical structures were detected.");
      }
    } catch (err) {
      console.error("Error processing image:", err);
      setError(`Error: ${err.response?.data?.detail || err.message || "Unable to process image"}`);
      setResults([]);
    } finally {
      setLoading(false);
      setProcessingStage("");
    }
  };

  // Clear all state
  const clearAll = () => {
    setFiles([]);
    setImageUrl("");
    setResults([]);
    setError(null);
    setReference("");
    setProcessingStage("");
  };

  // Get preview URL for uploaded file
  const filePreviewUrl = files.length > 0 ? URL.createObjectURL(files[0]) : null;
  const displayImageUrl = filePreviewUrl || (imageUrl.trim() ? imageUrl : null);

  // Cleanup object URL
  useEffect(() => {
    return () => {
      if (filePreviewUrl) {
        URL.revokeObjectURL(filePreviewUrl);
      }
    };
  }, [filePreviewUrl]);

  return (
    <div className="space-y-6">
      {/* Main Input Card */}
      <div className="bg-white dark:bg-slate-800 rounded-xl shadow-lg border border-slate-200 dark:border-slate-700 overflow-hidden">
        {/* Input Method Selector */}
        <div className="border-b border-slate-200 dark:border-slate-700 bg-slate-50 dark:bg-slate-900/50">
          <div className="flex">
            <button
              onClick={() => {
                setInputMethod("upload");
                setImageUrl("");
                setError(null);
              }}
              disabled={loading}
              className={`flex-1 px-6 py-4 text-sm font-medium transition-all duration-200 ${
                inputMethod === "upload"
                  ? "bg-white dark:bg-slate-800 text-blue-600 dark:text-blue-400 border-b-2 border-blue-600 dark:border-blue-400"
                  : "text-slate-600 dark:text-slate-400 hover:text-slate-900 dark:hover:text-slate-200 hover:bg-slate-100 dark:hover:bg-slate-800/50"
              }`}
            >
              <div className="flex items-center justify-center gap-2">
                <HiOutlineUpload className="h-5 w-5" />
                <span>Upload Image</span>
              </div>
            </button>

            <button
              onClick={() => {
                setInputMethod("url");
                setFiles([]);
                setError(null);
              }}
              disabled={loading}
              className={`flex-1 px-6 py-4 text-sm font-medium transition-all duration-200 ${
                inputMethod === "url"
                  ? "bg-white dark:bg-slate-800 text-blue-600 dark:text-blue-400 border-b-2 border-blue-600 dark:border-blue-400"
                  : "text-slate-600 dark:text-slate-400 hover:text-slate-900 dark:hover:text-slate-200 hover:bg-slate-100 dark:hover:bg-slate-800/50"
              }`}
            >
              <div className="flex items-center justify-center gap-2">
                <HiOutlineLink className="h-5 w-5" />
                <span>Use URL</span>
              </div>
            </button>
          </div>
        </div>

        {/* Input Area */}
        <div className="p-6 space-y-6">
          {/* Upload or URL Input */}
          {inputMethod === "upload" ? (
            <div
              {...getRootProps()}
              className={`relative border-2 border-dashed rounded-lg p-8 transition-all duration-200 cursor-pointer ${
                isDragActive && !isDragReject
                  ? "border-blue-500 bg-blue-50 dark:bg-blue-900/20"
                  : isDragReject
                    ? "border-red-500 bg-red-50 dark:bg-red-900/20"
                    : files.length > 0
                      ? "border-green-500 bg-green-50 dark:bg-green-900/20"
                      : "border-slate-300 dark:border-slate-600 hover:border-blue-400 dark:hover:border-blue-500 bg-slate-50 dark:bg-slate-900/50"
              } ${loading ? "opacity-50 cursor-not-allowed" : ""}`}
            >
              <input {...getInputProps()} />

              <div className="flex flex-col items-center justify-center text-center space-y-3">
                {files.length > 0 ? (
                  <>
                    <HiOutlineCheckCircle className="h-12 w-12 text-green-600 dark:text-green-400" />
                    <div>
                      <p className="text-sm font-medium text-green-700 dark:text-green-300">
                        File selected: {files[0].name}
                      </p>
                      <p className="text-xs text-slate-500 dark:text-slate-400 mt-1">
                        {(files[0].size / 1024).toFixed(2)} KB
                      </p>
                    </div>
                  </>
                ) : isDragActive ? (
                  <>
                    <HiOutlineUpload className="h-12 w-12 text-blue-600 dark:text-blue-400 animate-bounce" />
                    <p className="text-sm font-medium text-blue-700 dark:text-blue-300">
                      Drop your image here
                    </p>
                  </>
                ) : (
                  <>
                    <HiOutlinePhotograph className="h-12 w-12 text-slate-400 dark:text-slate-500" />
                    <div>
                      <p className="text-sm font-medium text-slate-700 dark:text-slate-300">
                        Click to upload or drag and drop
                      </p>
                      <p className="text-xs text-slate-500 dark:text-slate-400 mt-1">
                        PNG, JPG, or GIF (max 10MB)
                      </p>
                    </div>
                  </>
                )}
              </div>
            </div>
          ) : (
            <div className="space-y-2">
              <label
                htmlFor="image-url"
                className="block text-sm font-medium text-slate-700 dark:text-slate-300"
              >
                Image URL
              </label>
              <div className="relative">
                <HiOutlineLink className="absolute left-3 top-1/2 transform -translate-y-1/2 h-5 w-5 text-slate-400" />
                <input
                  id="image-url"
                  type="url"
                  value={imageUrl}
                  onChange={(e) => {
                    setImageUrl(e.target.value);
                    setFiles([]);
                    setError(null);
                    setResults([]);
                  }}
                  placeholder="https://example.com/structure.png"
                  disabled={loading}
                  className="w-full pl-10 pr-4 py-3 rounded-lg bg-white dark:bg-slate-900 border border-slate-300 dark:border-slate-600 text-slate-900 dark:text-white placeholder-slate-400 focus:outline-none focus:ring-2 focus:ring-blue-500 dark:focus:ring-blue-400 focus:border-transparent transition-all disabled:opacity-50 disabled:cursor-not-allowed"
                />
              </div>
            </div>
          )}

          {/* Image Preview */}
          <AnimatePresence>
            {displayImageUrl && (
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: "auto" }}
                exit={{ opacity: 0, height: 0 }}
                className="overflow-hidden"
              >
                <div className="relative rounded-lg overflow-hidden bg-slate-100 dark:bg-slate-900 border border-slate-200 dark:border-slate-700">
                  <img
                    src={displayImageUrl}
                    alt="Structure preview"
                    className="w-full h-auto max-h-96 object-contain"
                    onError={() => setError("Failed to load image. Please check the URL.")}
                  />
                  <div className="absolute top-2 right-2">
                    <button
                      onClick={(e) => {
                        e.stopPropagation();
                        clearAll();
                      }}
                      className="p-2 bg-red-600 hover:bg-red-700 text-white rounded-lg shadow-lg transition-colors"
                      title="Remove image"
                    >
                      <HiOutlineXCircle className="h-5 w-5" />
                    </button>
                  </div>
                </div>
              </motion.div>
            )}
          </AnimatePresence>

          {/* Reference Input */}
          <div className="space-y-2">
            <label
              htmlFor="reference"
              className="block text-sm font-medium text-slate-700 dark:text-slate-300"
            >
              Reference Name{" "}
              <span className="text-slate-500 dark:text-slate-400 font-normal">(Optional)</span>
            </label>
            <input
              id="reference"
              type="text"
              value={reference}
              onChange={(e) => setReference(e.target.value)}
              placeholder="e.g., Figure 1A, Compound X"
              disabled={loading}
              className="w-full px-4 py-3 rounded-lg bg-white dark:bg-slate-900 border border-slate-300 dark:border-slate-600 text-slate-900 dark:text-white placeholder-slate-400 focus:outline-none focus:ring-2 focus:ring-blue-500 dark:focus:ring-blue-400 focus:border-transparent transition-all disabled:opacity-50 disabled:cursor-not-allowed"
            />
          </div>

          {/* Controls */}
          <div className="flex flex-wrap items-center gap-4 pt-4 border-t border-slate-200 dark:border-slate-700">
            {/* Hand-Drawn Toggle */}
            <div className="flex items-center gap-3">
              <span className="text-sm font-medium text-slate-700 dark:text-slate-300">
                Hand-Drawn Model
              </span>
              <button
                type="button"
                onClick={() => setHandDrawn(!handDrawn)}
                disabled={loading}
                className={`relative inline-flex h-6 w-11 items-center rounded-full transition-all duration-200 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 dark:focus:ring-offset-slate-800 disabled:opacity-50 disabled:cursor-not-allowed ${
                  handDrawn ? "bg-green-600 dark:bg-green-500" : "bg-slate-300 dark:bg-slate-600"
                }`}
                role="switch"
                aria-checked={handDrawn}
                aria-label="Toggle hand-drawn model"
              >
                <span
                  className={`inline-block h-4 w-4 transform rounded-full bg-white shadow-lg transition-transform duration-200 ${
                    handDrawn ? "translate-x-6" : "translate-x-1"
                  }`}
                />
              </button>
            </div>

            <div className="flex-1" />

            {/* Action Buttons */}
            <div className="flex gap-3">
              {(files.length > 0 || imageUrl.trim() || results.length > 0 || error) && (
                <button
                  type="button"
                  onClick={clearAll}
                  disabled={loading}
                  className="px-5 py-2.5 rounded-lg font-medium bg-slate-200 hover:bg-slate-300 dark:bg-slate-700 dark:hover:bg-slate-600 text-slate-700 dark:text-slate-200 transition-colors focus:outline-none focus:ring-2 focus:ring-slate-500 focus:ring-offset-2 dark:focus:ring-offset-slate-800 disabled:opacity-50 disabled:cursor-not-allowed"
                >
                  <div className="flex items-center gap-2">
                    <HiOutlineRefresh className="h-4 w-4" />
                    <span>Reset</span>
                  </div>
                </button>
              )}

              <button
                type="button"
                onClick={handleProcessImage}
                disabled={loading || (!files.length && !imageUrl.trim())}
                className={`px-6 py-2.5 rounded-lg font-medium text-white transition-all focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 dark:focus:ring-offset-slate-800 ${
                  loading || (!files.length && !imageUrl.trim())
                    ? "bg-slate-400 dark:bg-slate-600 cursor-not-allowed"
                    : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-md hover:shadow-lg"
                }`}
              >
                {loading ? (
                  <div className="flex items-center gap-2">
                    <svg
                      className="animate-spin h-5 w-5 text-white"
                      xmlns="http://www.w3.org/2000/svg"
                      fill="none"
                      viewBox="0 0 24 24"
                    >
                      <circle
                        className="opacity-25"
                        cx="12"
                        cy="12"
                        r="10"
                        stroke="currentColor"
                        strokeWidth="4"
                      />
                      <path
                        className="opacity-75"
                        fill="currentColor"
                        d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                      />
                    </svg>
                    <span>Processing...</span>
                  </div>
                ) : (
                  "Extract Structures"
                )}
              </button>
            </div>
          </div>

          {/* Processing Stage Indicator */}
          <AnimatePresence>
            {loading && processingStage && (
              <motion.div
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                exit={{ opacity: 0, y: -10 }}
                className="flex items-center justify-center gap-2 text-sm text-blue-600 dark:text-blue-400"
              >
                <div className="flex gap-1">
                  <span
                    className="w-2 h-2 bg-blue-600 dark:bg-blue-400 rounded-full animate-bounce"
                    style={{ animationDelay: "0ms" }}
                  />
                  <span
                    className="w-2 h-2 bg-blue-600 dark:bg-blue-400 rounded-full animate-bounce"
                    style={{ animationDelay: "150ms" }}
                  />
                  <span
                    className="w-2 h-2 bg-blue-600 dark:bg-blue-400 rounded-full animate-bounce"
                    style={{ animationDelay: "300ms" }}
                  />
                </div>
                <span>{processingStage}</span>
              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </div>

      {/* Error Display */}
      <AnimatePresence>
        {error && !loading && (
          <motion.div
            initial={{ opacity: 0, y: -10 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -10 }}
            className="bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-lg p-4"
            role="alert"
          >
            <div className="flex items-start gap-3">
              <HiOutlineExclamationCircle className="h-6 w-6 text-red-600 dark:text-red-400 flex-shrink-0 mt-0.5" />
              <div className="flex-1">
                <h3 className="text-sm font-medium text-red-800 dark:text-red-300 mb-1">
                  Processing Error
                </h3>
                <p className="text-sm text-red-700 dark:text-red-200">{error}</p>
              </div>
            </div>
          </motion.div>
        )}
      </AnimatePresence>

      {/* Results Section */}
      <AnimatePresence>
        {results.length > 0 && !loading && (
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -20 }}
            className="space-y-4"
          >
            {/* Results Header */}
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-2">
                <HiOutlineCheckCircle className="h-6 w-6 text-green-600 dark:text-green-400" />
                <h3 className="text-xl font-semibold text-slate-900 dark:text-white">
                  Detected Structures
                </h3>
                <span className="px-2.5 py-0.5 rounded-full text-xs font-medium bg-green-100 dark:bg-green-900/30 text-green-700 dark:text-green-300">
                  {results.length} {results.length === 1 ? "structure" : "structures"}
                </span>
              </div>
            </div>

            {/* Results Grid */}
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
              {results.map((smilesResult, index) => (
                <motion.div
                  key={index}
                  initial={{ opacity: 0, scale: 0.95 }}
                  animate={{ opacity: 1, scale: 1 }}
                  transition={{ delay: index * 0.1 }}
                >
                  <MoleculeCard
                    smiles={smilesResult}
                    title={reference || `Structure ${index + 1}`}
                    showActions={true}
                    size="md"
                  />
                </motion.div>
              ))}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
};

export default OCRView;
