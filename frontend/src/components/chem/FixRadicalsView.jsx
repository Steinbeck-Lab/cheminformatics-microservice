// Description: This component allows users to input a SMILES string containing radicals and fix them using CDK. It displays the original structure, fixed structure, and statistics about the radicals.
import React, { useState, useEffect } from "react";
import {
  HiOutlineBeaker,
  HiOutlineInformationCircle,
  HiOutlineExclamationCircle,
  HiOutlineCheckCircle,
  HiOutlineClipboard,
  HiOutlineCheck,
} from "react-icons/hi";
import SMILESInput from "../common/SMILESInput";
import LoadingScreen from "../common/LoadingScreen";
import { fixRadicals } from "../../services/chemService";
import { generate2DDepictionEnhanced, generate2DDepiction } from "../../services/depictService";

const FixRadicalsView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null);
  const [copied, setCopied] = useState(false);
  const [originalDepiction, setOriginalDepiction] = useState(null);
  const [fixedDepiction, setFixedDepiction] = useState(null);
  const [depictionLoading, setDepictionLoading] = useState(false);

  // Generate depictions when result is available
  useEffect(() => {
    const generateDepictions = async () => {
      if (!result) return;

      setDepictionLoading(true);
      setOriginalDepiction(null);
      setFixedDepiction(null);
      
      try {
        // Use enhanced depiction for original molecule with radicals (perceive_radicals enabled)
        const originalSvg = await generate2DDepictionEnhanced(smiles, {
          width: 400,
          height: 400,
          perceive_radicals: true,
          annotate: "none",
          style: "cow",
        });
        
        // Remove namespace prefix to fix rendering
        const cleanedOriginalSvg = originalSvg.replace(/ns0:/g, '').replace('xmlns:ns0=', 'xmlns=');
        setOriginalDepiction(cleanedOriginalSvg);

        // Use standard depiction for fixed molecule
        const fixedSvg = await generate2DDepiction(result.fixed_smiles, {
          width: 400,
          height: 400,
          toolkit: "cdk",
        });
        
        // Remove namespace prefix to fix rendering
        const cleanedFixedSvg = fixedSvg.replace(/ns0:/g, '').replace('xmlns:ns0=', 'xmlns=');
        setFixedDepiction(cleanedFixedSvg);
      } catch (err) {
        console.error("Error generating depictions:", err);
        // Don't set error state, just log it - we can still show SMILES
      } finally {
        setDepictionLoading(false);
      }
    };

    generateDepictions();
  }, [result, smiles]);

  // Handle copy SMILES
  const handleCopySmiles = async () => {
    if (!result?.fixed_smiles) return;

    try {
      await navigator.clipboard.writeText(result.fixed_smiles);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch (err) {
      console.error("Failed to copy SMILES:", err);
      // Try fallback method
      try {
        const textArea = document.createElement("textarea");
        textArea.value = result.fixed_smiles;
        textArea.style.position = "fixed";
        textArea.style.left = "-999999px";
        document.body.appendChild(textArea);
        textArea.select();
        document.execCommand("copy");
        document.body.removeChild(textArea);
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      } catch (fallbackErr) {
        console.error("Fallback copy also failed:", fallbackErr);
      }
    }
  };

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      setResult(null);
      setOriginalDepiction(null);
      setFixedDepiction(null);
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);
    setOriginalDepiction(null);
    setFixedDepiction(null);

    try {
      const data = await fixRadicals(trimmedSmiles);
      setResult(data);

      // Show info message if no radicals were found
      if (data.radicals_detected === 0) {
        setError("No radicals detected in the input molecule.");
      }
    } catch (err) {
      console.error("Radical fixing error:", err);
      setError(`Error fixing radicals: ${err.message || "An unknown error occurred."}`);
      setResult(null);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Fix Radicals in Molecules
        </h2>

        {/* Info Box */}
        <div className="mb-4 p-4 bg-blue-50 dark:bg-blue-900 dark:bg-opacity-30 rounded-md border border-blue-200 dark:border-blue-700">
          <div className="flex items-start">
            <HiOutlineInformationCircle className="h-5 w-5 text-blue-500 dark:text-blue-400 mr-3 mt-0.5 flex-shrink-0" />
            <div className="text-sm text-blue-700 dark:text-blue-200">
              <p className="font-medium mb-1">About Radical Fixing</p>
              <p>
                This tool detects and fixes radical electrons (unpaired electrons) on atoms in molecular structures.
                It supports radicals on Carbon (C), Nitrogen (N), and Oxygen (O) atoms.
              </p>
              <p className="mt-2">
                <strong>Examples:</strong> [CH3] (methyl radical), C[CH2] (ethyl radical), [OH] (hydroxyl radical)
              </p>
            </div>
          </div>
        </div>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            required
            placeholder="Enter SMILES with radicals (e.g., [CH3], C[CH2], [OH])"
          />

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              disabled={!smiles.trim() || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !smiles.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm"
              }`}
            >
              <HiOutlineBeaker className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? "Fixing Radicals..." : "Fix Radicals"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Fixing radicals..." />}

      {/* Error/Info Display */}
      {error && !loading && (
        <div
          className={`p-4 rounded-md flex items-start shadow ${
            error.startsWith("No radicals")
              ? "bg-blue-50 dark:bg-blue-900 dark:bg-opacity-30 text-blue-700 dark:text-blue-200 border border-blue-300 dark:border-blue-700"
              : "bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700"
          }`}
          role={error.startsWith("No radicals") ? "status" : "alert"}
        >
          {error.startsWith("No radicals") ? (
            <HiOutlineInformationCircle
              className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-blue-500 dark:text-blue-400"
              aria-hidden="true"
            />
          ) : (
            <HiOutlineExclamationCircle
              className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
              aria-hidden="true"
            />
          )}
          <span>{error}</span>
        </div>
      )}

      {/* Results Display Section */}
      {result && !loading && (
        <div className="space-y-4">
          {/* Statistics Card */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-4 border-b border-gray-200 dark:border-gray-700 pb-3">
              Radical Fixing Statistics
            </h3>

            <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
              {/* Radicals Detected */}
              <div className="bg-gradient-to-br from-blue-50 to-blue-100 dark:from-blue-900 dark:to-blue-800 p-4 rounded-lg border border-blue-200 dark:border-blue-700">
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm font-medium text-blue-600 dark:text-blue-300">
                      Radicals Detected
                    </p>
                    <p className="text-2xl font-bold text-blue-700 dark:text-blue-200 mt-1">
                      {result.radicals_detected}
                    </p>
                  </div>
                  <HiOutlineInformationCircle className="h-8 w-8 text-blue-400 dark:text-blue-500" />
                </div>
              </div>

              {/* Radicals Fixed */}
              <div className="bg-gradient-to-br from-green-50 to-green-100 dark:from-green-900 dark:to-green-800 p-4 rounded-lg border border-green-200 dark:border-green-700">
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm font-medium text-green-600 dark:text-green-300">
                      Radicals Fixed
                    </p>
                    <p className="text-2xl font-bold text-green-700 dark:text-green-200 mt-1">
                      {result.radicals_fixed}
                    </p>
                  </div>
                  <HiOutlineCheckCircle className="h-8 w-8 text-green-400 dark:text-green-500" />
                </div>
              </div>

              {/* Success Rate */}
              <div className="bg-gradient-to-br from-purple-50 to-purple-100 dark:from-purple-900 dark:to-purple-800 p-4 rounded-lg border border-purple-200 dark:border-purple-700">
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm font-medium text-purple-600 dark:text-purple-300">
                      Success Rate
                    </p>
                    <p className="text-2xl font-bold text-purple-700 dark:text-purple-200 mt-1">
                      {result.radicals_detected > 0
                        ? Math.round((result.radicals_fixed / result.radicals_detected) * 100)
                        : 100}
                      %
                    </p>
                  </div>
                  <HiOutlineBeaker className="h-8 w-8 text-purple-400 dark:text-purple-500" />
                </div>
              </div>
            </div>

            {/* Note about unsupported radicals */}
            {result.radicals_detected > result.radicals_fixed && (
              <div className="mt-4 p-3 bg-yellow-50 dark:bg-yellow-900 dark:bg-opacity-30 rounded-md border border-yellow-200 dark:border-yellow-700">
                <div className="flex items-start">
                  <HiOutlineInformationCircle className="h-5 w-5 text-yellow-600 dark:text-yellow-400 mr-2 mt-0.5 flex-shrink-0" />
                  <p className="text-sm text-yellow-700 dark:text-yellow-200">
                    Note: Some radicals were detected but not fixed. Currently, only radicals on Carbon (C),
                    Nitrogen (N), and Oxygen (O) atoms are supported.
                  </p>
                </div>
              </div>
            )}
          </div>

          {/* Structures Comparison Card */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
            <div className="flex items-center justify-between mb-4 border-b border-gray-200 dark:border-gray-700 pb-3">
              <h3 className="text-lg font-semibold text-gray-900 dark:text-white">
                Structure Comparison
              </h3>
              
              {/* Copy Fixed SMILES Button */}
              <button
                onClick={handleCopySmiles}
                disabled={!result?.fixed_smiles}
                className={`flex items-center px-4 py-2 rounded-lg font-medium text-sm transition-all duration-200 ${
                  copied
                    ? "bg-green-500 text-white"
                    : "bg-blue-600 hover:bg-blue-700 text-white dark:bg-blue-500 dark:hover:bg-blue-600"
                } disabled:opacity-50 disabled:cursor-not-allowed`}
                title="Copy fixed SMILES to clipboard"
              >
                {copied ? (
                  <>
                    <HiOutlineCheck className="h-4 w-4 mr-2" />
                    Copied!
                  </>
                ) : (
                  <>
                    <HiOutlineClipboard className="h-4 w-4 mr-2" />
                    Copy Fixed SMILES
                  </>
                )}
              </button>
            </div>

            {depictionLoading && (
              <div className="flex justify-center items-center py-8">
                <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500"></div>
                <span className="ml-3 text-gray-600 dark:text-gray-400">Generating depictions...</span>
              </div>
            )}

            {!depictionLoading && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {/* Original Structure with Radicals */}
                <div>
                  <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-3 flex items-center">
                    <span className="inline-block w-2 h-2 bg-red-500 rounded-full mr-2"></span>
                    Original Structure (with radicals)
                  </h4>
                  
                  {/* Depiction Container */}
                  <div className="bg-white dark:bg-gray-900 rounded-lg border-2 border-gray-200 dark:border-gray-700 p-4 mb-3">
                    {originalDepiction ? (
                      <div 
                        className="w-full h-80 flex items-center justify-center [&>svg]:max-w-full [&>svg]:max-h-full [&>svg]:w-auto [&>svg]:h-auto"
                        dangerouslySetInnerHTML={{ __html: originalDepiction }}
                      />
                    ) : (
                      <div className="w-full h-80 flex items-center justify-center bg-gray-100 dark:bg-gray-800 rounded">
                        <p className="text-gray-500 dark:text-gray-400 text-sm">
                          Depiction unavailable
                        </p>
                      </div>
                    )}
                  </div>

                  {/* SMILES Display */}
                  <div className="bg-gray-100 dark:bg-gray-900 p-3 rounded-md border border-gray-200 dark:border-gray-700 font-mono text-xs overflow-x-auto shadow-sm">
                    <div className="flex items-start justify-between">
                      <pre className="text-gray-700 dark:text-gray-300 whitespace-pre-wrap break-all flex-1">
                        {smiles}
                      </pre>
                    </div>
                  </div>
                </div>

                {/* Fixed Structure */}
                <div>
                  <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-3 flex items-center">
                    <span className="inline-block w-2 h-2 bg-green-500 rounded-full mr-2"></span>
                    Fixed Structure
                  </h4>
                  
                  {/* Depiction Container */}
                  <div className="bg-white dark:bg-gray-900 rounded-lg border-2 border-green-200 dark:border-green-700 p-4 mb-3">
                    {fixedDepiction ? (
                      <div 
                        className="w-full h-80 flex items-center justify-center [&>svg]:max-w-full [&>svg]:max-h-full [&>svg]:w-auto [&>svg]:h-auto"
                        dangerouslySetInnerHTML={{ __html: fixedDepiction }}
                      />
                    ) : (
                      <div className="w-full h-80 flex items-center justify-center bg-gray-100 dark:bg-gray-800 rounded">
                        <p className="text-gray-500 dark:text-gray-400 text-sm">
                          Depiction unavailable
                        </p>
                      </div>
                    )}
                  </div>

                  {/* SMILES Display */}
                  <div className="bg-gray-100 dark:bg-gray-900 p-3 rounded-md border border-gray-200 dark:border-gray-700 font-mono text-xs overflow-x-auto shadow-sm">
                    <div className="flex items-start justify-between">
                      <pre className="text-gray-700 dark:text-gray-300 whitespace-pre-wrap break-all flex-1">
                        {result.fixed_smiles}
                      </pre>
                    </div>
                  </div>
                </div>
              </div>
            )}

            {/* Comparison Note */}
            {smiles === result.fixed_smiles && result.radicals_detected === 0 && (
              <div className="mt-4 p-3 bg-gray-50 dark:bg-gray-700 rounded-md border border-gray-200 dark:border-gray-600">
                <p className="text-sm text-gray-600 dark:text-gray-300 text-center">
                  The structures are identical - no radicals were detected in the input.
                </p>
              </div>
            )}

            {smiles !== result.fixed_smiles && result.radicals_fixed > 0 && (
              <div className="mt-4 p-3 bg-green-50 dark:bg-green-900 dark:bg-opacity-30 rounded-md border border-green-200 dark:border-green-700">
                <div className="flex items-start">
                  <HiOutlineCheckCircle className="h-5 w-5 text-green-600 dark:text-green-400 mr-2 mt-0.5 flex-shrink-0" />
                  <p className="text-sm text-green-700 dark:text-green-200">
                    Successfully fixed {result.radicals_fixed} radical{result.radicals_fixed !== 1 ? "s" : ""} by
                    adding implicit hydrogens.
                  </p>
                </div>
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
};

export default FixRadicalsView;
