// Description : This component allows users to calculate molecular descriptors from a SMILES string.
import React, { useState } from "react";
// Ensure all used icons are imported
import {
  HiOutlineCalculator,
  HiOutlineInformationCircle,
  HiOutlineExclamationCircle, // Added for error display
} from "react-icons/hi";
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import LoadingScreen from "../common/LoadingScreen";
// Assuming this service is configured correctly
import { calculateNPLikeness } from "../../services/chemService"; // Assuming this service exists

const NPlikenessView = () => {
  const [smiles, setSmiles] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [npScore, setNpScore] = useState(null); // Store the score (null initially)

  // Handle form submission
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string.");
      setNpScore(null); // Clear previous score
      return;
    }

    setLoading(true);
    setError(null);
    setNpScore(null); // Clear previous score before calculation

    try {
      // Call the service function
      const score = await calculateNPLikeness(trimmedSmiles);
      // Ensure score is a number before setting state
      if (typeof score === "number" && !isNaN(score)) {
        setNpScore(score);
      } else {
        // Handle cases where API might not return a valid number
        console.warn("Received non-numeric score:", score);
        throw new Error("Received an invalid score format from the server.");
      }
    } catch (err) {
      console.error("NP-likeness calculation error:", err); // Log the error
      setError(`Error calculating NP-likeness: ${err.message || "An unknown error occurred."}`);
      setNpScore(null); // Ensure score is null on error
    } finally {
      setLoading(false);
    }
  };

  // Determine the Tailwind CSS color class based on the NP-likeness score
  // Adjusted colors for better light/dark mode contrast
  const getNpScoreColor = (score) => {
    if (score === null || isNaN(score)) return "text-gray-500 dark:text-gray-400"; // Default/invalid color
    if (score >= 1) return "text-green-600 dark:text-green-400";
    if (score >= 0) return "text-blue-600 dark:text-blue-400";
    if (score >= -1) return "text-amber-600 dark:text-yellow-400"; // Using amber for light mode yellow contrast
    if (score >= -2) return "text-orange-600 dark:text-orange-400";
    return "text-red-600 dark:text-red-400"; // Score < -2
  };

  // Get a textual interpretation of the score
  const getNpScoreInterpretation = (score) => {
    if (score === null || isNaN(score)) return "Enter a SMILES string to calculate the score.";
    if (score >= 1) return "Very high natural product likeness - likely to be a natural product";
    if (score >= 0) return "High natural product likeness";
    if (score >= -1) return "Moderate natural product likeness";
    if (score >= -2) return "Low natural product likeness";
    return "Very low natural product likeness - likely to be synthetic";
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Natural Product Likeness Score
        </h2>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          {/* SMILES Input */}
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            required
            // Pass theme props if needed
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
              <HiOutlineCalculator className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? "Calculating..." : "Calculate NP Score"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Calculating NP-likeness score..." />}

      {/* Error Display */}
      {error &&
        !loading && ( // Show error only if not loading
          <div
            className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow"
            role="alert"
          >
            <HiOutlineExclamationCircle
              className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
              aria-hidden="true"
            />
            <span>{error}</span>
          </div>
        )}

      {/* Results Display Section */}
      {/* Show only if score exists (is a number) and not loading and no error */}
      {typeof npScore === "number" && !isNaN(npScore) && !loading && !error && (
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
          {/* Results Header */}
          <div className="flex items-center justify-between mb-4 border-b border-gray-200 dark:border-gray-700 pb-3">
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white">
              NP-likeness Results
            </h3>
            {/* Info Button (Consider adding tooltip/modal functionality) */}
            <button
              className="p-1 text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white rounded-full focus:outline-none focus:ring-2 focus:ring-blue-500"
              title="About NP-likeness scores (Not implemented)"
              onClick={() =>
                alert(
                  "Info about NP-likeness score ranges:\n\n> 1: Very High NP-like\n0 to 1: High NP-like\n-1 to 0: Moderate NP-like\n-2 to -1: Low NP-like\n< -2: Very Low NP-like (Synthetic)"
                )
              } // Simple alert for now
            >
              <HiOutlineInformationCircle className="h-6 w-6" aria-hidden="true" />
            </button>
          </div>

          {/* Results Content */}
          <div className="space-y-8">
            {" "}
            {/* Increased spacing */}
            {/* Score Display */}
            <div className="text-center">
              <div className={`text-5xl font-bold mb-2 ${getNpScoreColor(npScore)}`}>
                {/* Format score to 2 decimal places */}
                {npScore.toFixed(2)}
              </div>
              {/* Interpretation Text */}
              <p className="text-sm text-gray-600 dark:text-gray-400">
                {getNpScoreInterpretation(npScore)}
              </p>
            </div>
            {/* Score Gauge Visualization */}
            <div className="relative pt-6 px-2">
              {" "}
              {/* Added padding */}
              {/* Gauge Background Bar */}
              <div className="overflow-hidden h-2.5 mb-4 flex rounded-full bg-gray-200 dark:bg-gray-700 shadow-inner">
                {/* Gradient Segments */}
                {/* These gradients provide a visual representation */}
                <div className="h-full w-[30%] bg-gradient-to-r from-red-500 to-orange-500"></div>{" "}
                {/* -5 to -2 */}
                <div className="h-full w-[10%] bg-gradient-to-r from-orange-500 to-amber-500"></div>{" "}
                {/* -2 to -1 */}
                <div className="h-full w-[10%] bg-gradient-to-r from-amber-500 to-blue-500"></div>{" "}
                {/* -1 to 0 */}
                <div className="h-full w-[10%] bg-gradient-to-r from-blue-500 to-green-400"></div>{" "}
                {/* 0 to 1 */}
                <div className="h-full w-[40%] bg-gradient-to-r from-green-400 to-green-600"></div>{" "}
                {/* 1 to 5 */}
              </div>
              {/* Marker for the current score */}
              {/* Calculate position based on -5 to 5 range (total range of 10) */}
              <div
                className="absolute bottom-4 w-4 h-4 bg-white dark:bg-gray-200 rounded-full ring-2 ring-gray-500 dark:ring-gray-400 transform -translate-x-1/2 shadow-lg transition-all duration-500 ease-out"
                style={{
                  // Normalize score to 0-1 range, then scale to 100%
                  left: `${Math.max(0, Math.min(100, ((npScore + 5) / 10) * 100))}%`,
                }}
                title={`Score: ${npScore.toFixed(2)}`} // Tooltip for exact score
              ></div>
              {/* Labels for Gauge */}
              {/* Positioned below the gauge */}
              <div className="flex justify-between text-xs mt-1">
                <span className="text-red-600 dark:text-red-400">-5</span>
                <span className="text-orange-600 dark:text-orange-400">-2</span>
                <span className="text-amber-600 dark:text-yellow-400">-1</span>
                <span className="text-blue-600 dark:text-blue-400">0</span>
                <span className="text-green-600 dark:text-green-400">1</span>
                <span className="text-green-700 dark:text-green-300">5</span>
              </div>
            </div>
            {/* Informational Box about the score */}
            <div className="bg-gray-50 dark:bg-gray-900 p-4 rounded-lg text-sm text-gray-700 dark:text-gray-300 border border-gray-200 dark:border-gray-700 shadow-sm">
              <p>
                <strong>About NP-likeness Score:</strong> This score measures molecular similarity
                to known natural products, ranging from -5 (synthetic-like) to +5 (natural
                product-like).
              </p>
              <p className="mt-2 text-xs text-gray-500 dark:text-gray-400">
                Calculation based on RDKit implementation of the method by Ertl et al., J.
                Cheminform. 1:8 (2008).
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Initial State Message (Optional) */}
      {/* Show only if no score, not loading, and no error */}
      {npScore === null && !loading && !error && (
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 text-center shadow">
          <p className="text-gray-600 dark:text-gray-300">
            Enter a SMILES string to calculate its Natural Product-likeness score.
          </p>
        </div>
      )}
    </div>
  );
};

export default NPlikenessView;
