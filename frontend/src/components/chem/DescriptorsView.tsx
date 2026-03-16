// Description: This component allows users to input a SMILES string and calculate molecular descriptors using different toolkits (RDKit, CDK) and formats (JSON, HTML). It handles loading states, errors, and displays results in a user-friendly manner.
import { useState } from "react";
// Ensure all used icons are imported
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import { getDescriptors } from "../../services/chemService"; // Assuming this service exists
import DOMPurify from "dompurify"; // Import DOMPurify for sanitizing HTML
import { AlertCircle, Calculator, FileBarChart, Loader2 } from "lucide-react";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { EmptyState } from "@/components/feedback/EmptyState";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { AddToCompareButton } from "../common/AddToCompareButton";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";

const DescriptorsView = () => {
  const [smiles, setSmiles] = useState("");
  const [descriptors, setDescriptors] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [toolkit, setToolkit] = useState("rdkit"); // Default toolkit
  const [format, setFormat] = useState("json"); // Default format

  // Handle calculation request
  const handleCalculateDescriptors = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim();
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string."); // Set error if input is empty
      return;
    }

    setLoading(true);
    setError(null);
    setDescriptors(null); // Clear previous results

    try {
      // Call the service function
      const result = await getDescriptors(trimmedSmiles, toolkit, format);
      setDescriptors(result); // Store the result (could be JSON object or HTML string)
    } catch (err) {
      console.error("Descriptor calculation error:", err); // Log the actual error
      // Set a user-friendly error message
      setError(getErrorMessage("chem", err));
    } finally {
      setLoading(false);
    }
  };

  // Render the results, either as a table (JSON) or raw HTML
  const renderDescriptorsTable = () => {
    if (!descriptors) return null;

    // Handle HTML format response
    if (format === "html" && typeof descriptors === "string") {
      // Apply base styling for the container, actual table style comes from API
      // Import DOMPurify at the top of your file
      // import DOMPurify from 'dompurify';

      return (
        <div
          className="prose prose-sm dark:prose-invert max-w-none bg-white dark:bg-gray-800 p-4 rounded-sm border border-gray-200 dark:border-gray-700 shadow-xs"
          // Sanitize HTML before rendering to prevent XSS attacks
          dangerouslySetInnerHTML={{ __html: DOMPurify.sanitize(descriptors) }}
        />
      );
    }

    // Handle JSON format response (assuming 'descriptors' is an object)
    if (format === "json" && typeof descriptors === "object" && descriptors !== null) {
      const entries = Object.entries(descriptors);
      if (entries.length === 0) {
        return <p className="text-gray-500 dark:text-gray-400">No descriptors returned.</p>;
      }

      return (
        // Table container with overflow and styling
        <div className="overflow-x-auto border border-gray-200 dark:border-gray-700 rounded-md shadow-xs">
          <table className="w-full min-w-[400px] border-collapse">
            {/* Table Header */}
            <thead className="bg-gray-50 dark:bg-gray-700">
              <tr>
                <th
                  scope="col"
                  className="px-4 py-2 text-left text-xs font-medium uppercase tracking-wider text-gray-500 dark:text-gray-400"
                >
                  Property
                </th>
                <th
                  scope="col"
                  className="px-4 py-2 text-right text-xs font-medium uppercase tracking-wider text-gray-500 dark:text-gray-400"
                >
                  Value
                </th>
              </tr>
            </thead>
            {/* Table Body */}
            <tbody className="bg-white dark:bg-gray-800 divide-y divide-gray-200 dark:divide-gray-700">
              {entries.map(([key, value]) => (
                <tr
                  key={key}
                  className="hover:bg-gray-50 dark:hover:bg-gray-700/50 transition-colors duration-150"
                >
                  {/* Property Name */}
                  <td className="px-4 py-2 text-sm text-gray-700 dark:text-gray-300 font-medium whitespace-nowrap">
                    {key}
                  </td>
                  {/* Property Value */}
                  <td className="px-4 py-2 text-right text-sm text-gray-800 dark:text-blue-300 font-mono">
                    {/* Handle complex values (objects/arrays) by stringifying */}
                    {typeof value === "object" && value !== null
                      ? JSON.stringify(value)
                      : String(value)}{" "}
                    {/* Convert other types to string */}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      );
    }

    // Fallback if format is unexpected or data is not in expected type
    return (
      <p className="text-red-600 dark:text-red-400">
        Could not render descriptors. Unexpected data format received.
      </p>
    );
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 rounded-lg shadow-md dark:shadow-lg p-6">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Molecular Descriptors Calculator
        </h2>

        {/* Form */}
        <form onSubmit={handleCalculateDescriptors} className="space-y-4">
          {/* SMILES Input */}
          <SMILESInput
            value={smiles}
            onChange={setSmiles}
            required
            // Pass theme props if needed
          />

          {/* Options Grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 pt-2">
            {/* Toolkit Selection */}
            <div>
              <label
                htmlFor="toolkit-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Toolkit
              </label>
              <Select value={toolkit} onValueChange={setToolkit}>
                <SelectTrigger id="toolkit-select" className="w-full">
                  <SelectValue placeholder="Select toolkit" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="rdkit">RDKit</SelectItem>
                  <SelectItem value="cdk">CDK</SelectItem>
                  <SelectItem value="all">All (CDK + RDKit)</SelectItem>
                </SelectContent>
              </Select>
            </div>

            {/* Format Selection */}
            <div>
              <label
                htmlFor="format-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Output Format
              </label>
              <Select value={format} onValueChange={setFormat}>
                <SelectTrigger id="format-select" className="w-full">
                  <SelectValue placeholder="Select format" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="json">JSON (Table)</SelectItem>
                  <SelectItem value="html">HTML (Raw)</SelectItem>
                </SelectContent>
              </Select>
            </div>
          </div>

          {/* Submit Button */}
          <div className="pt-2">
            <Button
              type="submit"
              disabled={!smiles.trim() || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-hidden focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !smiles.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-xs"
              }`}
            >
              {loading ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Calculating...
                </>
              ) : (
                <>
                  <Calculator className="mr-2 h-5 w-5" aria-hidden="true" />
                  Calculate Descriptors
                </>
              )}
            </Button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && !descriptors && <ToolSkeleton variant="descriptors" />}

      {/* Error Display */}
      {error && !loading && (
        <GlassErrorCard
          message={error}
          onRetry={() => {
            setError(null);
            document.getElementById("smiles-input")?.focus();
          }}
        />
      )}

      {/* Results Display Section */}
      {/* Show only if descriptors exist and not loading */}
      {descriptors && !loading && (
        <div className="bg-white dark:bg-gray-800 rounded-lg shadow-md dark:shadow-lg overflow-hidden">
          {/* Results Header */}
          <div className="bg-gray-50 dark:bg-gray-900 px-6 py-4 flex justify-between items-center border-b border-gray-200 dark:border-gray-700">
            <h3 className="text-lg font-medium text-gray-900 dark:text-white flex items-center">
              <FileBarChart
                className="mr-3 h-5 w-5 text-blue-600 dark:text-blue-400 shrink-0"
                aria-hidden="true"
              />
              Molecular Descriptors
            </h3>
            <div className="flex items-center gap-2">
              {smiles.trim() && (
                <AddToCompareButton
                  smiles={smiles.trim()}
                  title="Descriptors molecule"
                  descriptors={
                    format === "json" && typeof descriptors === "object" && descriptors !== null
                      ? descriptors
                      : undefined
                  }
                  sourceToolId="descriptors"
                />
              )}
              {/* Toolkit Info */}
              <div className="text-sm text-gray-500 dark:text-gray-400 bg-gray-100 dark:bg-gray-700 px-2 py-0.5 rounded-sm">
                Toolkit: {toolkit === "all" ? "RDKit + CDK" : toolkit.toUpperCase()} | Format:{" "}
                {format.toUpperCase()}
              </div>
            </div>
          </div>

          {/* Results Content Area */}
          <div className="p-6">{renderDescriptorsTable()}</div>
        </div>
      )}

      {/* Initial State Message (Optional) */}
      {!descriptors && !loading && !error && (
        <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 text-center shadow-sm">
          <p className="text-gray-600 dark:text-gray-300">
            Enter a SMILES string and select options to calculate molecular descriptors.
          </p>
        </div>
      )}
    </div>
  );
};
export default DescriptorsView;
