// Description: Tanimoto similarity calculator component
import React, { useState } from "react";
// Ensure all used icons are imported
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { useAppContext } from "../../context/AppContext"; // Assuming context provides addRecentMolecule
import api from "../../services/api"; // Assuming api service is configured
import { AlertCircle, Calculator, Info, Loader2 } from "lucide-react";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { EmptyState } from "@/components/feedback/EmptyState";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";

const TanimotoView = () => {
  const [smiles1, setSmiles1] = useState("");
  const [smiles2, setSmiles2] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null); // Stores the similarity score (number or null)
  const [fingerprinter, setFingerprinter] = useState("ECFP"); // Default fingerprinter
  const [toolkit, setToolkit] = useState("rdkit"); // Default toolkit
  const [nBits, setNBits] = useState(2048); // Default nBits
  const [radius, setRadius] = useState(2); // Default radius (for ECFP/MAPC)
  const { addRecentMolecule } = useAppContext();
  const [infoVisible, setInfoVisible] = useState(false); // State for toggling info box

  // Function to calculate Tanimoto similarity via API
  const calculateTanimotoSimilarity = async () => {
    const trimmedSmiles1 = smiles1.trim();
    const trimmedSmiles2 = smiles2.trim();

    if (!trimmedSmiles1 || !trimmedSmiles2) {
      throw new Error("Please provide both SMILES strings");
    }

    // Format the combined SMILES string correctly for the API parameter
    const combinedSmiles = `${trimmedSmiles1},${trimmedSmiles2}`;

    // Make the API call directly
    const response = await api.get(`/chem/tanimoto`, {
      params: {
        smiles: combinedSmiles,
        toolkit,
        fingerprinter,
        nBits,
        radius,
      },
    });

    // Assuming response.data is the similarity score (number)
    if (typeof response.data !== "number") {
      throw new Error("Received invalid similarity score from API.");
    }
    return response.data;
  };

  // Function to handle the calculate button click
  const handleCalculate = async (e) => {
    // Prevent default form submission if called from a form event
    if (e && typeof e.preventDefault === "function") {
      e.preventDefault();
    }

    const trimmedSmiles1 = smiles1.trim();
    const trimmedSmiles2 = smiles2.trim();

    if (!trimmedSmiles1 || !trimmedSmiles2) {
      setError("Please provide both SMILES strings");
      setResult(null); // Clear previous result
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null); // Clear previous result before calculation

    try {
      const similarity = await calculateTanimotoSimilarity();
      setResult(similarity);

      // Add molecules to recent list on successful calculation
      addRecentMolecule({
        smiles: trimmedSmiles1,
        name: "Molecule 1 (Tanimoto)", // Add context
        timestamp: new Date().toISOString(),
      });
      addRecentMolecule({
        smiles: trimmedSmiles2,
        name: "Molecule 2 (Tanimoto)", // Add context
        timestamp: new Date().toISOString(),
      });
    } catch (err) {
      console.error("Tanimoto calculation error:", err);
      setError(getErrorMessage("chem", err));
      setResult(null); // Ensure result is null on error
    } finally {
      setLoading(false);
    }
  };

  // Function to get Tailwind color class based on similarity value
  // Adjusted colors for light/dark mode contrast
  const getSimilarityColor = (value) => {
    if (value === null || isNaN(value)) return "text-gray-500 dark:text-gray-400";
    if (value >= 0.85) return "text-green-600 dark:text-green-400"; // High
    if (value >= 0.7) return "text-lime-600 dark:text-lime-400"; // Good (using lime for differentiation)
    if (value >= 0.4) return "text-amber-600 dark:text-yellow-400"; // Moderate
    return "text-red-600 dark:text-red-400"; // Low
  };

  // Function to get background color class for the gauge bar
  const getGaugeBgColor = (value) => {
    if (value === null || isNaN(value)) return "bg-gray-300 dark:bg-gray-600";
    if (value >= 0.85) return "bg-green-500 dark:bg-green-500";
    if (value >= 0.7) return "bg-lime-500 dark:bg-lime-500";
    if (value >= 0.4) return "bg-yellow-400 dark:bg-yellow-400"; // Yellow is often visible enough
    return "bg-red-500 dark:bg-red-500";
  };

  // Function to get textual description based on similarity value
  const getSimilarityDescription = (value) => {
    if (value === null || isNaN(value)) return "Enter SMILES and calculate similarity.";
    if (value >= 0.85) return "High similarity - very close structural match";
    if (value >= 0.7) return "Good similarity - significant structural overlap";
    if (value >= 0.4) return "Moderate similarity - some shared features";
    return "Low similarity - structurally different compounds";
  };

  // Get appropriate fingerprinter options based on selected toolkit
  const getFingerprinterOptions = () => {
    if (toolkit === "rdkit") {
      return [
        { value: "ECFP", label: "ECFP (Extended Connectivity)" },
        { value: "RDKit", label: "RDKit Topological" },
        { value: "Atompairs", label: "Atom Pairs" },
        { value: "MACCS", label: "MACCS Keys" },
        { value: "MAPC", label: "MAPC (MinHashed Atom-Pair Chiral)" },
      ];
    } else {
      // cdk toolkit
      return [
        { value: "ECFP", label: "ECFP (Extended Connectivity)" },
        { value: "PubChem", label: "PubChem Fingerprints" },
      ];
    }
  };

  // Handle toolkit change - reset fingerprinter if needed
  const handleToolkitChange = (newToolkit) => {
    setToolkit(newToolkit);
    // No need for currentOptions variable
    const newOptions = getFingerprinterOptions(newToolkit); // Get options for the *new* toolkit

    // Check if the current fingerprinter is valid in the new toolkit's options
    const isValid = newOptions.some((option) => option.value === fingerprinter);

    // If not valid, reset to the default (e.g., ECFP) or the first available option
    if (!isValid) {
      setFingerprinter(newOptions[0]?.value || "ECFP"); // Reset to first option or default ECFP
    }
  };

  return (
    // Main container
    <div className="space-y-6 p-4 md:p-6">
      {/* Header and Info Toggle */}
      <div className="flex justify-between items-start mb-4">
        <h2 className="text-2xl font-bold text-gray-800 dark:text-blue-400">
          Tanimoto Similarity Calculator
        </h2>
        <Button
          variant="ghost"
          onClick={() => setInfoVisible(!infoVisible)}
          className="p-2 text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white rounded-full focus:outline-hidden focus:ring-2 focus:ring-blue-500"
          aria-label="Toggle Information"
          aria-expanded={infoVisible}
        >
          <Info className="h-6 w-6" />
        </Button>
      </div>

      {/* Toggleable Info Box */}
      {infoVisible && (
        // Info box styling for light/dark mode
        <div className="bg-blue-50 dark:bg-gray-800 border border-blue-200 dark:border-gray-700 rounded-lg p-4 mb-6 shadow-xs">
          <h3 className="text-lg font-semibold text-blue-800 dark:text-blue-300 mb-2">
            About Tanimoto Similarity
          </h3>
          <p className="text-gray-700 dark:text-gray-300 mb-2">
            The Tanimoto coefficient (Tc), ranging from 0 to 1, measures structural similarity
            between molecules based on their fingerprints. 1 means identical (based on the
            fingerprint), 0 means no common features.
          </p>
          <div className="mt-3 text-sm">
            <p className="mb-1 font-medium text-gray-600 dark:text-gray-400">
              Interpretation guide:
            </p>
            {/* Interpretation list styling */}
            <ul className="list-disc pl-5 space-y-1 text-gray-600 dark:text-gray-300">
              {/* Adjusted colors for light/dark */}
              <li>
                <span className="font-medium text-green-600 dark:text-green-400">≥ 0.85:</span> High
                similarity
              </li>
              <li>
                <span className="font-medium text-lime-600 dark:text-lime-400">0.7 - 0.85:</span>{" "}
                Good similarity
              </li>
              <li>
                <span className="font-medium text-amber-600 dark:text-yellow-400">0.4 - 0.7:</span>{" "}
                Moderate similarity
              </li>
              <li>
                <span className="font-medium text-red-600 dark:text-red-400">&lt; 0.4:</span> Low
                similarity
              </li>
            </ul>
          </div>
        </div>
      )}

      {/* Input Molecules Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        {/* Molecule 1 Input Card */}
        <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
          <SMILESInput
            value={smiles1}
            onChange={setSmiles1}
            label="Molecule 1 SMILES"
            required={true}
            // Pass theme props if needed
          />
        </div>
        {/* Molecule 2 Input Card */}
        <div className="bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
          <SMILESInput
            value={smiles2}
            onChange={setSmiles2}
            label="Molecule 2 SMILES"
            required={true}
            // Pass theme props if needed
          />
        </div>
      </div>

      {/* Options and Calculate Button Grid */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 items-end">
        {/* Options Column 1 Card */}
        <div className="space-y-4 bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700 lg:col-span-1">
          <div>
            <label
              htmlFor="toolkit-select-tanimoto"
              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
            >
              Toolkit
            </label>
            {/* Select styling */}
            <Select value={toolkit} onValueChange={handleToolkitChange}>
              <SelectTrigger id="toolkit-select-tanimoto" className="w-full">
                <SelectValue placeholder="Select toolkit" />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="rdkit">RDKit</SelectItem>
                <SelectItem value="cdk">CDK</SelectItem>
              </SelectContent>
            </Select>
          </div>
          <div>
            <label
              htmlFor="fingerprinter-select"
              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
            >
              Fingerprinter
            </label>
            {/* Select styling */}
            <Select value={fingerprinter} onValueChange={setFingerprinter}>
              <SelectTrigger id="fingerprinter-select" className="w-full">
                <SelectValue placeholder="Select fingerprinter" />
              </SelectTrigger>
              <SelectContent>
                {getFingerprinterOptions().map((option) => (
                  <SelectItem key={option.value} value={option.value}>
                    {option.label}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>
        </div>

        {/* Options Column 2 Card */}
        <div className="space-y-4 bg-white dark:bg-gray-800 p-4 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700 lg:col-span-1">
          <div className="grid grid-cols-2 gap-4">
            <div>
              <label
                htmlFor="nbits-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Bits
              </label>
              {/* Select styling */}
              <Select value={String(nBits)} onValueChange={(v) => setNBits(Number(v))}>
                <SelectTrigger id="nbits-select" className="w-full">
                  <SelectValue placeholder="Select bits" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="1024">1024</SelectItem>
                  <SelectItem value="2048">2048</SelectItem>
                  <SelectItem value="4096">4096</SelectItem>
                </SelectContent>
              </Select>
            </div>
            <div>
              <label
                htmlFor="radius-select"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                Radius
              </label>
              {/* Select styling with disabled state */}
              <Select
                value={String(radius)}
                onValueChange={(v) => setRadius(Number(v))}
                disabled={fingerprinter !== "ECFP" && fingerprinter !== "MAPC"}
              >
                <SelectTrigger
                  id="radius-select"
                  className={cn(
                    "w-full",
                    fingerprinter !== "ECFP" && fingerprinter !== "MAPC" && "opacity-50"
                  )}
                >
                  <SelectValue placeholder="Select radius" />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="2">2 (ECFP4)</SelectItem>
                  <SelectItem value="3">3 (ECFP6)</SelectItem>
                  <SelectItem value="4">4 (ECFP8)</SelectItem>
                </SelectContent>
              </Select>
              {/* Info text for disabled state */}
              {fingerprinter !== "ECFP" && fingerprinter !== "MAPC" && (
                <p className="text-xs text-gray-400 dark:text-gray-500 mt-1">
                  Radius only applies to ECFP/MAPC.
                </p>
              )}
            </div>
          </div>
        </div>

        {/* Calculate Button Column */}
        <div className="lg:col-span-1">
          {/* Button styling */}
          <Button
            onClick={handleCalculate} // Use onClick for button
            disabled={loading || !smiles1.trim() || !smiles2.trim()}
            className={`w-full px-6 py-3 rounded-lg font-medium flex items-center justify-center transition-colors duration-200 focus:outline-hidden focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
              loading || !smiles1.trim() || !smiles2.trim()
                ? "bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed" // Disabled state
                : "bg-blue-600 text-white hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-xs" // Enabled state
            }`}
          >
            {loading ? (
              <>
                <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                Calculating...
              </>
            ) : (
              <>
                <Calculator className="mr-2 h-5 w-5" />
                Calculate Similarity
              </>
            )}
          </Button>
        </div>
      </div>

      {/* Loading State */}
      {loading && !result && <ToolSkeleton variant="general" />}

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
      {typeof result === "number" && !isNaN(result) && !loading && (
        <div className="space-y-6 pt-6">
          {/* Score and Gauge Card */}
          <div className="bg-white dark:bg-gray-800 rounded-lg p-6 border border-gray-200 dark:border-gray-700 shadow-md dark:shadow-lg">
            <h3 className="text-xl font-semibold text-gray-800 dark:text-blue-300 mb-4">
              Similarity Result
            </h3>
            <div className="flex flex-col items-center mb-6">
              {/* Score */}
              <div
                className={`text-5xl font-bold mb-2 whitespace-nowrap ${getSimilarityColor(result)}`}
              >
                {result.toFixed(4)}
              </div>
              {/* Interpretation */}
              <div className="text-gray-600 dark:text-gray-400 text-center text-sm">
                {getSimilarityDescription(result)}
              </div>
              {/* Calculation Details */}
              <div className="text-gray-500 dark:text-gray-500 text-xs mt-2">
                Using {fingerprinter} fingerprints ({toolkit.toUpperCase()})
                {(fingerprinter === "ECFP" || fingerprinter === "MAPC") &&
                  ` - Radius: ${radius}, Bits: ${nBits}`}
                {fingerprinter !== "ECFP" && fingerprinter !== "MAPC" && ` - Bits: ${nBits}`}
              </div>
            </div>

            {/* Gauge Bar */}
            <div className="w-full h-3 bg-gray-200 dark:bg-gray-700 rounded-full overflow-hidden shadow-inner">
              <div
                className={`h-full rounded-full transition-all duration-500 ease-out ${getGaugeBgColor(result)}`}
                style={{ width: `${Math.max(1, result * 100)}%` }}
              ></div>
            </div>
            {/* Gauge Labels */}
            <div className="flex justify-between text-xs text-gray-500 dark:text-gray-400 mt-1 px-1">
              <span>0.0</span>
              <span>0.5</span>
              <span>1.0</span>
            </div>
          </div>

          {/* Molecule Comparison Cards */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            {/* Ensure MoleculeCard is theme-aware */}
            <MoleculeCard smiles={smiles1} title="Molecule 1" size="lg" />
            <MoleculeCard smiles={smiles2} title="Molecule 2" size="lg" />
          </div>
        </div>
      )}
    </div>
  );
};
export default TanimotoView;
