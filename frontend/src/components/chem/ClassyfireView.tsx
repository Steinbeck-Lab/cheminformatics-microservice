// Description: ClassyFire classification view component
import { useState, useEffect, useCallback } from "react";
// Ensure all used icons are imported
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { useAppContext } from "../../context/AppContext"; // Assuming AppContext provides apiConfig and addRecentMolecule
import { AlertCircle, Check, Clipboard, RefreshCw, Search, Loader2 } from "lucide-react";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { EmptyState } from "@/components/feedback/EmptyState";
import { getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";

const POLLING_INTERVAL = 5000; // 5 seconds

const ClassyfireView = () => {
  const [smiles, setSmiles] = useState("");
  const [jobId, setJobId] = useState(null);
  const [jobStatus, setJobStatus] = useState(null);
  const [classificationResults, setClassificationResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [polling, setPolling] = useState(false);
  const [error, setError] = useState(null);
  const [activeTab, setActiveTab] = useState("classify");
  const [manualJobId, setManualJobId] = useState("");
  const { addRecentMolecule, apiConfig } = useAppContext(); // Get necessary functions/data from context
  const [copied, setCopied] = useState(false);

  // Submit classification request
  const handleSubmitClassification = async (e) => {
    e.preventDefault();
    const trimmedSmiles = smiles.trim(); // Trim input
    if (!trimmedSmiles) {
      setError("Please enter a SMILES string");
      return;
    }

    setLoading(true);
    setError(null);
    setJobId(null);
    setJobStatus(null);
    setClassificationResults(null);

    try {
      // Use trimmedSmiles for the API call
      const response = await fetch(
        `${apiConfig.baseUrl}/chem/classyfire/classify?smiles=${encodeURIComponent(trimmedSmiles)}`
      );

      if (!response.ok) {
        // Try to get error details from response body
        let errorMsg = `Error ${response.status}: ${response.statusText}`;
        try {
          const errorData = await response.json();
          errorMsg = errorData.detail || errorMsg; // Use detail if available
        } catch (jsonError) {
          // Ignore if response is not JSON
        }
        throw new Error(errorMsg);
      }

      const data = await response.json();

      if (data && data.id) {
        setJobId(data.id);
        setJobStatus("submitted");
        setPolling(true); // Start polling after successful submission

        // Add molecule to recent list (ensure addRecentMolecule is provided by context)
        addRecentMolecule({
          smiles: trimmedSmiles, // Use trimmed SMILES
          name: "ClassyFire Molecule", // Or derive a name if possible
          timestamp: new Date().toISOString(),
        });
      } else {
        throw new Error("Invalid response from server: Missing job ID.");
      }
    } catch (err) {
      console.error("Classification submission error:", err); // Log the actual error
      setError(getErrorMessage("chem", err));
      setJobId(null); // Clear Job ID on error
      setJobStatus(null);
    } finally {
      setLoading(false);
    }
  };

  // Fetch job results using useCallback
  const fetchClassificationResults = useCallback(
    async (id) => {
      // No need to setPolling(true) here, it's managed by the caller/useEffect
      try {
        const response = await fetch(`${apiConfig.baseUrl}/chem/classyfire/${id}/result`);

        if (!response.ok) {
          // Handle potential errors like 404 Not Found if ID is invalid
          let errorMsg = `Error ${response.status}: ${response.statusText}`;
          try {
            const errorData = await response.json();
            errorMsg = errorData.detail || errorMsg;
          } catch (jsonError) {
            // Ignore if response is not JSON
          }
          // If status is 404, provide a more specific message
          if (response.status === 404) {
            errorMsg = `Job ID "${id}" not found. Please check the ID and try again.`;
            setJobId(null); // Clear invalid job ID
            setPolling(false); // Stop polling for invalid ID
          }
          throw new Error(errorMsg);
        }

        const data = await response.json();

        if (data) {
          if (data.classification_status === "Done") {
            setClassificationResults(data);
            setJobStatus("completed");
            setPolling(false); // Stop polling only when successfully completed
          } else {
            // Update status but keep polling active
            setJobStatus(data.classification_status || "processing");
            // Ensure polling remains true if status is not 'Done'
            setPolling(true);
          }
        } else {
          // Handle cases where response is OK but data is missing/unexpected
          throw new Error("Received empty or invalid results data.");
        }
      } catch (err) {
        console.error(`Error fetching results for Job ID ${id}:`, err); // Log the actual error
        setError(getErrorMessage("chem", err));
        setPolling(false); // Stop polling on error
        setJobStatus("error"); // Set a specific error status
      }
      // Removed finally block that might prematurely set polling to false
    },
    [apiConfig.baseUrl]
  ); // Dependencies for useCallback

  // Poll for results
  useEffect(() => {
    let intervalId = null;

    if (polling && jobId) {
      // Define the polling function inside useEffect or ensure it's stable
      const poll = () => {
        console.log(`Polling for Job ID: ${jobId}, Status: ${jobStatus}`);
        fetchClassificationResults(jobId);
      };

      // Fetch immediately when polling starts or jobId changes while polling is true
      poll();

      // Set interval for subsequent polls
      intervalId = setInterval(poll, POLLING_INTERVAL);
    }

    // Cleanup function to clear the interval when component unmounts,
    // or when polling stops (polling becomes false or jobId becomes null)
    return () => {
      if (intervalId) {
        console.log(`Clearing interval for Job ID: ${jobId}`);
        clearInterval(intervalId);
      }
    };
    // Rerun effect if polling status or jobId changes
  }, [polling, jobId, fetchClassificationResults, jobStatus]); // Added jobStatus dependency to log current status

  // Handle manual job ID submission
  const handleSubmitJobId = (e) => {
    e.preventDefault();
    const trimmedJobId = manualJobId.trim(); // Trim whitespace
    if (!trimmedJobId) {
      setError("Please enter a valid job ID");
      return;
    }
    setJobId(trimmedJobId); // Use trimmed ID
    setPolling(true); // Start polling for the manually entered ID
    setError(null);
    setClassificationResults(null); // Clear previous results
    setJobStatus("fetching"); // Set initial status for manual fetch
    // fetchClassificationResults(trimmedJobId); // Initial fetch is now handled by the useEffect
  };

  // Handle copy to clipboard
  const handleCopy = (text) => {
    if (!navigator.clipboard) {
      setError("Clipboard API not available in this browser.");
      return;
    }
    navigator.clipboard
      .writeText(text)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000); // Reset copied state after 2 seconds
      })
      .catch((err) => {
        console.error("Failed to copy Job ID:", err);
        setError("Failed to copy Job ID to clipboard."); // Inform user
      });
  };

  // Function to manually trigger a refresh when not polling
  const handleManualRefresh = () => {
    if (jobId && !polling) {
      setError(null); // Clear previous errors
      setJobStatus("fetching"); // Show fetching status
      setPolling(true); // Restart polling
      // fetchClassificationResults(jobId); // Fetch is handled by useEffect
    }
  };

  return (
    // Main container with spacing
    <div className="space-y-6 p-4 md:p-6">
      {/* Input section card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg">
        {/* Card Title */}
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          ClassyFire Chemical Classification
        </h2>

        {/* Tabs */}
        <div className="mb-6">
          <div className="flex border-b border-gray-200 dark:border-gray-700">
            {/* Classify Tab Button */}
            <Button
              variant="ghost"
              className={`h-auto px-4 py-2 font-medium text-sm rounded-none focus:outline-hidden transition-colors duration-150 ${
                activeTab === "classify"
                  ? "text-sky-700 dark:text-sky-400 border-b-2 border-sky-600 dark:border-sky-400"
                  : "bg-transparent text-slate-500 dark:text-slate-400 hover:text-slate-700 dark:hover:text-white hover:border-b-2 hover:border-slate-300 dark:hover:border-slate-500"
              }`}
              onClick={() => setActiveTab("classify")}
              aria-current={activeTab === "classify" ? "page" : undefined}
            >
              Classify Molecule
            </Button>
            {/* Check Job ID Tab Button */}
            <Button
              variant="ghost"
              className={`h-auto px-4 py-2 font-medium text-sm rounded-none focus:outline-hidden transition-colors duration-150 ${
                activeTab === "jobid"
                  ? "text-sky-700 dark:text-sky-400 border-b-2 border-sky-600 dark:border-sky-400"
                  : "bg-transparent text-slate-500 dark:text-slate-400 hover:text-slate-700 dark:hover:text-white hover:border-b-2 hover:border-slate-300 dark:hover:border-slate-500"
              }`}
              onClick={() => setActiveTab("jobid")}
              aria-current={activeTab === "jobid" ? "page" : undefined}
            >
              Check Job Status
            </Button>
          </div>
        </div>

        {/* Tab Content: Classify */}
        {activeTab === "classify" ? (
          <form onSubmit={handleSubmitClassification} className="space-y-4">
            {/* Use the themed SMILESInput component */}
            <SMILESInput
              value={smiles}
              onChange={setSmiles} // Pass the setter function directly
              required
              // Add any necessary props for theme handling if SMILESInput needs them
            />

            {/* Submit Button */}
            <Button
              type="submit"
              disabled={!smiles.trim() || loading} // Disable if input is empty/whitespace or loading
              className={`w-full px-4 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-hidden focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !smiles.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-xs"
              }`}
            >
              {loading ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  Submitting...
                </>
              ) : (
                <>
                  <Search className="mr-2 h-5 w-5" aria-hidden="true" />
                  Classify Molecule
                </>
              )}
            </Button>
          </form>
        ) : (
          /* Tab Content: Check Job ID */
          <form onSubmit={handleSubmitJobId} className="space-y-4">
            <div>
              <label
                htmlFor="jobIdInput"
                className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
              >
                ClassyFire Job ID
              </label>
              <Input
                id="jobIdInput"
                type="text"
                value={manualJobId}
                onChange={(e) => setManualJobId(e.target.value)}
                placeholder="Enter job ID..."
                className="w-full"
                required
                aria-required="true"
              />
            </div>

            {/* Submit Button */}
            <Button
              type="submit"
              // Disable if no manualJobId, or if currently loading a new classification,
              // or if currently polling *this specific* job ID.
              disabled={!manualJobId.trim() || loading || (polling && jobId === manualJobId.trim())}
              className={`w-full px-4 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-hidden focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !manualJobId.trim() || loading || (polling && jobId === manualJobId.trim())
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-xs"
              }`}
            >
              <Search className="mr-2 h-5 w-5" aria-hidden="true" />
              {/* Indicate if currently checking this specific ID */}
              {polling && jobId === manualJobId.trim()
                ? "Checking..."
                : "Get Classification Results"}
            </Button>
          </form>
        )}
      </div>

      {/* Loading Indicator (only for initial submission) */}
      {loading && !classificationResults && <ToolSkeleton variant="molecule" />}

      {/* Error Message */}
      {error && !loading && (
        <GlassErrorCard
          message={error}
          onRetry={() => {
            setError(null);
            document.getElementById("smiles-input")?.focus();
          }}
        />
      )}

      {/* Job Status / Polling Indicator */}
      {/* Show this block if we have a job ID, are not in the initial loading state, and don't have final results yet */}
      {jobId && !loading && !classificationResults && (
        <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 flex items-start shadow-sm">
          {/* Spinner/Icon */}
          <div className="mr-4 shrink-0 pt-1">
            {/* Show spinner if polling OR in the initial 'fetching' or 'submitted' state */}
            {polling || ["fetching", "submitted", "processing"].includes(jobStatus) ? (
              <div
                className="animate-spin h-5 w-5 border-2 border-blue-500 dark:border-blue-400 border-t-transparent rounded-full"
                aria-label="Processing..."
              ></div>
            ) : (
              // Show refresh icon if polling stopped (e.g., error, or completed but results somehow cleared)
              <RefreshCw
                className="h-5 w-5 text-blue-600 dark:text-blue-400" // Non-interactive, just an indicator
                aria-hidden="true"
              />
            )}
          </div>
          {/* Status Text and Actions */}
          <div className="flex-1">
            <h3 className="text-md font-medium text-blue-800 dark:text-blue-300 mb-1">
              {/* Display more informative status */}
              Classification Job{" "}
              <span className="font-mono bg-gray-100 dark:bg-gray-700 px-1 rounded-sm">
                {jobId}
              </span>{" "}
              - Status:{" "}
              <span className="font-semibold capitalize">{jobStatus || "Waiting..."}</span>
            </h3>
            <p className="text-gray-600 dark:text-gray-300 text-sm">
              {jobStatus === "error"
                ? "An error occurred. Check error message above."
                : polling
                  ? "Classification is processing. Results will automatically update."
                  : "Processing stopped. Click refresh to check status again."}
            </p>
            {/* Action buttons */}
            <div className="mt-3 flex items-center text-sm space-x-4">
              <Button
                variant="ghost"
                onClick={() => handleCopy(jobId.toString())}
                className="flex items-center text-blue-600 dark:text-blue-400 hover:text-blue-800 dark:hover:text-blue-300"
              >
                <Clipboard className="mr-1 h-4 w-4" aria-hidden="true" />
                {copied ? "Copied!" : "Copy Job ID"}
              </Button>
              {/* Show Refresh button only if not actively polling and status is not completed/error */}
              {!polling && jobStatus !== "completed" && jobStatus !== "error" && (
                <Button
                  variant="ghost"
                  onClick={handleManualRefresh}
                  className="flex items-center text-blue-600 dark:text-blue-400 hover:text-blue-800 dark:hover:text-blue-300"
                >
                  <RefreshCw className="mr-1 h-4 w-4" aria-hidden="true" />
                  Refresh Status
                </Button>
              )}
            </div>
          </div>
        </div>
      )}

      {/* Results Section */}
      {classificationResults && (
        <div className="space-y-6">
          {/* Results Card */}
          <div className="bg-white dark:bg-gray-800 rounded-lg overflow-hidden shadow-md dark:shadow-lg">
            {/* Results Header */}
            <div className="p-4 bg-gray-50 dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700 flex justify-between items-center">
              <h3 className="text-lg font-semibold text-gray-800 dark:text-white">
                ClassyFire Classification Results
              </h3>
              {/* Status Badge */}
              <div className="flex items-center text-green-600 dark:text-green-400 text-sm font-medium">
                <Check className="mr-1 h-4 w-4" aria-hidden="true" />
                {classificationResults.classification_status || "Completed"}
              </div>
            </div>

            {/* Results Body */}
            {classificationResults.entities && classificationResults.entities.length > 0 ? (
              <div className="p-6 space-y-8">
                {classificationResults.entities.map((entity, index) => (
                  // Added role="region" and aria-labelledby for better structure
                  <section
                    key={entity.smiles || index}
                    className="border-b border-gray-200 dark:border-gray-700 pb-8 last:border-b-0 last:pb-0"
                    aria-labelledby={`entity-heading-${index}`}
                  >
                    {/* Structure and Classification Grid */}
                    <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
                      {/* Molecule Structure Card */}
                      <div className="lg:col-span-1">
                        {/* Ensure MoleculeCard is accessible and theme-aware */}
                        <MoleculeCard
                          smiles={entity.smiles}
                          title={`Input Structure ${index + 1}`}
                          size="sm" // Adjust size as needed
                          // Pass theme prop if needed by MoleculeCard
                        />
                      </div>

                      {/* Classification Details */}
                      <div className="lg:col-span-2 space-y-4">
                        <h4
                          id={`entity-heading-${index}`}
                          className="text-lg font-medium text-gray-800 dark:text-blue-400 mb-3"
                        >
                          Chemical Ontology
                        </h4>
                        {/* Render classification levels */}
                        {["kingdom", "superclass", "class_", "direct_parent"].map((level) => {
                          const levelData = entity[level];
                          const levelName =
                            level === "class_"
                              ? "Class"
                              : level.replace("_", " ").replace("direct ", ""); // Clean up label
                          return (
                            levelData && (
                              // Ontology Item Card
                              <div
                                key={level}
                                className="bg-gray-50 dark:bg-gray-900 p-3 rounded-lg border border-gray-200 dark:border-gray-700 shadow-xs"
                              >
                                <h5 className="text-xs font-semibold uppercase tracking-wider text-gray-500 dark:text-gray-400 mb-1">
                                  {levelName}
                                </h5>
                                <p className="text-base font-medium text-gray-800 dark:text-blue-300">
                                  {levelData.name}
                                </p>
                                {levelData.description && (
                                  <p className="text-sm text-gray-600 dark:text-gray-400 mt-1">
                                    {levelData.description}
                                  </p>
                                )}
                              </div>
                            )
                          );
                        })}
                        {/* Handle case where no levels are present */}
                        {!entity.kingdom &&
                          !entity.superclass &&
                          !entity.class_ &&
                          !entity.direct_parent && (
                            <p className="text-gray-500 dark:text-gray-400 italic">
                              No classification levels found.
                            </p>
                          )}
                      </div>
                    </div>

                    {/* Additional Details Section */}
                    <div className="mt-4 pt-4 border-t border-gray-200 dark:border-gray-700">
                      <h4 className="text-md font-semibold text-gray-700 dark:text-gray-300 mb-3">
                        Additional Details
                      </h4>
                      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                        {/* Molecular Framework */}
                        {entity.molecular_framework && (
                          <div className="bg-gray-50 dark:bg-gray-900 p-3 rounded-lg border border-gray-200 dark:border-gray-700 shadow-xs">
                            <h5 className="text-xs font-semibold uppercase tracking-wider text-gray-500 dark:text-gray-400 mb-1">
                              Molecular Framework
                            </h5>
                            <p className="text-sm text-gray-800 dark:text-blue-300">
                              {entity.molecular_framework}
                            </p>
                          </div>
                        )}

                        {/* Substituents */}
                        {entity.substituents && entity.substituents.length > 0 && (
                          <div className="bg-gray-50 dark:bg-gray-900 p-3 rounded-lg border border-gray-200 dark:border-gray-700 shadow-xs">
                            <h5 className="text-xs font-semibold uppercase tracking-wider text-gray-500 dark:text-gray-400 mb-1">
                              Substituents
                            </h5>
                            <div className="flex flex-wrap gap-2 mt-1">
                              {entity.substituents.map((substituent, idx) => (
                                // Substituent Tag
                                <span
                                  key={idx}
                                  className="inline-block bg-blue-100 dark:bg-blue-900/50 text-blue-800 dark:text-blue-300 px-2 py-1 rounded-full text-xs font-medium"
                                >
                                  {substituent}
                                </span>
                              ))}
                            </div>
                          </div>
                        )}

                        {/* Classification Lineage (Ancestors) */}
                        {entity.ancestors && entity.ancestors.length > 0 && (
                          <div className="bg-gray-50 dark:bg-gray-900 p-3 rounded-lg border border-gray-200 dark:border-gray-700 md:col-span-2 shadow-xs">
                            <h5 className="text-xs font-semibold uppercase tracking-wider text-gray-500 dark:text-gray-400 mb-1">
                              Classification Lineage
                            </h5>
                            {/* Using an ordered list for semantic hierarchy */}
                            <ol className="flex flex-wrap gap-x-1 items-center text-sm text-gray-600 dark:text-gray-400 mt-1 list-none p-0">
                              {entity.ancestors.map((ancestor, idx) => (
                                <li key={idx} className="flex items-center">
                                  <span className="text-gray-700 dark:text-gray-300">
                                    {ancestor}
                                  </span>
                                  {/* Separator */}
                                  {idx < entity.ancestors.length - 1 && (
                                    <span
                                      className="text-gray-400 dark:text-gray-500 mx-1"
                                      aria-hidden="true"
                                    >
                                      &gt;
                                    </span>
                                  )}
                                </li>
                              ))}
                            </ol>
                          </div>
                        )}
                      </div>
                    </div>
                  </section> // End of region for entity
                ))}
              </div>
            ) : (
              // No results message
              <div className="p-6 text-center">
                <p className="text-gray-500 dark:text-gray-400">
                  No classification entities found in the results.
                </p>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Initial State / Info Box */}
      {/* Show only if not loading, no error, no job submitted/pending, and no results */}
      {/* Using the detailed about section with light/dark mode styling */}
      {!classificationResults && !loading && !error && !jobId && (
        <div
          className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 shadow-sm"
          role="complementary"
        >
          <h3 className="text-lg font-medium text-blue-800 dark:text-blue-300 mb-3">
            About ClassyFire
          </h3>
          <p className="text-gray-700 dark:text-gray-300 mb-4">
            ClassyFire is a web-based application for automated structural classification of
            chemical entities. It uses a comprehensive, computable chemical taxonomy developed by
            the Wishart Lab at the University of Alberta.
          </p>
          <p className="text-gray-700 dark:text-gray-300">
            ClassyFire's classification system organizes compounds into:
          </p>
          {/* List styling */}
          <ul className="list-disc list-inside mt-2 space-y-1 text-gray-600 dark:text-gray-300 marker:text-blue-500 dark:marker:text-blue-400">
            <li>Kingdoms (highest level)</li>
            <li>Superclasses</li>
            <li>Classes</li>
            <li>Subclasses (and additional levels)</li>
            <li>Direct Parents (lowest level)</li>
          </ul>
          <p className="mt-4 text-xs text-gray-500 dark:text-gray-400">
            Note: Classification may take several minutes to complete. You can save the Job ID to
            check results later.
          </p>
        </div>
      )}
    </div>
  );
};
export default ClassyfireView;
