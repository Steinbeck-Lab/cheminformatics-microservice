// Description: This component provides an interface to lookup chemical structures by name, formula, or identifiers using the PubChem database
import React, { useState } from "react";
import {
  HiOutlineSearch,
  HiOutlineExclamationCircle,
  HiOutlineInformationCircle,
} from "react-icons/hi";
import LoadingScreen from "../common/LoadingScreen";
import MoleculeCard from "../common/MoleculeCard";
import SMILESDisplay from "../common/SMILESDisplay";
import { useAppContext } from "../../context/AppContext"; // For adding to recent molecules

// Import the service function from chemService
import { lookupPubChem } from "../../services/chemService";

const PubChemLookupView = () => {
  const [identifier, setIdentifier] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null);
  const { addRecentMolecule } = useAppContext();

  // Examples for each identifier type
  const examples = [
    {
      name: "Chemical Name",
      value: "aspirin",
      description: "Common or IUPAC name",
    },
    { name: "CID", value: "2244", description: "PubChem Compound ID" },
    {
      name: "InChI",
      value: "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
      description: "International Chemical Identifier",
    },
    {
      name: "InChIKey",
      value: "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
      description: "Hashed InChI",
    },
    {
      name: "CAS",
      value: "50-78-2",
      description: "Chemical Abstracts Service registry number",
    },
    { name: "Formula", value: "C9H8O4", description: "Molecular formula" },
    {
      name: "SMILES",
      value: "CC(=O)OC1=CC=CC=C1C(=O)O",
      description: "Simplified molecular-input line-entry system",
    },
  ];

  const handleSubmit = async (e) => {
    e.preventDefault();

    // Trim the identifier and check if it's empty
    const trimmedIdentifier = identifier.trim();
    if (!trimmedIdentifier) {
      setError("Please enter a chemical identifier");
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);

    try {
      const data = await lookupPubChem(trimmedIdentifier);
      setResult(data);

      // Add to recent molecules if lookup was successful
      if (data.success && data.canonical_smiles) {
        addRecentMolecule({
          smiles: data.canonical_smiles,
          name: trimmedIdentifier,
          timestamp: new Date().toISOString(),
        });
      }
    } catch (err) {
      console.error("PubChem lookup error:", err);
      setError(`${err.message}`);
    } finally {
      setLoading(false);
    }
  };

  const handleUseExample = (exampleValue) => {
    setIdentifier(exampleValue);
  };

  return (
    <div className="space-y-6">
      {/* Input Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">
          Chemical Structure Finder
        </h2>

        {/* Form */}
        <form onSubmit={handleSubmit} className="space-y-4">
          <div>
            <label
              htmlFor="identifier-input"
              className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1"
            >
              Chemical Identifier
            </label>
            <div className="relative">
              <input
                id="identifier-input"
                type="text"
                value={identifier}
                onChange={(e) => setIdentifier(e.target.value)}
                placeholder="Enter name, CID, InChI, InChIKey, CAS, formula, or SMILES..."
                className="w-full px-4 py-2 pl-10 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500"
                required
              />
              <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                <HiOutlineSearch className="h-5 w-5 text-gray-400 dark:text-gray-500" />
              </div>
            </div>
          </div>

          {/* Examples */}
          <div>
            <p className="text-sm text-gray-500 dark:text-gray-400 mb-2">Examples:</p>
            <div className="flex flex-wrap gap-2">
              {examples.map((example, index) => (
                <button
                  key={index}
                  type="button"
                  onClick={() => handleUseExample(example.value)}
                  className="inline-flex items-center px-2.5 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-indigo-500 dark:focus:ring-offset-gray-800"
                  title={`${example.name}: ${example.description}`}
                >
                  {example.name}
                </button>
              ))}
            </div>
          </div>

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              disabled={!identifier.trim() || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${
                !identifier.trim() || loading
                  ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                  : "bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm"
              }`}
            >
              <HiOutlineSearch className="mr-2 h-5 w-5" />
              {loading ? "Searching..." : "Find Structure"}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Searching PubChem..." />}

      {/* Error Display */}
      {error && !loading && (
        <div
          className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow-sm"
          role="alert"
        >
          <HiOutlineExclamationCircle
            className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
            aria-hidden="true"
          />
          <span>{error}</span>
        </div>
      )}

      {/* Results Display */}
      {result && !loading && (
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
          <h3 className="text-lg font-semibold text-gray-800 dark:text-white mb-4">
            Lookup Results
          </h3>

          {/* Success or Failure Message */}
          <div
            className={`p-4 mb-6 rounded-md flex items-start ${
              result.success
                ? "bg-green-50 dark:bg-green-900 dark:bg-opacity-20 text-green-700 dark:text-green-200 border border-green-200 dark:border-green-800"
                : "bg-yellow-50 dark:bg-yellow-900 dark:bg-opacity-20 text-yellow-700 dark:text-yellow-200 border border-yellow-200 dark:border-yellow-800"
            }`}
          >
            <HiOutlineInformationCircle
              className={`h-5 w-5 mr-3 flex-shrink-0 mt-0.5 ${
                result.success
                  ? "text-green-500 dark:text-green-400"
                  : "text-yellow-500 dark:text-yellow-400"
              }`}
              aria-hidden="true"
            />
            <div className="flex-1">
              <p className="text-sm font-medium">
                {result.success
                  ? "Successfully retrieved structure from PubChem"
                  : "No structure found in PubChem for this identifier"}
              </p>
              <p className="text-sm mt-1">
                Input identifier: <span className="font-mono">{result.input}</span> (detected as{" "}
                {result.input_type})
              </p>
              {result.success && result.cids && (
                <div className="mt-3 space-y-2">
                  <p className="text-sm font-medium">
                    PubChem CID{result.cids.length > 1 ? "s" : ""}:{" "}
                    {result.cids.slice(0, 5).join(", ")}
                    {result.cids.length > 5 && ` and ${result.cids.length - 5} more`}
                  </p>
                  <div className="flex flex-wrap gap-2">
                    {result.pubchem_links &&
                      result.pubchem_links.slice(0, 3).map((link, index) => (
                        <a
                          key={index}
                          href={link}
                          target="_blank"
                          rel="noopener noreferrer"
                          className="inline-flex items-center px-3 py-1 text-xs font-medium rounded-full bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300 hover:bg-blue-200 dark:hover:bg-blue-800 transition-colors duration-200"
                        >
                          View CID {result.cids[index]} in PubChem
                          <svg
                            className="w-3 h-3 ml-1"
                            fill="none"
                            stroke="currentColor"
                            viewBox="0 0 24 24"
                          >
                            <path
                              strokeLinecap="round"
                              strokeLinejoin="round"
                              strokeWidth={2}
                              d="M10 6H6a2 2 0 00-2 2v10a2 2 0 002 2h10a2 2 0 002-2v-4M14 4h6m0 0v6m0-6L10 14"
                            />
                          </svg>
                        </a>
                      ))}
                    {result.pubchem_links && result.pubchem_links.length > 3 && (
                      <span className="inline-flex items-center px-3 py-1 text-xs font-medium rounded-full bg-gray-100 dark:bg-gray-700 text-gray-600 dark:text-gray-400">
                        +{result.pubchem_links.length - 3} more compounds
                      </span>
                    )}
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* SMILES Result */}
          {result.success && result.canonical_smiles && (
            <>
              <div className="mb-6">
                <div className="flex justify-between items-center mb-2">
                  <h4 className="text-md font-medium text-gray-700 dark:text-gray-300">
                    Canonical SMILES
                  </h4>
                </div>
                <SMILESDisplay smiles={result.canonical_smiles} label="" showDownload={true} />
              </div>

              {/* Molecule Visualization */}
              <div className="mt-4">
                <h4 className="text-md font-medium text-gray-700 dark:text-gray-300 mb-3">
                  Structure Visualization
                </h4>
                <div className="flex justify-center">
                  <MoleculeCard smiles={result.canonical_smiles} title={result.input} size="lg" />
                </div>
              </div>
            </>
          )}
        </div>
      )}

      {/* Information Box */}
      <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-4 text-sm shadow">
        <h4 className="font-semibold text-blue-800 dark:text-blue-300 mb-2 flex items-center">
          <HiOutlineInformationCircle className="h-5 w-5 mr-2 text-blue-500 dark:text-blue-400" />
          About Chemical Structure Finder
        </h4>
        <p className="text-gray-700 dark:text-gray-300 mb-2">
          This tool helps you find chemical structures using a variety of inputs. Just enter what
          you know about the compound:
        </p>
        <ul className="list-disc list-inside space-y-1 text-gray-600 dark:text-gray-400 ml-2">
          <li>Chemical names (common names, trade names, or IUPAC)</li>
          <li>PubChem CIDs (Compound IDs)</li>
          <li>InChI or InChIKey strings</li>
          <li>CAS registry numbers</li>
          <li>Molecular formulas</li>
          <li>SMILES notations</li>
        </ul>
        <p className="text-gray-700 dark:text-gray-300 mt-2">
          The service automatically detects what you've entered and searches the PubChem database,
          which contains over 100 million chemical compounds.
        </p>
      </div>
    </div>
  );
};

export default PubChemLookupView;
