// Description: This component combines PubChem lookup with 2D and 3D visualization
import React, { useState } from 'react';
import {
  HiOutlineSearch,
  HiOutlineExclamationCircle,
  HiOutlineInformationCircle,
  HiOutlineDocumentDuplicate,
  HiOutlineCheck,
  HiOutlineBeaker,
  HiOutlineChartSquareBar,
  HiOutlineCube
} from 'react-icons/hi';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faAtom, faCircle } from '@fortawesome/free-solid-svg-icons';
import LoadingScreen from '../common/LoadingScreen';
import SMILESDisplay from '../common/SMILESDisplay';
import { useAppContext } from '../../context/AppContext'; // For adding to recent molecules
import { lookupPubChem } from '../../services/chemService';

// Import visualization components
import MoleculeDepiction2D from './MoleculeDepiction2D';
import MoleculeDepiction3D from './MoleculeDepiction3D';

// Animated Atom Component
const AnimatedAtom = () => {
  return (
    <div className="w-full h-full flex items-center justify-center relative">
      {/* Main nucleus */}
      <div className="text-indigo-500 dark:text-indigo-400 z-10">
        <FontAwesomeIcon icon={faAtom} className="text-6xl animate-pulse" />
      </div>
      
      {/* Orbiting electrons */}
      <div className="absolute inset-0 w-full h-full">
        {/* First orbit */}
        <div className="absolute inset-0 w-full h-full animate-spin-slow">
          <div className="absolute top-0 left-1/2 transform -translate-x-1/2 -translate-y-1/2">
            <FontAwesomeIcon icon={faCircle} className="text-blue-500 text-xs" />
          </div>
          <div className="absolute bottom-0 left-1/2 transform -translate-x-1/2 translate-y-1/2">
            <FontAwesomeIcon icon={faCircle} className="text-blue-500 text-xs" />
          </div>
        </div>
        
        {/* Second orbit - rotated */}
        <div className="absolute inset-0 w-full h-full animate-spin-reverse">
          <div className="absolute top-1/2 left-0 transform -translate-y-1/2 -translate-x-1/2">
            <FontAwesomeIcon icon={faCircle} className="text-green-500 text-xs" />
          </div>
          <div className="absolute top-1/2 right-0 transform -translate-y-1/2 translate-x-1/2">
            <FontAwesomeIcon icon={faCircle} className="text-green-500 text-xs" />
          </div>
        </div>
        
        {/* Third orbit - tilted */}
        <div className="absolute inset-0 w-full h-full rotate-45 animate-spin-medium">
          <div className="absolute top-0 right-0 transform -translate-y-1/2 translate-x-1/2">
            <FontAwesomeIcon icon={faCircle} className="text-purple-500 text-xs" />
          </div>
          <div className="absolute bottom-0 left-0 transform translate-y-1/2 -translate-x-1/2">
            <FontAwesomeIcon icon={faCircle} className="text-purple-500 text-xs" />
          </div>
        </div>
      </div>
      
      {/* Pulse ring animation */}
      <div className="absolute inset-0 rounded-full border-2 border-indigo-200 dark:border-indigo-800 animate-ping opacity-75"></div>
    </div>
  );
};

const StructureVisualizerView = () => {
  // Search state
  const [identifier, setIdentifier] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null);
  const [copySuccess, setCopySuccess] = useState(false);
  const { addRecentMolecule } = useAppContext();
  
  // Visualization settings
  const [depictToolkit, setDepictToolkit] = useState('rdkit');
  const [vis3DToolkit, setVis3DToolkit] = useState('openbabel');

  // Examples for each identifier type
  const examples = [
    { name: 'Caffeine', value: 'caffeine', description: 'Common stimulant' },
    { name: 'Aspirin', value: 'aspirin', description: 'Common pain reliever' },
    { name: 'Ibuprofen', value: 'ibuprofen', description: 'NSAID pain reliever' },
    { name: 'Paracetamol', value: 'paracetamol', description: 'Also known as acetaminophen' },
    { name: 'CAS: 50-78-2', value: '50-78-2', description: 'CAS for aspirin' },
    { name: 'Formula: C8H10N4O2', value: 'C8H10N4O2', description: 'Caffeine formula' },
  ];

  const handleSubmit = async (e) => {
    e.preventDefault();
    
    // Trim the identifier and check if it's empty
    const trimmedIdentifier = identifier.trim();
    if (!trimmedIdentifier) {
      setError('Please enter a chemical identifier');
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
          name: data.name || trimmedIdentifier,
          timestamp: new Date().toISOString()
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

  const copyToClipboard = (text) => {
    navigator.clipboard.writeText(text).then(
      () => {
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
      },
      (err) => {
        console.error('Failed to copy text:', err);
      }
    );
  };

  return (
    <div className="flex flex-col gap-6 p-4 md:p-6 bg-gradient-to-br from-sky-50 to-indigo-100 dark:from-gray-900 dark:to-slate-900 min-h-screen">
      {/* Header with animated background */}
      <div className="bg-gradient-to-r from-blue-600 to-purple-700 dark:from-blue-700 dark:to-purple-900 rounded-xl shadow-xl overflow-hidden relative">
        {/* Animated background elements */}
        <div className="absolute inset-0 overflow-hidden opacity-10">
          <div className="absolute left-0 top-0 w-40 h-40 rounded-full bg-white transform -translate-x-1/2 -translate-y-1/2"></div>
          <div className="absolute right-20 bottom-0 w-64 h-64 rounded-full bg-white transform translate-x-1/2 translate-y-1/2"></div>
          <div className="absolute left-1/3 top-1/2 w-32 h-32 rounded-full bg-white transform -translate-y-1/2"></div>
        </div>

        <div className="relative p-8 md:p-10 flex flex-col md:flex-row items-center justify-between gap-6">
          {/* Title and description */}
          <div className="text-white space-y-2 text-center md:text-left">
            <h1 className="text-3xl md:text-4xl font-bold">Chemical Structure Visualizer</h1>
            <p className="text-blue-100 dark:text-blue-200 opacity-90 max-w-xl">
              Search, visualize, and explore chemical structures in 2D and 3D
            </p>
          </div>
          
          {/* Decorative icon */}
          <div className="hidden md:flex items-center justify-center bg-white/10 backdrop-blur-sm p-6 rounded-full w-24 h-24 border border-white/20 shadow-lg">
            <HiOutlineBeaker className="w-12 h-12 text-white" />
          </div>
        </div>
      </div>

      {/* Main content area */}
      <div className="grid grid-cols-1 xl:grid-cols-12 gap-6">
        {/* Search panel - left side on larger screens */}
        <div className="xl:col-span-4 xl:row-span-2 flex flex-col gap-4">
          {/* Search Card */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-xl border border-gray-100 dark:border-gray-700 backdrop-blur-sm flex-grow">
            <form onSubmit={handleSubmit} className="space-y-5">
              <div>
                <div className="flex items-center justify-between mb-2">
                  <label htmlFor="identifier-input" className="block text-sm font-medium text-gray-700 dark:text-gray-300">
                    Chemical Identifier
                  </label>
                  <span className="text-xs text-gray-500 dark:text-gray-400">Enter name, CAS, formula, SMILES...</span>
                </div>
                <div className="relative">
                  <input
                    id="identifier-input"
                    type="text"
                    value={identifier}
                    onChange={(e) => setIdentifier(e.target.value)}
                    placeholder="Search for a chemical compound..."
                    className="w-full px-4 py-3 pl-10 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-indigo-500 focus:border-indigo-500 transition-all duration-200"
                    required
                  />
                  <div className="absolute inset-y-0 left-0 pl-3 flex items-center pointer-events-none">
                    <HiOutlineSearch className="h-5 w-5 text-gray-400 dark:text-gray-500" />
                  </div>
                </div>
              </div>

              {/* Examples */}
              <div>
                <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                  Quick Examples:
                </p>
                <div className="flex flex-wrap gap-2">
                  {examples.map((example, index) => (
                    <button
                      key={index}
                      type="button"
                      onClick={() => handleUseExample(example.value)}
                      className="inline-flex items-center px-3 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded-full text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-indigo-500 dark:focus:ring-offset-gray-800 transition-all duration-200"
                      title={`${example.name}: ${example.description}`}
                    >
                      {example.name}
                    </button>
                  ))}
                </div>
              </div>

              {/* Action Buttons */}
              <div className="pt-4">
                {/* Submit Button */}
                <button
                  type="submit"
                  disabled={!identifier.trim() || loading}
                  className={`w-full relative overflow-hidden px-6 py-3 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 shadow-lg ${
                    !identifier.trim() || loading
                      ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                      : 'bg-gradient-to-r from-blue-600 to-indigo-700 hover:from-blue-700 hover:to-indigo-800 transform hover:-translate-y-1'
                  }`}
                >
                  <HiOutlineSearch className="mr-2 h-5 w-5" />
                  {loading ? 'Searching...' : 'Search & Visualize'}
                  {!loading && !identifier.trim() ? null : (
                    <span className="absolute inset-0 overflow-hidden rounded-lg">
                      <span className="absolute inset-0 rounded-lg bg-gradient-to-r from-blue-400 to-indigo-500 opacity-0 transition-opacity hover:opacity-20"></span>
                    </span>
                  )}
                </button>
                
                {/* Toolkit Selectors */}
                <div className="grid grid-cols-2 gap-4 mt-4">
                  <div className="bg-gray-50 dark:bg-gray-700/50 p-3 rounded-lg border border-gray-200 dark:border-gray-600">
                    <div className="flex items-center mb-2">
                      <HiOutlineChartSquareBar className="h-4 w-4 text-indigo-500 dark:text-indigo-400 mr-2" />
                      <label htmlFor="depict-toolkit" className="text-sm font-medium text-gray-700 dark:text-gray-300">
                        2D Toolkit
                      </label>
                    </div>
                    <select
                      id="depict-toolkit"
                      value={depictToolkit}
                      onChange={(e) => setDepictToolkit(e.target.value)}
                      className="w-full bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white text-sm px-3 py-2 focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                    >
                      <option value="rdkit">RDKit</option>
                      <option value="cdk">CDK</option>
                    </select>
                  </div>
                  
                  <div className="bg-gray-50 dark:bg-gray-700/50 p-3 rounded-lg border border-gray-200 dark:border-gray-600">
                    <div className="flex items-center mb-2">
                      <HiOutlineCube className="h-4 w-4 text-indigo-500 dark:text-indigo-400 mr-2" />
                      <label htmlFor="vis3d-toolkit" className="text-sm font-medium text-gray-700 dark:text-gray-300">
                        3D Toolkit
                      </label>
                    </div>
                    <select
                      id="vis3d-toolkit"
                      value={vis3DToolkit}
                      onChange={(e) => setVis3DToolkit(e.target.value)}
                      className="w-full bg-white dark:bg-gray-800 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white text-sm px-3 py-2 focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
                    >
                      <option value="openbabel">Open Babel</option>
                      <option value="rdkit">RDKit</option>
                    </select>
                  </div>
                </div>
              </div>
            </form>
          </div>

          {/* Information Box */}
          <div className="bg-gradient-to-br from-blue-50 to-indigo-50 dark:from-blue-900/30 dark:to-indigo-900/30 border border-blue-200 dark:border-blue-800/50 rounded-xl p-5 text-sm shadow-lg">
            <h4 className="font-bold text-blue-800 dark:text-blue-300 mb-3 flex items-center">
              <HiOutlineInformationCircle className="h-5 w-5 mr-2 text-blue-500 dark:text-blue-400" />
              About This Tool
            </h4>
            <div className="space-y-3 text-gray-700 dark:text-gray-300">
              <p>
                This visualizer helps you explore chemical structures using data from PubChem, one of the world's largest collections of freely accessible chemical information.
              </p>
              <div>
                <h5 className="font-medium mb-1 text-gray-800 dark:text-gray-200">Features:</h5>
                <ul className="list-disc list-inside space-y-1 pl-1 text-gray-600 dark:text-gray-400">
                  <li>Search by name, CAS number, formula, InChI, or SMILES</li>
                  <li>View precise 2D structural diagrams</li>
                  <li>Explore interactive 3D molecular models</li>
                  <li>Switch between rendering engines for different views</li>
                </ul>
              </div>
            </div>
          </div>
        </div>

        {/* Main content area - right side / bottom */}
        <div className="xl:col-span-8 space-y-6">
          {/* Loading State */}
          {loading && <LoadingScreen text="Searching and preparing visualizations..." />}

          {/* Error Display */}
          {error && !loading && (
            <div className="p-5 rounded-xl bg-red-50 dark:bg-red-900/30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700/50 flex items-start shadow-lg" role="alert">
              <HiOutlineExclamationCircle className="h-6 w-6 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
              <div>
                <h3 className="font-medium text-lg text-red-800 dark:text-red-300 mb-1">Search Error</h3>
                <p>{error}</p>
                <p className="mt-2 text-sm">Try using a different identifier or check your spelling.</p>
              </div>
            </div>
          )}

          {/* Results Display */}
          {result && !loading && (
            <div className="space-y-6 animate-fadeIn">
              {/* Status Card */}
              <div className={`p-4 rounded-xl flex items-start shadow-lg ${
                result.success 
                  ? 'bg-gradient-to-br from-green-50 to-emerald-50 dark:from-green-900/20 dark:to-emerald-900/20 text-green-700 dark:text-green-200 border border-green-200 dark:border-green-800/50'
                  : 'bg-gradient-to-br from-yellow-50 to-amber-50 dark:from-yellow-900/20 dark:to-amber-900/20 text-yellow-700 dark:text-yellow-200 border border-yellow-200 dark:border-yellow-800/50'
              }`}>
                <HiOutlineInformationCircle 
                  className={`h-6 w-6 mr-3 flex-shrink-0 mt-0.5 ${
                    result.success 
                      ? 'text-green-500 dark:text-green-400' 
                      : 'text-yellow-500 dark:text-yellow-400'
                  }`} 
                  aria-hidden="true" 
                />
                <div>
                  <h3 className={`font-medium text-lg ${
                    result.success 
                      ? 'text-green-800 dark:text-green-300' 
                      : 'text-yellow-800 dark:text-yellow-300'
                  } mb-1`}>
                    {result.success 
                      ? `Found: ${result.name || result.input}` 
                      : 'No structure found in PubChem'}
                  </h3>
                  <p className="text-sm">
                    Input identifier: <span className="font-mono bg-white/50 dark:bg-black/20 px-1.5 py-0.5 rounded">{result.input}</span>
                    <span className="ml-2 px-2 py-0.5 rounded-full text-xs bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300">
                      {result.input_type}
                    </span>
                  </p>
                </div>
              </div>

              {/* SMILES Result & Visualizations */}
              {result.success && result.canonical_smiles && (
                <div className="bg-white dark:bg-gray-800 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700 overflow-hidden">
                  {/* SMILES Header */}
                  <div className="bg-gradient-to-r from-gray-50 to-gray-100 dark:from-gray-800 dark:to-gray-750 px-6 py-4 border-b border-gray-200 dark:border-gray-700">
                    <div className="flex justify-between items-center">
                      <h3 className="text-lg font-bold text-gray-800 dark:text-white flex items-center">
                        <span className="bg-indigo-100 dark:bg-indigo-900/50 p-1.5 rounded-lg mr-2">
                          <svg className="w-5 h-5 text-indigo-700 dark:text-indigo-400" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
                            <path d="M14.5 9C14.5 9 13.7609 8 11 8C8.23908 8 6.5 9.5 6.5 13.5C6.5 17.5 7.92954 19 11 19C14.0705 19 14.5 17 14.5 17" stroke="currentColor" strokeWidth="2" strokeLinecap="round"/>
                            <circle cx="17" cy="7" r="3" stroke="currentColor" strokeWidth="2"/>
                            <circle cx="17" cy="17" r="3" stroke="currentColor" strokeWidth="2"/>
                            <path d="M14 7H6.5C4.567 7 3 8.567 3 10.5V19" stroke="currentColor" strokeWidth="2" strokeLinecap="round"/>
                          </svg>
                        </span>
                        Canonical SMILES
                      </h3>
                      <button
                        onClick={() => copyToClipboard(result.canonical_smiles)}
                        className="flex items-center text-sm bg-white dark:bg-gray-700 px-3 py-1.5 rounded-lg border border-gray-200 dark:border-gray-600 text-gray-700 dark:text-gray-300 hover:bg-gray-50 dark:hover:bg-gray-600 transition-colors shadow-sm"
                      >
                        {copySuccess ? (
                          <>
                            <HiOutlineCheck className="h-4 w-4 mr-1.5 text-green-500" />
                            <span>Copied</span>
                          </>
                        ) : (
                          <>
                            <HiOutlineDocumentDuplicate className="h-4 w-4 mr-1.5" />
                            <span>Copy SMILES</span>
                          </>
                        )}
                      </button>
                    </div>
                  </div>
                  
                  {/* Content */}
                  <div className="p-6 space-y-6">
                    {/* SMILES Display */}
                    <div className="bg-gray-50 dark:bg-gray-900/50 p-4 rounded-lg border border-gray-200 dark:border-gray-700 overflow-x-auto">
                      <SMILESDisplay 
                        smiles={result.canonical_smiles} 
                        label="" 
                        showDownload={true}
                      />
                    </div>

                    {/* Visualizations Grid */}
                    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                      <MoleculeDepiction2D 
                        smiles={result.canonical_smiles} 
                        title={`2D Structure: ${result.name || result.input}`}
                        toolkit={depictToolkit} 
                      />
                      <MoleculeDepiction3D 
                        smiles={result.canonical_smiles} 
                        title={`3D Structure: ${result.name || result.input}`}
                        toolkit={vis3DToolkit} 
                      />
                    </div>
                  </div>
                </div>
              )}
            </div>
          )}

          {/* Initial State or No Results */}
          {!result && !loading && !error && (
            <div className="bg-white dark:bg-gray-800 p-8 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700 flex flex-col items-center justify-center text-center h-full min-h-[400px]">
              <div className="w-24 h-24 mb-6 relative">
                <AnimatedAtom />
              </div>
              
              <h3 className="text-2xl font-bold text-gray-800 dark:text-white mb-2">Begin Your Search</h3>
              <p className="text-gray-500 dark:text-gray-400 max-w-md mb-6">
                Enter a chemical name, formula, CAS number, or SMILES to search PubChem and visualize molecular structures
              </p>
              
              <div className="flex flex-wrap gap-2 justify-center">
                <button
                  onClick={() => handleUseExample('caffeine')}
                  className="px-4 py-2 bg-indigo-50 dark:bg-indigo-900/30 text-indigo-700 dark:text-indigo-300 rounded-lg hover:bg-indigo-100 dark:hover:bg-indigo-900/50 transition-colors"
                >
                  Try Caffeine
                </button>
                <button
                  onClick={() => handleUseExample('aspirin')}
                  className="px-4 py-2 bg-indigo-50 dark:bg-indigo-900/30 text-indigo-700 dark:text-indigo-300 rounded-lg hover:bg-indigo-100 dark:hover:bg-indigo-900/50 transition-colors"
                >
                  Try Aspirin
                </button>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

// Add the custom atom animation CSS
const styles = `
@keyframes fadeIn {
  from { opacity: 0; transform: translateY(10px); }
  to { opacity: 1; transform: translateY(0); }
}
.animate-fadeIn {
  animation: fadeIn 0.5s ease-out forwards;
}

@keyframes spin-slow {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(360deg); }
}
@keyframes spin-reverse {
  0% { transform: rotate(0deg); }
  100% { transform: rotate(-360deg); }
}
@keyframes spin-medium {
  0% { transform: rotate(45deg); }
  100% { transform: rotate(405deg); }
}
.animate-spin-slow {
  animation: spin-slow 8s linear infinite;
}
.animate-spin-reverse {
  animation: spin-reverse 10s linear infinite;
}
.animate-spin-medium {
  animation: spin-medium 12s linear infinite;
}
`;

// Add the styles to the document
if (typeof document !== 'undefined') {
  const styleSheet = document.createElement('style');
  styleSheet.type = 'text/css';
  styleSheet.innerText = styles;
  document.head.appendChild(styleSheet);
}

export default StructureVisualizerView;