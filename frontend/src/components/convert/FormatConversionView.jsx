// Description: This component handles the format conversion between different chemical notations.
import React, { useState } from 'react';
// Ensure all used icons are imported
import {
  HiOutlineSwitchHorizontal,
  HiOutlineClipboard,
  HiOutlineCheck,
  HiOutlineExclamationCircle,
  HiOutlineArrowRight
} from 'react-icons/hi';
// Assuming these components are correctly implemented and styled for dark/light mode
import SMILESInput from '../common/SMILESInput';
import LoadingScreen from '../common/LoadingScreen';
// Assuming this service is configured correctly
import convertService from '../../services/convertService';

// Input format options configuration
const INPUT_FORMAT_OPTIONS = [
  { id: 'smiles', label: 'SMILES' },
  { id: 'iupac', label: 'IUPAC Name' },
  { id: 'selfies', label: 'SELFIES' }
];

// Output format options configuration
const OUTPUT_FORMAT_OPTIONS = [
  { id: 'smiles', label: 'SMILES', method: null },
  { id: 'canonicalsmiles', label: 'Canonical SMILES', method: 'generateCanonicalSMILES' },
  { id: 'inchi', label: 'InChI', method: 'generateInChI' },
  { id: 'inchikey', label: 'InChI Key', method: 'generateInChIKey' },
  { id: 'cxsmiles', label: 'CXSMILES', method: 'generateCXSMILES' },
  { id: 'selfies', label: 'SELFIES', method: 'generateSELFIES' }
];

// Toolkit options configuration
const TOOLKIT_OPTIONS = [
  { id: 'cdk', label: 'CDK (Chemistry Development Kit)' },
  { id: 'rdkit', label: 'RDKit' },
  { id: 'openbabel', label: 'OpenBabel' }
];

// Converter options for IUPAC
const IUPAC_CONVERTER_OPTIONS = [
  { id: 'opsin', label: 'OPSIN' }
];

const FormatConversionView = () => {
  const [input, setInput] = useState('');
  const [inputFormat, setInputFormat] = useState('smiles');
  const [outputFormat, setOutputFormat] = useState('canonicalsmiles');
  const [toolkit, setToolkit] = useState('cdk');
  const [iupacConverter, setIupacConverter] = useState('opsin');
  const [result, setResult] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [copied, setCopied] = useState(false);

  // When input format changes, automatically update output format if needed
  const handleInputFormatChange = (format) => {
    setInputFormat(format);

    // If switching to IUPAC or SELFIES, automatically set output to SMILES
    if (format === 'iupac' || format === 'selfies') {
      setOutputFormat('smiles');
    }
  };

  // Handle form submission for conversion
  const handleSubmit = async (e) => {
    e.preventDefault();
    const trimmedInput = input.trim();
    if (!trimmedInput) {
      setError('Please enter input data.');
      setResult('');
      return;
    }

    setLoading(true);
    setError(null);
    setResult('');

    try {
      let convertedResult;

      // Handle IUPAC to SMILES or SELFIES to SMILES conversion
      if (inputFormat !== 'smiles') {
        // First convert IUPAC or SELFIES to SMILES
        const smiles = await convertService.generateSMILES(
          trimmedInput,
          inputFormat,
          inputFormat === 'iupac' ? iupacConverter : undefined
        );

        // If the output is SMILES, we're done
        if (outputFormat === 'smiles') {
          convertedResult = smiles;
        } else {
          // Otherwise, convert SMILES to the target format
          const formatOption = OUTPUT_FORMAT_OPTIONS.find(option => option.id === outputFormat);
          if (!formatOption || !formatOption.method) {
            throw new Error(`Unsupported output format: ${outputFormat}`);
          }

          const method = convertService[formatOption.method];
          if (typeof method !== 'function') {
            throw new Error(`Conversion function not available for format: ${outputFormat}`);
          }

          // Convert the SMILES to the target format
          convertedResult = await method(smiles, toolkit);
        }
      } else {
        // Direct SMILES conversion to target format
        if (outputFormat === 'smiles') {
          // Just return the input if output is also SMILES
          convertedResult = trimmedInput;
        } else {
          const formatOption = OUTPUT_FORMAT_OPTIONS.find(option => option.id === outputFormat);
          if (!formatOption || !formatOption.method) {
            throw new Error(`Unsupported output format: ${outputFormat}`);
          }

          const method = convertService[formatOption.method];
          if (typeof method !== 'function') {
            throw new Error(`Conversion function not available for format: ${outputFormat}`);
          }

          convertedResult = await method(trimmedInput, toolkit);
        }
      }

      // Handle cases where conversion might return null/undefined/empty
      if (!convertedResult) {
        setError(`Conversion resulted in empty output.`);
        setResult('');
      } else {
        let finalResult = String(convertedResult);

        // Remove surrounding double quotes if present
        if (finalResult.length >= 2 && finalResult.startsWith('"') && finalResult.endsWith('"')) {
          finalResult = finalResult.substring(1, finalResult.length - 1);
        }

        setResult(finalResult);
      }
    } catch (err) {
      console.error("Conversion failed:", err);
      setError(`Conversion failed: ${err.message || 'An unknown error occurred.'}`);
      setResult('');
    } finally {
      setLoading(false);
    }
  };

  // Handle copying the result to clipboard
  const handleCopyResult = () => {
    if (!result || !navigator.clipboard) return;

    navigator.clipboard.writeText(result)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      })
      .catch(err => {
        console.error('Failed to copy result:', err);
        setError('Failed to copy result to clipboard.');
      });
  };

  // Determine if toolkit selection should be shown based on input/output format
  const showToolkitSelection = inputFormat === 'smiles' &&
    outputFormat !== 'selfies'; // Some conversions don't need toolkit selection

  // Determine if IUPAC converter selection should be shown
  const showIupacConverterSelection = inputFormat === 'iupac';

  return (
    <div className="space-y-6 p-4 md:p-6">
      {/* Input and Options Card */}
      <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
        <h2 className="text-xl font-semibold text-gray-800 dark:text-blue-400 mb-4">Format Conversion</h2>

        <form onSubmit={handleSubmit} className="space-y-4">
          {/* Input Format Selection */}
          <div>
            <label htmlFor="input-format-select" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
              Input Format
            </label>
            <select
              id="input-format-select"
              value={inputFormat}
              onChange={(e) => handleInputFormatChange(e.target.value)}
              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
            >
              {INPUT_FORMAT_OPTIONS.map((option) => (
                <option key={option.id} value={option.id}>
                  {option.label}
                </option>
              ))}
            </select>
          </div>

          {/* IUPAC Converter Selection (conditionally shown) */}
          {showIupacConverterSelection && (
            <div>
              <label htmlFor="iupac-converter-select" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                IUPAC Converter
              </label>
              <select
                id="iupac-converter-select"
                value={iupacConverter}
                onChange={(e) => setIupacConverter(e.target.value)}
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                {IUPAC_CONVERTER_OPTIONS.map((option) => (
                  <option key={option.id} value={option.id}>
                    {option.label}
                  </option>
                ))}
              </select>
              <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
                OPSIN is used to convert IUPAC names to SMILES
              </p>
            </div>
          )}

          {/* Input Field using SMILESInput component */}
          <div>
            <SMILESInput
              value={input}
              onChange={setInput}
              label={inputFormat === 'smiles' ? 'SMILES Input' :
                inputFormat === 'iupac' ? 'IUPAC Name' : 'SELFIES Input'}
              placeholder={inputFormat === 'smiles' ? 'Enter SMILES notation...' :
                inputFormat === 'iupac' ? 'Enter IUPAC chemical name...' :
                  'Enter SELFIES notation...'}
              required
            />

            {inputFormat === 'iupac' && (
              <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
                Example: 1,3,7-trimethylpurine-2,6-dione (caffeine)
              </p>
            )}
            {inputFormat === 'selfies' && (
              <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
                Example: [C][N][C][=Branch1][C][=O][N][=Branch2][C][=Branch1][C][=O][N][Ring1][C]
              </p>
            )}
          </div>

          {/* Conversion Direction Indicator */}
          <div className="flex items-center justify-center py-4">
            <div className="flex-grow h-px bg-gray-200 dark:bg-gray-700"></div>
            <div className="mx-4 bg-gray-100 dark:bg-gray-700 p-2 rounded-full ring-1 ring-gray-300 dark:ring-gray-600">
              <HiOutlineArrowRight className="h-6 w-6 text-blue-600 dark:text-blue-400" />
            </div>
            <div className="flex-grow h-px bg-gray-200 dark:bg-gray-700"></div>
          </div>

          {/* Output Format Selection */}
          <div>
            <label htmlFor="output-format-select" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
              Output Format
            </label>
            <select
              id="output-format-select"
              value={outputFormat}
              onChange={(e) => setOutputFormat(e.target.value)}
              className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              disabled={inputFormat !== 'smiles'} // Disable selection if input is IUPAC or SELFIES
            >
              {OUTPUT_FORMAT_OPTIONS.map((option) => (
                <option key={option.id} value={option.id}>
                  {option.label}
                </option>
              ))}
            </select>
            {(inputFormat === 'iupac' || inputFormat === 'selfies') && (
              <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
                {inputFormat === 'iupac' ? 'IUPAC names' : 'SELFIES'} can only be converted to SMILES format
              </p>
            )}
          </div>

          {/* Toolkit Selection (conditionally shown) */}
          {showToolkitSelection && (
            <div>
              <label htmlFor="toolkit-select-convert" className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                Toolkit
              </label>
              <select
                id="toolkit-select-convert"
                value={toolkit}
                onChange={(e) => setToolkit(e.target.value)}
                className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-md text-gray-900 dark:text-white focus:ring-indigo-500 focus:border-indigo-500 shadow-sm"
              >
                {TOOLKIT_OPTIONS.map((option) => (
                  <option key={option.id} value={option.id}>
                    {option.label}
                  </option>
                ))}
              </select>
              <p className="mt-1 text-xs text-gray-500 dark:text-gray-400">
                Note: Toolkit support may vary for different format conversions.
              </p>
            </div>
          )}

          {/* Submit Button */}
          <div className="pt-2">
            <button
              type="submit"
              disabled={!input.trim() || loading}
              className={`w-full sm:w-auto px-6 py-2 rounded-lg text-white font-medium flex items-center justify-center transition-colors duration-200 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-blue-500 ${!input.trim() || loading
                  ? 'bg-gray-400 dark:bg-gray-600 cursor-not-allowed'
                  : 'bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 shadow-sm'
                }`}
            >
              <HiOutlineSwitchHorizontal className="mr-2 h-5 w-5" aria-hidden="true" />
              {loading ? 'Converting...' : 'Convert Format'}
            </button>
          </div>
        </form>
      </div>

      {/* Loading State */}
      {loading && <LoadingScreen text="Converting format..." />}

      {/* Error Display */}
      {error && !loading && (
        <div className="p-4 rounded-md bg-red-50 dark:bg-red-900 dark:bg-opacity-30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700 flex items-start shadow" role="alert">
          <HiOutlineExclamationCircle className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400" aria-hidden="true" />
          <span>{error}</span>
        </div>
      )}

      {/* Results Display Section */}
      {result && !loading && (
        <div className="bg-white dark:bg-gray-800 p-6 rounded-lg shadow-md dark:shadow-lg border border-gray-200 dark:border-gray-700">
          {/* Results Header */}
          <div className="flex justify-between items-center mb-3 border-b border-gray-200 dark:border-gray-700 pb-2">
            <h3 className="text-lg font-semibold text-gray-900 dark:text-white">Result</h3>
            {/* Copy Button */}
            <button
              onClick={handleCopyResult}
              className={`p-1.5 rounded-md transition-colors focus:outline-none focus:ring-1 focus:ring-blue-500 ${copied
                  ? 'text-green-500 dark:text-green-500'
                  : 'text-gray-500 dark:text-gray-400 hover:text-gray-800 dark:hover:text-white hover:bg-gray-100 dark:hover:bg-gray-700'
                }`}
              title={copied ? "Copied!" : "Copy result to clipboard"}
              aria-label={copied ? "Result Copied" : "Copy Result"}
            >
              {copied ? (
                <HiOutlineCheck className="h-5 w-5" />
              ) : (
                <HiOutlineClipboard className="h-5 w-5" />
              )}
            </button>
          </div>

          {/* Result Display Box */}
          <div className="mt-2 p-3 bg-gray-100 dark:bg-gray-900 rounded-md font-mono text-sm overflow-x-auto border border-gray-200 dark:border-gray-700 shadow-sm">
            <pre className="whitespace-pre-wrap break-all text-gray-700 dark:text-gray-300">{result}</pre>
          </div>

          {/* Conversion Info Text */}
          <div className="mt-4 text-xs text-gray-500 dark:text-gray-400">
            Converted from {INPUT_FORMAT_OPTIONS.find(o => o.id === inputFormat)?.label || inputFormat.toUpperCase()}
            to {OUTPUT_FORMAT_OPTIONS.find(o => o.id === outputFormat)?.label || outputFormat.toUpperCase()}
            {showToolkitSelection && ` using ${TOOLKIT_OPTIONS.find(o => o.id === toolkit)?.label || toolkit}`}
            {showIupacConverterSelection && ` with ${IUPAC_CONVERTER_OPTIONS.find(o => o.id === iupacConverter)?.label || iupacConverter}`}.
          </div>
        </div>
      )}

      {/* Initial State Message */}
      {!result && !loading && !error && (
        <div className="bg-blue-50 dark:bg-blue-900 dark:bg-opacity-20 border border-blue-200 dark:border-blue-800 rounded-lg p-6 text-center shadow">
          <p className="text-gray-600 dark:text-gray-300">
            Enter input data and select options to perform format conversion.
          </p>
        </div>
      )}
    </div>
  );
};

export default FormatConversionView;