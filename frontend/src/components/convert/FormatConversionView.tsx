import { useEffect, useRef, useState, type ChangeEvent, type FormEvent } from "react";
import { ArrowLeftRight, Loader2, Upload } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MoleculeDepiction2D from "../depict/MoleculeDepiction2D";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { getErrorMessage } from "@/lib/error-messages";
import convertService from "../../services/convertService";
import XYZGridResult from "./XYZGridResult";
import type { XYZBatchConversionResult } from "@/types/api";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Seg, Select, Switch } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

type InputFormat = "smiles" | "iupac" | "selfies" | "molsdf" | "cdx" | "xyz";
type OutputFormat =
  | "smiles"
  | "canonicalsmiles"
  | "inchi"
  | "inchikey"
  | "cxsmiles"
  | "selfies"
  | "smarts"
  | "iupac"
  | "mol"
  | "sdf";
type Toolkit = "cdk" | "rdkit" | "openbabel";

const INPUT_FORMAT_OPTIONS: { id: InputFormat; label: string }[] = [
  { id: "smiles", label: "SMILES" },
  { id: "iupac", label: "IUPAC Name" },
  { id: "selfies", label: "SELFIES" },
  { id: "molsdf", label: "MOL/SDF Block" },
  { id: "cdx", label: "ChemDraw CDX/CDXML" },
  { id: "xyz", label: "XYZ Coordinates" },
];

const OUTPUT_FORMAT_OPTIONS: { id: OutputFormat; label: string }[] = [
  { id: "smiles", label: "SMILES" },
  { id: "canonicalsmiles", label: "Canonical SMILES" },
  { id: "inchi", label: "InChI" },
  { id: "inchikey", label: "InChI Key" },
  { id: "cxsmiles", label: "CXSMILES" },
  { id: "selfies", label: "SELFIES" },
  { id: "smarts", label: "SMARTS" },
  { id: "iupac", label: "IUPAC Name" },
  { id: "mol", label: "MOL Block" },
  { id: "sdf", label: "SDF" },
];

// Output formats reachable from XYZ input.
const XYZ_OUTPUT_FORMAT_IDS = new Set<OutputFormat>([
  "canonicalsmiles",
  "inchi",
  "inchikey",
  "mol",
  "sdf",
]);

const TOOLKIT_OPTIONS: { id: Toolkit; label: string; short: string }[] = [
  { id: "cdk", label: "CDK (Chemistry Development Kit)", short: "CDK" },
  { id: "rdkit", label: "RDKit", short: "RDKit" },
  { id: "openbabel", label: "OpenBabel", short: "OpenBabel" },
];

const IUPAC_CONVERTER_OPTIONS = [{ id: "opsin", label: "OPSIN" }];

// Outputs with no toolkit choice
const NO_TOOLKIT_OUTPUTS = new Set<OutputFormat>(["selfies", "smarts", "iupac"]);

// Orthonym naming tiers, strictest first
const IUPAC_TIER_OPTIONS = [
  { value: "pin", label: "PIN" },
  { value: "valid", label: "Valid" },
  { value: "best-effort", label: "Best effort" },
];

const labelOf = <T extends string>(list: { id: T; label: string }[], id: T) =>
  list.find((o) => o.id === id)?.label ?? id.toUpperCase();

// Detect whether input text is SMILES, SELFIES, or IUPAC name (local fallback)
const detectInputFormatLocal = (text: string): string | null => {
  const trimmed = text.trim();
  if (!trimmed) return null;

  // SELFIES: entirely composed of [token] sequences
  if (/^(\[[\w@=#+\-/\\%;.,^*]+\])+$/.test(trimmed)) return "selfies";

  // IUPAC: word-like name — only lowercase letters, digits, spaces, hyphens, commas,
  // parentheses; and no SMILES ring-closure patterns (letter immediately adjacent to digit)
  if (
    /^[a-z0-9][a-z0-9 ,\-().]*$/.test(trimmed) &&
    !/[a-z]\d/.test(trimmed) &&
    !/\d[a-z]/.test(trimmed)
  ) {
    return "iupac";
  }

  // SMILES: no whitespace, only valid SMILES characters
  if (!/\s/.test(trimmed) && /^[A-Za-z0-9[\]()=#@+\-/.\\%:*]+$/.test(trimmed)) return "smiles";

  // Otherwise assume IUPAC name
  return "iupac";
};

const mapBackendDetectedFormat = (format: string): string => {
  if (format === "molblock") return "molsdf";
  if (format === "xyz") return "xyz";
  if (format === "inchi") return "smiles";
  return format;
};

// Ensure a molblock is in the shape the backend expects
const formatMolblockForBackend = (molblock: string) => {
  // Normalize line endings, remove trailing spaces from each line and trim
  let formatted = molblock.replace(/\r\n/g, "\n").replace(/\r/g, "\n");
  formatted = formatted
    .split("\n")
    .map((line) => line.trimEnd())
    .join("\n")
    .trim();

  // V2000/V3000 molblocks have a 3-line header: name, program/timestamp, comment.
  // CDK-style molblocks omit the name line (start with program line), so we
  // prepend an empty name line.  Molblocks that already include a name line
  // (e.g. Actelion MolfileCreator) must NOT get an extra blank line.
  const lines = formatted.split("\n");
  const countsIdx = lines.findIndex((l) => /V[23]000/.test(l));
  if (countsIdx >= 0 && countsIdx <= 2) {
    // Counts line found too early → name line is missing → add blank name
    formatted = "\n" + formatted;
  }

  // Ensure it ends with M  END followed by a newline
  if (!formatted.endsWith("M  END\n") && formatted.endsWith("M  END")) {
    formatted += "\n";
  }
  return formatted;
};

/** Format conversion: between identifiers (SMILES, InChI, SELFIES, …) and structure files. */
const FormatConversionView = () => {
  const [input, setInput] = useState("");
  const [inputFormat, setInputFormat] = useState<InputFormat>("smiles");
  const [outputFormat, setOutputFormat] = useState<OutputFormat>("canonicalsmiles");
  const [toolkit, setToolkit] = useState<Toolkit>("cdk");
  const [iupacConverter, setIupacConverter] = useState("opsin");
  const [iupacTier, setIupacTier] = useState("pin");
  const [result, setResult] = useState("");
  // Title and "converted from … to …" line, fixed when the result arrives
  const [resultLabel, setResultLabel] = useState("");
  const [resultSummary, setResultSummary] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [autoDetected, setAutoDetected] = useState(false);
  const detectTimeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  useEffect(() => {
    return () => {
      if (detectTimeoutRef.current) {
        clearTimeout(detectTimeoutRef.current);
      }
    };
  }, []);
  // Molecular structure display
  const [smilesForStructure, setSmilesForStructure] = useState("");
  const [showStructure, setShowStructure] = useState(false);
  // MOL/SDF file upload
  const [uploadedFilename, setUploadedFilename] = useState("");
  const fileInputRef = useRef<HTMLInputElement>(null);
  // CDX/CDXML file upload
  const [cdxFile, setCdxFile] = useState<File | null>(null);
  const [cdxFilename, setCdxFilename] = useState("");
  const cdxFileInputRef = useRef<HTMLInputElement>(null);
  // XYZ file upload + bond-perception parameters
  const [xyzFilename, setXyzFilename] = useState("");
  const [xyzCharge, setXyzCharge] = useState(0);
  const [xyzUseHueckel, setXyzUseHueckel] = useState(false);
  const xyzFileInputRef = useRef<HTMLInputElement>(null);
  const [xyzBatchResult, setXyzBatchResult] = useState<XYZBatchConversionResult | null>(null);

  // Read an uploaded MOL/SDF file into the input
  const handleFileUpload = (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    setError(null);
    setUploadedFilename("");

    const fileName = file.name.toLowerCase();
    if (!fileName.endsWith(".mol") && !fileName.endsWith(".sdf")) {
      setError("Please upload a valid .mol or .sdf file.");
      return;
    }

    setUploadedFilename(file.name);

    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        const fileContent = String(event.target?.result ?? "")
          .replace(/\r\n/g, "\n")
          .replace(/\r/g, "\n")
          .trim();

        if (!fileContent.includes("M  END")) {
          setError(
            'Invalid file format: The uploaded file does not appear to be a valid MOL/SDF file (missing "M  END").'
          );
          setInput("");
          return;
        }

        const molBlockCount = (fileContent.match(/M {2}END/g) || []).length;
        if (molBlockCount > 1) {
          setError(
            "Multiple molecules detected in file. Please upload a file containing only one molecule."
          );
          setInput("");
          return;
        }

        setInput(fileContent);
        setError(null);
      } catch (err) {
        console.error("File reading error:", err);
        setError("Failed to process the uploaded file. Please ensure it's a valid MOL/SDF file.");
        setInput("");
      }
    };
    reader.onerror = () => {
      setError("Failed to read the uploaded file. Please ensure it's a valid text file.");
      setInput("");
    };
    reader.readAsText(file);
  };

  const handleClearFile = () => {
    setUploadedFilename("");
    setInput("");
    if (fileInputRef.current) {
      fileInputRef.current.value = "";
    }
  };

  const handleCdxFileUpload = (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    setError(null);
    const name = file.name.toLowerCase();
    if (!name.endsWith(".cdx") && !name.endsWith(".cdxml")) {
      setError("Please upload a valid .cdx or .cdxml file.");
      return;
    }
    setCdxFile(file);
    setCdxFilename(file.name);
  };

  const handleClearCdxFile = () => {
    setCdxFile(null);
    setCdxFilename("");
    if (cdxFileInputRef.current) {
      cdxFileInputRef.current.value = "";
    }
  };

  // Read an XYZ file as text into the input
  const handleXyzFileUpload = (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    setError(null);
    if (!file.name.toLowerCase().endsWith(".xyz")) {
      setError("Please upload a valid .xyz file.");
      return;
    }
    setXyzFilename(file.name);
    const reader = new FileReader();
    reader.onload = (event) => {
      const content = event.target?.result;
      setInput(typeof content === "string" ? content : "");
    };
    reader.onerror = () => setError("Failed to read the uploaded XYZ file.");
    reader.readAsText(file);
  };

  const handleClearXyzFile = () => {
    setXyzFilename("");
    setInput("");
    setXyzBatchResult(null);
    if (xyzFileInputRef.current) {
      xyzFileInputRef.current.value = "";
    }
  };

  // A new input format may lock the output format and narrow the toolkits
  const handleInputFormatChange = (format: InputFormat) => {
    setInputFormat(format);
    setAutoDetected(false);

    if (format !== "molsdf") handleClearFile();
    if (format !== "cdx") handleClearCdxFile();
    if (format !== "xyz") handleClearXyzFile();

    if (format === "iupac" || format === "selfies" || format === "molsdf") {
      setOutputFormat("smiles");
      // MOL/SDF supports CDK and RDKit only
      if (format === "molsdf" && toolkit === "openbabel") setToolkit("cdk");
    } else if (format === "cdx") {
      setOutputFormat("mol");
    } else if (format === "xyz") {
      // XYZ: default to canonical SMILES; cap toolkit to rdkit/openbabel.
      setOutputFormat("canonicalsmiles");
      if (toolkit === "cdk") setToolkit("rdkit");
    }
  };

  const applyDetectedFormat = (detected: string | null) => {
    if (!detected) {
      setAutoDetected(false);
      return;
    }
    const mapped = mapBackendDetectedFormat(detected) as InputFormat;
    if (mapped && mapped !== inputFormat) {
      setInputFormat(mapped);
      setAutoDetected(true);
      if (mapped === "iupac" || mapped === "selfies") {
        setOutputFormat("smiles");
      } else if (mapped === "smiles") {
        setOutputFormat("canonicalsmiles");
      } else if (mapped === "molsdf") {
        setOutputFormat("smiles");
        if (toolkit === "openbabel") setToolkit("cdk");
      } else if (mapped === "xyz") {
        setOutputFormat("canonicalsmiles");
        if (toolkit === "cdk") setToolkit("rdkit");
      }
    }
  };

  // Text input with auto-format detection (backend first, local fallback)
  const handleInputChange = (value: string) => {
    setInput(value);
    if (inputFormat === "molsdf") return;

    if (!value.trim()) {
      setAutoDetected(false);
      return;
    }

    if (detectTimeoutRef.current) {
      clearTimeout(detectTimeoutRef.current);
    }

    detectTimeoutRef.current = setTimeout(async () => {
      try {
        const { detected_format } = await convertService.detectInputFormat(value);
        applyDetectedFormat(detected_format);
      } catch {
        applyDetectedFormat(detectInputFormatLocal(value));
      }
    }, 300);
  };

  const handleOutputFormatChange = (format: OutputFormat) => {
    setOutputFormat(format);
    // SMARTS only supports RDKit
    if (format === "smarts") {
      setToolkit("rdkit");
    }
  };

  // Convert a SMILES string to the chosen output format
  const convertSmilesTo = async (smiles: string, options: { autoDetect: boolean }) => {
    switch (outputFormat) {
      case "smiles":
        return smiles;
      case "mol":
        // MOL Block output: generate 2D coordinates from SMILES
        return convertService.generate2DCoordinates(smiles, toolkit, options);
      case "canonicalsmiles":
        return convertService.generateCanonicalSMILES(smiles, toolkit, options);
      case "inchi":
        return convertService.generateInChI(smiles, toolkit, options);
      case "inchikey":
        return convertService.generateInChIKey(smiles, toolkit, options);
      case "cxsmiles":
        return convertService.generateCXSMILES(smiles, toolkit, options);
      case "selfies":
        return convertService.generateSELFIES(smiles, options);
      case "smarts":
        return convertService.generateSMARTS(smiles, toolkit, options);
      default:
        throw new Error(`Unsupported output format: ${outputFormat}`);
    }
  };

  // Toolkit choice applies to SMILES, MOL/SDF and XYZ input, not to SELFIES/SMARTS/IUPAC output
  const showToolkitSelection =
    (inputFormat === "smiles" || inputFormat === "molsdf" || inputFormat === "xyz") &&
    !NO_TOOLKIT_OUTPUTS.has(outputFormat);

  const showIupacConverterSelection = inputFormat === "iupac";

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const trimmedInput = input.trim();
    if (inputFormat !== "cdx" && !trimmedInput) {
      setError("Please enter input data.");
      setResult("");
      return;
    }

    if (inputFormat === "cdx" && !cdxFile) {
      setError("Please upload a .cdx or .cdxml file.");
      return;
    }

    setLoading(true);
    setError(null);
    setResult("");
    setXyzBatchResult(null);
    setSmilesForStructure("");
    setShowStructure(false);

    try {
      let convertedResult: string | undefined;
      let iupacTierNote = "";
      let smilesForDisplay = "";
      const convertOptions = { autoDetect: autoDetected };

      if (inputFormat === "xyz") {
        // XYZ → SMILES / Canonical SMILES / InChI / InChIKey / MOL / SDF
        const xyzToolkit = toolkit === "cdk" ? "rdkit" : toolkit;
        const batch = await convertService.convertXYZ(trimmedInput, {
          charge: xyzCharge,
          useHueckel: xyzUseHueckel,
          toolkit: xyzToolkit,
        });

        setXyzBatchResult(batch);

        // Multi-frame: the grid renders everything; suppress the single-result panel.
        if (batch.summary.total > 1) {
          setResult("");
          setSmilesForStructure("");
          setShowStructure(false);
          return;
        }

        // Single frame: pick the chosen output.
        const only = batch.structures[0];
        if (!only || !only.success) {
          throw new Error(only?.error || "XYZ conversion failed");
        }
        smilesForDisplay = only.canonicalsmiles;
        switch (outputFormat) {
          case "smiles":
          case "canonicalsmiles":
            convertedResult = only.canonicalsmiles;
            break;
          case "inchi":
            convertedResult = only.inchi;
            break;
          case "inchikey":
            convertedResult = only.inchikey;
            break;
          case "mol":
            convertedResult = only.molblock;
            break;
          case "sdf":
            convertedResult = batch.sdf;
            break;
          default:
            throw new Error(`Unsupported output format for XYZ input: ${outputFormat}.`);
        }
      } else if (inputFormat === "cdx" && cdxFile) {
        // CDX/CDXML → MOL Block
        const molblock = await convertService.convertCDXToMol(cdxFile);
        convertedResult = molblock;
        try {
          smilesForDisplay = await convertService.molblockToSMILES(molblock, "rdkit");
        } catch {
          // Structure preview is optional
        }
      } else {
        // SMILES directly; MOL/SDF, IUPAC and SELFIES via SMILES first
        let smiles: string;
        if (inputFormat === "molsdf") {
          smiles = await convertService.molblockToSMILES(
            formatMolblockForBackend(trimmedInput),
            toolkit
          );
        } else if (inputFormat === "smiles") {
          smiles = trimmedInput;
        } else {
          smiles = await convertService.generateSMILES(
            trimmedInput,
            inputFormat,
            inputFormat === "iupac" ? iupacConverter : undefined
          );
        }
        smilesForDisplay = smiles;
        if (outputFormat === "iupac") {
          const named = await convertService.generateIUPACName(smiles, iupacTier, convertOptions);
          convertedResult = named.name;
          iupacTierNote = ` with Orthonym (${named.tier.replace(/_/g, " ")})`;
        } else {
          convertedResult = await convertSmilesTo(smiles, convertOptions);
        }
      }

      if (!convertedResult) {
        setError(`Conversion resulted in empty output.`);
        setResult("");
        setSmilesForStructure("");
        setShowStructure(false);
      } else {
        let finalResult = String(convertedResult);

        // Remove surrounding double quotes if present
        if (finalResult.length >= 2 && finalResult.startsWith('"') && finalResult.endsWith('"')) {
          finalResult = finalResult.substring(1, finalResult.length - 1);
        }

        setResult(finalResult);
        setResultLabel(labelOf(OUTPUT_FORMAT_OPTIONS, outputFormat));
        setResultSummary(
          `Converted from ${labelOf(INPUT_FORMAT_OPTIONS, inputFormat)} to ${labelOf(
            OUTPUT_FORMAT_OPTIONS,
            outputFormat
          )}${showToolkitSelection ? ` using ${labelOf(TOOLKIT_OPTIONS, toolkit)}` : ""}${
            showIupacConverterSelection
              ? ` with ${labelOf(IUPAC_CONVERTER_OPTIONS, iupacConverter)}`
              : ""
          }${iupacTierNote}.`
        );

        if (smilesForDisplay && smilesForDisplay.trim()) {
          // Clean up SMILES string - remove quotes and trim
          let cleanedSmiles = smilesForDisplay.trim();
          if (cleanedSmiles.startsWith('"') && cleanedSmiles.endsWith('"')) {
            cleanedSmiles = cleanedSmiles.substring(1, cleanedSmiles.length - 1);
          }
          setSmilesForStructure(cleanedSmiles);
          setShowStructure(true);
        }
      }
    } catch (err) {
      console.error("Conversion failed:", err);
      setError(getErrorMessage("convert", err));
      setResult("");
    } finally {
      setLoading(false);
    }
  };

  const outputLocked = inputFormat !== "smiles" && inputFormat !== "xyz";
  const outputChoices = OUTPUT_FORMAT_OPTIONS.filter((option) => {
    if (inputFormat === "xyz") return XYZ_OUTPUT_FORMAT_IDS.has(option.id);
    // SDF output is only meaningful for XYZ input today.
    return option.id !== "sdf";
  });
  const toolkitChoices = TOOLKIT_OPTIONS.filter((option) => {
    // MOL/SDF only supports CDK and RDKit
    if (inputFormat === "molsdf") return option.id === "cdk" || option.id === "rdkit";
    // XYZ supports RDKit (xyz2mol) and OpenBabel; CDK lacks bond perception.
    if (inputFormat === "xyz") return option.id === "rdkit" || option.id === "openbabel";
    return true;
  });

  const xyzMulti = !!xyzBatchResult && xyzBatchResult.summary.total > 1 && !loading;
  const done = !!result && !loading;
  const title = done
    ? resultLabel
    : xyzMulti && xyzBatchResult
      ? `${xyzBatchResult.summary.total} structures`
      : "Results";
  const focusId =
    inputFormat === "molsdf"
      ? "molsdf-input"
      : inputFormat === "xyz"
        ? "xyz-input"
        : "smiles-input";

  useRevealOnResult("fc-results", loading, !!result, !!xyzBatchResult, !!error);

  return (
    <div className="pc-bento">
      <Cell
        id="fc-input"
        title="Structure input"
        className="pc-t-input"
        meta={autoDetected && <span className="pc-tag text-primary">Auto-detected</span>}
      >
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <Select
            label="Input format"
            value={inputFormat}
            onChange={(v) => handleInputFormatChange(v as InputFormat)}
            options={INPUT_FORMAT_OPTIONS.map((o) => [o.id, o.label])}
          />

          {showIupacConverterSelection && (
            <div className="grid gap-2">
              <Select
                label="IUPAC converter"
                value={iupacConverter}
                onChange={setIupacConverter}
                options={IUPAC_CONVERTER_OPTIONS.map((o) => [o.id, o.label])}
              />
              <p className="pc-note">OPSIN is used to convert IUPAC names to SMILES.</p>
            </div>
          )}

          {inputFormat === "xyz" ? (
            <div className="grid gap-2">
              <label htmlFor="xyz-input" className="pc-lbl">
                XYZ block
              </label>
              <textarea
                id="xyz-input"
                className="pc-well resize-y"
                value={input}
                onChange={(e) => setInput(e.target.value)}
                placeholder={
                  "Paste XYZ block here...\n\nExample (water):\n3\nwater\nO    0.0000    0.0000    0.0000\nH    0.7572    0.5860    0.0000\nH   -0.7572    0.5860    0.0000"
                }
                rows={10}
                spellCheck={false}
                required
                aria-required="true"
              />
              <p className="pc-note">
                Paste plain-text XYZ coordinates (atom count, comment, then element x y z) or upload
                a .xyz file.
              </p>
              {xyzFilename && (
                <div className="pc-row">
                  <span className="pc-tag">{xyzFilename}</span>
                  <span className="pc-note">Loaded</span>
                  <Button
                    type="button"
                    variant="destructive"
                    size="sm"
                    onClick={handleClearXyzFile}
                  >
                    Clear
                  </Button>
                </div>
              )}
              <input
                ref={xyzFileInputRef}
                id="xyz-file-upload"
                type="file"
                accept=".xyz,text/plain"
                onChange={handleXyzFileUpload}
                className="sr-only"
                tabIndex={-1}
                aria-hidden="true"
              />
            </div>
          ) : inputFormat === "cdx" ? (
            <div className="grid gap-2">
              <span className="pc-lbl">ChemDraw file</span>
              <div className="pc-row">
                <Button
                  type="button"
                  variant="outline"
                  onClick={() => cdxFileInputRef.current?.click()}
                >
                  <Upload />
                  {cdxFilename ? "Choose another file" : "Choose CDX or CDXML file"}
                </Button>
                {cdxFilename && (
                  <>
                    <span className="pc-tag">{cdxFilename}</span>
                    <span className="pc-note">Ready</span>
                    <Button
                      type="button"
                      variant="destructive"
                      size="sm"
                      onClick={handleClearCdxFile}
                    >
                      Clear
                    </Button>
                  </>
                )}
              </div>
              <p className="pc-note">
                A ChemDraw binary (.cdx) or XML (.cdxml) file; the output is a V2000 MOL block.
              </p>
              <input
                ref={cdxFileInputRef}
                id="cdx-file-upload"
                type="file"
                accept=".cdx,.cdxml"
                onChange={handleCdxFileUpload}
                className="sr-only"
                tabIndex={-1}
                aria-hidden="true"
              />
            </div>
          ) : inputFormat === "molsdf" ? (
            <div className="grid gap-2">
              <label htmlFor="molsdf-input" className="pc-lbl">
                MOL/SDF block
              </label>
              <textarea
                id="molsdf-input"
                className="pc-well resize-y"
                value={input}
                onChange={(e) => setInput(e.target.value)}
                placeholder="Paste MOL or SDF block here...&#10;&#10;Example:&#10;  CDK     09012308392D&#10;&#10;  2  1  0  0  0  0  0  0  0  0999 V2000&#10;    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0&#10;    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0&#10;  1  2  1  0  0  0  0&#10;M  END"
                rows={12}
                spellCheck={false}
                required
                aria-required="true"
              />
              <p className="pc-note">
                Paste or upload a MOL or SDF block; for SDF (contains $$$$) only the first molecule
                is processed.
              </p>
              {uploadedFilename && (
                <div className="pc-row">
                  <span className="pc-tag">{uploadedFilename}</span>
                  <span className="pc-note">Loaded</span>
                  <Button type="button" variant="destructive" size="sm" onClick={handleClearFile}>
                    Clear
                  </Button>
                </div>
              )}
              <input
                ref={fileInputRef}
                id="mol-file-upload"
                type="file"
                accept=".mol,.sdf"
                onChange={handleFileUpload}
                className="sr-only"
                tabIndex={-1}
                aria-hidden="true"
              />
            </div>
          ) : (
            <>
              <SMILESInput
                value={input}
                onChange={handleInputChange}
                label={
                  inputFormat === "smiles"
                    ? "SMILES"
                    : inputFormat === "iupac"
                      ? "IUPAC name"
                      : "SELFIES"
                }
                placeholder={
                  inputFormat === "smiles"
                    ? "Enter SMILES notation..."
                    : inputFormat === "iupac"
                      ? "Enter IUPAC chemical name..."
                      : "Enter SELFIES notation..."
                }
                required
              />
              {inputFormat === "iupac" && (
                <p className="pc-note">Example: 1,3,7-trimethylpurine-2,6-dione (caffeine)</p>
              )}
              {inputFormat === "selfies" && (
                <p className="pc-note break-all">
                  Example: [C][N][C][=Branch1][C][=O][N][=Branch2][C][=Branch1][C][=O][N][Ring1][C]
                </p>
              )}
            </>
          )}

          <div className="pc-row pc-actions">
            <Button
              type="submit"
              disabled={(inputFormat === "cdx" ? !cdxFile : !input.trim()) || loading}
            >
              {loading ? <Loader2 className="animate-spin" /> : <ArrowLeftRight />}
              {loading ? "Converting…" : "Convert format"}
            </Button>
            {inputFormat === "molsdf" && (
              <Button type="button" variant="outline" onClick={() => fileInputRef.current?.click()}>
                <Upload />
                Upload MOL/SDF file
              </Button>
            )}
            {inputFormat === "xyz" && (
              <Button
                type="button"
                variant="outline"
                onClick={() => xyzFileInputRef.current?.click()}
              >
                <Upload />
                Upload XYZ file
              </Button>
            )}
          </div>
        </form>
      </Cell>

      <Cell id="fc-opts" title="Options" className="pc-t-side">
        <div className="pc-tool-form">
          <div className="grid gap-2">
            <Select
              label="Output format"
              value={outputFormat}
              disabled={outputLocked}
              onChange={(v) => handleOutputFormatChange(v as OutputFormat)}
              options={outputChoices.map((o) => [o.id, o.label])}
            />
            {(inputFormat === "iupac" || inputFormat === "selfies" || inputFormat === "molsdf") && (
              <p className="pc-note">
                {inputFormat === "iupac"
                  ? "IUPAC names"
                  : inputFormat === "selfies"
                    ? "SELFIES"
                    : "MOL/SDF blocks"}{" "}
                can only be converted to SMILES.
              </p>
            )}
            {inputFormat === "cdx" && (
              <p className="pc-note">
                ChemDraw CDX/CDXML files are converted to a V2000 MOL block.
              </p>
            )}
            {outputFormat === "smarts" && (
              <p className="pc-note">
                SMARTS describes molecular patterns for substructure searching; it is only available
                with RDKit.
              </p>
            )}
          </div>

          {outputFormat === "iupac" && (
            <div className="pc-opt-group">
              <span className="pc-lbl" id="fc-iupac-tier">
                Naming tier
              </span>
              <Seg
                labelId="fc-iupac-tier"
                full
                value={iupacTier}
                onChange={setIupacTier}
                options={IUPAC_TIER_OPTIONS}
              />
              <p className="pc-note">
                Orthonym names are checked by OPSIN. PIN gives only preferred IUPAC names.
              </p>
            </div>
          )}

          {showToolkitSelection && (
            <div className="pc-opt-group">
              <span className="pc-lbl" id="fc-toolkit">
                Toolkit
              </span>
              <Seg
                labelId="fc-toolkit"
                full
                value={toolkit}
                onChange={setToolkit}
                options={toolkitChoices.map((o) => ({ value: o.id, label: o.short }))}
              />
              <p className="pc-note">
                {inputFormat === "molsdf"
                  ? "MOL/SDF conversion supports CDK and RDKit."
                  : "Toolkit support may vary for different format conversions."}
              </p>
            </div>
          )}

          {inputFormat === "xyz" && (
            <div className="pc-opt-group">
              <span className="pc-lbl">Bond perception</span>
              <label className="pc-field">
                <span>Net charge</span>
                <input
                  id="xyz-charge"
                  type="number"
                  className="pc-input disabled:cursor-not-allowed disabled:opacity-50"
                  min={-10}
                  max={10}
                  step={1}
                  value={xyzCharge}
                  onChange={(e) => {
                    const v = parseInt(e.target.value, 10);
                    setXyzCharge(Number.isFinite(v) ? v : 0);
                  }}
                  disabled={toolkit === "openbabel"}
                />
              </label>
              <p className="pc-note">
                e.g. -1 for acetate. Used by RDKit (xyz2mol); OpenBabel ignores charge.
              </p>
              <Switch
                checked={xyzUseHueckel}
                onChange={setXyzUseHueckel}
                disabled={toolkit === "openbabel"}
                label="Use extended Hückel (RDKit only)"
                hint="Slower; better for unusual valences"
              />
              <p className="pc-note">
                Bond orders are perceived from 3D coordinates. CDK is not supported (its core
                distribution lacks XYZ→bond perception).
              </p>
            </div>
          )}
        </div>
      </Cell>

      <Cell
        id="fc-results"
        title={title}
        className="pc-t-wide"
        meta={done && <CopyButton text={result} label="Copy result" />}
      >
        {loading && <ToolSkeleton variant="conversion" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document.getElementById(focusId)?.focus();
            }}
          />
        )}
        {!result && !loading && !error && !xyzMulti && (
          <p className="pc-empty">Enter a structure and choose the formats to convert it.</p>
        )}
        {xyzMulti && xyzBatchResult && <XYZGridResult result={xyzBatchResult} />}
        {done && (
          <div className={`grid gap-4 ${showStructure ? "lg:grid-cols-2" : ""}`}>
            <div className="grid content-start gap-2">
              <pre className="max-h-96 overflow-auto rounded-xl border border-border bg-muted p-3 font-mono text-sm break-all whitespace-pre-wrap text-foreground">
                {result}
              </pre>
              <p className="pc-note">{resultSummary}</p>
            </div>
            {showStructure && smilesForStructure && (
              <div className="grid content-start gap-2">
                <MoleculeDepiction2D
                  smiles={smilesForStructure}
                  title="Generated Structure"
                  toolkit="cdk"
                  showCIP={true}
                />
                <p className="pc-note">
                  Structure generated from SMILES:{" "}
                  <code className="pc-mono break-all">{smilesForStructure}</code>
                </p>
              </div>
            )}
          </div>
        )}
      </Cell>
    </div>
  );
};

export default FormatConversionView;
