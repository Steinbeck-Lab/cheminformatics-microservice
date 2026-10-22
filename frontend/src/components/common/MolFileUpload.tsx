// Description: A reusable React component for uploading and converting MOL/SDF files to SMILES
import { useState, useRef, useId, type ChangeEvent, type DragEvent } from "react";
import convertService from "../../services/convertService";
import { AlertCircle, CheckCircle, FileText, Upload, X } from "lucide-react";
import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";

interface MolFileUploadProps {
  onConversionSuccess?: (smiles: string, molblock: string, filename: string) => void;
  onConversionError?: (error: string) => void;
  toolkit?: string;
  className?: string;
  showMolblock?: boolean;
  allowMultipleMolecules?: boolean;
}

/**
 * MolFileUpload Component
 *
 * A reusable component that allows users to upload MOL/SDF files and convert them to SMILES using CDK or RDKit.
 *
 * @param {Object} props
 * @param {Function} props.onConversionSuccess - Callback when conversion is successful, receives (smiles, molblock, filename)
 * @param {Function} props.onConversionError - Callback when conversion fails, receives (error)
 * @param {string} props.toolkit - Toolkit to use for conversion: "cdk" (default) or "rdkit"
 * @param {string} props.className - Additional CSS classes for the container
 * @param {boolean} props.showMolblock - Whether to display the MOL block content (default: false)
 * @param {boolean} props.allowMultipleMolecules - Whether to allow SDF files with multiple molecules (default: false)
 */
const MolFileUpload = ({
  onConversionSuccess,
  onConversionError,
  toolkit = "cdk",
  className = "",
  showMolblock = false,
  allowMultipleMolecules = false,
}: MolFileUploadProps) => {
  const [molblock, setMolblock] = useState("");
  const [filename, setFilename] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [smiles, setSmiles] = useState("");
  const [dragging, setDragging] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const inputId = useId();

  // Helper function to format molblock for backend
  const formatMolblockForBackend = (molblock: string) => {
    // Normalize line endings first
    let formatted = molblock.replace(/\r\n/g, "\n").replace(/\r/g, "\n");

    // Clean up - remove trailing spaces from each line and trim
    formatted = formatted
      .split("\n")
      .map((line: string) => line.trimEnd())
      .join("\n")
      .trim();

    // Ensure proper molblock format - starts with newline
    if (!formatted.startsWith("\n")) {
      formatted = "\n" + formatted;
    }

    // Ensure it ends with M  END followed by a newline
    if (!formatted.endsWith("M  END\n")) {
      if (formatted.endsWith("M  END")) {
        formatted += "\n";
      }
    }

    return formatted;
  };

  // Read, validate and convert one file (from the picker or a drop)
  const processFile = (file: File | undefined) => {
    if (!file) return;

    // Reset state
    setError(null);
    setSmiles("");
    setMolblock("");
    setFilename("");

    // Validate file extension
    const fileName = file.name.toLowerCase();
    if (!fileName.endsWith(".mol") && !fileName.endsWith(".sdf")) {
      const errorMsg = "Please upload a valid .mol or .sdf file.";
      setError(errorMsg);
      if (onConversionError) onConversionError(errorMsg);
      return;
    }

    setFilename(file.name);
    setLoading(true);

    const reader = new FileReader();

    reader.onload = async (event) => {
      try {
        let fileContent = String(event.target?.result ?? "");

        // Normalize line endings and trim
        fileContent = fileContent.replace(/\r\n/g, "\n").replace(/\r/g, "\n").trim();

        // Validate molblock format
        if (!fileContent.includes("M  END")) {
          throw new Error(
            'Invalid file format: The uploaded file does not appear to be a valid MOL/SDF file (missing "M  END").'
          );
        }

        // Check for multiple molecules if not allowed
        if (!allowMultipleMolecules) {
          const molBlockCount = (fileContent.match(/M {2}END/g) || []).length;
          if (molBlockCount > 1) {
            throw new Error(
              "Multiple molecules detected in file. This component only supports single molecules. Please upload a file containing only one molecule."
            );
          }
        }

        // Handle SDF terminator - extract first molecule
        if (fileContent.includes("$$$$")) {
          fileContent = (fileContent.split("$$$$")[0] ?? "").trim();
        }

        setMolblock(fileContent);

        // Format and convert to SMILES
        const formattedMolblock = formatMolblockForBackend(fileContent);
        const convertedSmiles = await convertService.molblockToSMILES(formattedMolblock, toolkit);

        // Clean up SMILES (remove quotes if present)
        const cleanSmiles = convertedSmiles.replace(/^"|"$/g, "").trim();

        setSmiles(cleanSmiles);
        setError(null);

        // Call success callback
        if (onConversionSuccess) {
          onConversionSuccess(cleanSmiles, fileContent, file.name);
        }
      } catch (err) {
        console.error("File conversion error:", err);
        const errorMsg =
          (err instanceof Error && err.message) || "Failed to convert MOL/SDF file to SMILES.";
        setError(errorMsg);
        if (onConversionError) onConversionError(errorMsg);
      } finally {
        setLoading(false);
      }
    };

    reader.onerror = () => {
      const errorMsg = "Failed to read the uploaded file. Please ensure it's a valid text file.";
      setError(errorMsg);
      setLoading(false);
      if (onConversionError) onConversionError(errorMsg);
    };

    reader.readAsText(file);
  };

  const handleFileUpload = (e: ChangeEvent<HTMLInputElement>) => {
    processFile(e.target.files?.[0]);
  };

  const handleDragOver = (e: DragEvent<HTMLLabelElement>) => {
    e.preventDefault();
    if (!loading) setDragging(true);
  };

  const handleDrop = (e: DragEvent<HTMLLabelElement>) => {
    e.preventDefault();
    setDragging(false);
    if (!loading) processFile(e.dataTransfer.files?.[0]);
  };

  // Clear/reset the component
  const handleClear = () => {
    setMolblock("");
    setFilename("");
    setSmiles("");
    setError(null);
    if (fileInputRef.current) {
      fileInputRef.current.value = "";
    }
  };

  return (
    <div className={cn("grid gap-3", className)}>
      {/* Dropzone */}
      <label
        htmlFor={inputId}
        onDragOver={handleDragOver}
        onDragLeave={() => setDragging(false)}
        onDrop={handleDrop}
        className={cn(
          "flex cursor-pointer flex-col items-center justify-center gap-1.5 rounded-2xl border border-dashed border-[var(--line-2)] px-4 py-6 text-center transition-colors",
          "hover:border-primary/45 focus-within:outline-2 focus-within:outline-offset-2 focus-within:outline-ring",
          dragging && "border-primary bg-[var(--inset)]",
          loading && "cursor-not-allowed opacity-60"
        )}
      >
        <input
          ref={fileInputRef}
          id={inputId}
          type="file"
          accept=".mol,.sdf"
          onChange={handleFileUpload}
          className="sr-only"
          disabled={loading}
        />
        <Upload className="size-5 text-muted-foreground" aria-hidden="true" />
        <span className="text-sm font-medium text-foreground">
          {filename || "Drop a MOL/SDF file here or click to browse"}
        </span>
        <span className="pc-note">Accepts .mol and .sdf</span>
      </label>

      {loading && (
        <p className="pc-note" role="status">
          Converting to SMILES...
        </p>
      )}

      {error && !loading && (
        <div role="alert" className="pc-err on m-0 items-center">
          <AlertCircle className="pc-i size-4 shrink-0" aria-hidden="true" />
          <span className="flex-1">{error}</span>
          <Button
            variant="ghost"
            size="icon-xs"
            onClick={() => setError(null)}
            aria-label="Dismiss error"
          >
            <X aria-hidden="true" />
          </Button>
        </div>
      )}

      {smiles && !loading && !error && (
        <div className="grid gap-3">
          <div className="grid gap-2 rounded-2xl border border-[var(--line)] bg-[var(--inset)] p-3">
            <div className="flex items-center gap-2">
              <CheckCircle className="size-4 shrink-0 text-primary" aria-hidden="true" />
              <p className="m-0 flex-1 text-sm font-semibold text-foreground">
                Conversion successful
              </p>
              <Button
                variant="ghost"
                size="icon-sm"
                onClick={handleClear}
                title="Clear and upload new file"
                aria-label="Clear and upload new file"
              >
                <X aria-hidden="true" />
              </Button>
            </div>
            <p className="pc-well pc-mono m-0 min-h-0 py-2 text-sm break-all">{smiles}</p>
          </div>

          {showMolblock && molblock && (
            <div className="grid gap-2">
              <p className="m-0 flex items-center gap-2 text-xs font-semibold text-muted-foreground">
                <FileText className="size-4" aria-hidden="true" />
                MOL block content
              </p>
              <pre className="pc-well pc-mono m-0 max-h-60 overflow-auto text-xs">{molblock}</pre>
            </div>
          )}
        </div>
      )}

      {filename && !loading && (
        <p className="m-0 flex items-center gap-2 text-sm text-muted-foreground">
          <FileText className="size-4 shrink-0" aria-hidden="true" />
          <span className="truncate">{filename}</span>
        </p>
      )}
    </div>
  );
};
export default MolFileUpload;
