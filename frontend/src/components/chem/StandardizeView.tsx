import { useRef, useState, type ChangeEvent, type FormEvent } from "react";
import { FileBarChart, Loader2, Upload } from "lucide-react";
import SMILESInput from "../common/SMILESInput";
import MoleculeCard from "../common/MoleculeCard";
import { generate2DCoordinates } from "../../services/convertService";
import { ToolSkeleton } from "@/components/feedback/ToolSkeleton";
import { GlassErrorCard } from "@/components/feedback/GlassErrorCard";
import { errorText, getErrorMessage } from "@/lib/error-messages";
import { Button } from "@/components/ui/button";
import { Cell, CopyButton, Seg } from "@/components/common/porcelain-ui";
import { useRevealOnResult } from "@/hooks/useRevealOnResult";

const API_BASE_URL = import.meta.env.VITE_API_URL || "https://dev.api.naturalproducts.net/latest";

type InputMethod = "molblock" | "smiles";

/** What /chem/standardize returns: the four known fields plus any extra properties. */
interface StandardizeData {
  standardized_mol: string;
  canonical_smiles: string;
  inchi: string;
  inchikey: string;
  [key: string]: unknown;
}

const KNOWN_KEYS = ["standardized_mol", "canonical_smiles", "inchi", "inchikey"];

const MOLBLOCK_PLACEHOLDER = `Paste your molblock here, or upload a .mol/.sdf file below...

Example format (simple molecule):

  CDK     08302311362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.9743    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    1.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3248    2.8125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6238    0.5625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  END`;

/** Normalize a molblock for the backend: line endings, trailing spaces, a missing name line. */
const formatMolblockForBackend = (molblock: string) => {
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
    formatted = "\n" + formatted;
  }

  // Ensure it ends with M  END followed by a newline
  if (!formatted.endsWith("M  END\n") && formatted.endsWith("M  END")) {
    formatted += "\n";
  }
  return formatted;
};

/** Standardize: one molecule (molblock or SMILES) through the ChEMBL curation pipeline. */
const StandardizeView = () => {
  const [inputMethod, setInputMethod] = useState<InputMethod>("molblock");
  const [molblock, setMolblock] = useState("");
  const [smiles, setSmiles] = useState("");
  const [standardizedData, setStandardizedData] = useState<StandardizeData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const fileRef = useRef<HTMLInputElement>(null);

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();

    let molblockToStandardize = "";

    if (inputMethod === "smiles") {
      const trimmedSmiles = smiles.trim();
      if (!trimmedSmiles) {
        setError("Please enter a SMILES string.");
        setStandardizedData(null);
        return;
      }

      if (
        trimmedSmiles.includes(".") &&
        trimmedSmiles.split(".").filter((s) => s.trim().length > 0).length > 1
      ) {
        setError(
          "Please enter only one molecule at a time. Multiple molecules (separated by dots) are not supported."
        );
        setStandardizedData(null);
        return;
      }

      try {
        setLoading(true);
        setError(null);
        // Convert SMILES to molblock first
        const rawMolblock = await generate2DCoordinates(trimmedSmiles, "rdkit");
        if (!rawMolblock || !rawMolblock.includes("M  END")) {
          throw new Error("Received invalid response from the server. Please try again later.");
        }
        molblockToStandardize = formatMolblockForBackend(rawMolblock);
      } catch (conversionError) {
        console.error("SMILES to molblock conversion error:", conversionError);
        setError(
          `Error converting SMILES to molblock: ${errorText(conversionError, "") || "Invalid SMILES string."}`
        );
        setStandardizedData(null);
        setLoading(false);
        return;
      }
    } else {
      const trimmedMolblock = molblock.trim();
      if (!trimmedMolblock) {
        setError("Please paste or upload a molblock.");
        setStandardizedData(null);
        return;
      }

      if (!trimmedMolblock.includes("M  END")) {
        setError(
          'Invalid molblock format: Missing "M  END" terminator. Please ensure you have pasted a complete molblock.'
        );
        setStandardizedData(null);
        return;
      }

      const molBlockCount = (trimmedMolblock.match(/M {2}END/g) || []).length;
      if (molBlockCount > 1) {
        setError(
          "Please provide only one molecule at a time. Multiple molecules in a single input are not supported."
        );
        setStandardizedData(null);
        return;
      }

      molblockToStandardize = formatMolblockForBackend(trimmedMolblock);
      setLoading(true);
      setError(null);
    }

    setStandardizedData(null);

    console.log("Sending molblock to standardize:", {
      inputMethod,
      molblockLength: molblockToStandardize.length,
      molblockPreview: molblockToStandardize.substring(0, 200) + "...",
      endsWithNewline: molblockToStandardize.endsWith("\n"),
      containsMEnd: molblockToStandardize.includes("M  END"),
    });

    try {
      const response = await fetch(`${API_BASE_URL}/chem/standardize`, {
        method: "POST",
        headers: {
          "Content-Type": "text/plain; charset=utf-8",
          Accept: "application/json",
        },
        body: molblockToStandardize,
      });

      if (!response.ok) {
        let errorMsg = `Error ${response.status}: ${response.statusText}`;
        try {
          const errorData = await response.json();
          errorMsg = errorData.detail || errorMsg;
        } catch {
          // Ignore if response is not JSON
        }
        throw new Error(errorMsg);
      }

      const data = (await response.json()) as StandardizeData;
      setStandardizedData(data);
    } catch (err) {
      console.error("Standardization error:", err);
      const rawMsg = errorText(err, "");

      if (rawMsg.includes("cannot access local variable")) {
        setError(
          "Error parsing molblock: The molblock format appears to be invalid or corrupted. Please ensure you have pasted a complete, valid molblock with proper structure and exactly one molecule."
        );
      } else if (rawMsg.includes("422") && inputMethod === "molblock") {
        setError(
          "Invalid molblock format: Please ensure your molblock is in proper MDL MOL format and contains exactly one valid molecule structure."
        );
      } else if (
        rawMsg.includes("Unable to parse molblock") ||
        rawMsg.includes("Failure parsing molblock")
      ) {
        setError(
          "Unable to parse molblock: Please check that your molblock is in valid MDL MOL format and try again."
        );
      } else if (rawMsg.includes("500") || rawMsg.includes("Internal Server Error")) {
        setError(
          "Server processing error: There was an issue processing your molblock. Please verify the format is correct and contains a single valid molecule."
        );
      } else {
        setError(getErrorMessage("chem", err));
      }
      setStandardizedData(null);
    } finally {
      setLoading(false);
    }
  };

  // Read an uploaded .mol/.sdf file into the molblock field
  const handleFileUpload = (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    e.target.value = "";
    if (!file) return;

    if (!file.name.toLowerCase().endsWith(".mol") && !file.name.toLowerCase().endsWith(".sdf")) {
      setError("Please upload a valid .mol or .sdf file.");
      return;
    }

    const reader = new FileReader();
    reader.onload = (event) => {
      const fileContent = String(event.target?.result ?? "")
        .replace(/\r\n/g, "\n")
        .replace(/\r/g, "\n")
        .trim();

      if (!fileContent.includes("M  END")) {
        setError(
          'Invalid file format: The uploaded file does not appear to be a valid molblock (missing "M  END").'
        );
        return;
      }

      const molBlockCount = (fileContent.match(/M {2}END/g) || []).length;
      if (molBlockCount > 1) {
        setError(
          "Multiple molecules detected in file. Please upload a file containing only one molecule."
        );
        return;
      }

      setMolblock(fileContent);
      setError(null);
    };
    reader.onerror = () => {
      console.error("File reading error:", reader.error);
      setError("Failed to read the uploaded file. Please ensure it's a valid text file.");
    };
    reader.readAsText(file);
  };

  const done = !!standardizedData && !loading;
  const extra = standardizedData
    ? Object.entries(standardizedData).filter(([key]) => !KNOWN_KEYS.includes(key))
    : [];
  const ids: [string, string][] = standardizedData
    ? [
        ["Canonical SMILES", standardizedData.canonical_smiles],
        ["InChI", standardizedData.inchi],
        ["InChI Key", standardizedData.inchikey],
      ]
    : [];

  useRevealOnResult("sd-results", loading, !!standardizedData, !!error);

  return (
    <div className="pc-bento">
      <Cell id="sd-input" title="Structure input" className="pc-t-input">
        <form onSubmit={handleSubmit} className="pc-tool-form">
          <div className="pc-opt-group">
            <span className="pc-lbl" id="sd-method">
              Input method
            </span>
            <Seg
              labelId="sd-method"
              full
              value={inputMethod}
              onChange={setInputMethod}
              options={[
                { value: "molblock", label: "Molblock" },
                { value: "smiles", label: "SMILES" },
              ]}
            />
          </div>

          {inputMethod === "smiles" ? (
            <>
              <SMILESInput
                value={smiles}
                onChange={setSmiles}
                placeholder="Enter SMILES string (e.g., CN1C=NC2=C1C(=O)N(C(=O)N2C)C)"
                required
              />
              <p className="pc-note">
                One SMILES string; it is converted to a molblock before standardization.
              </p>
            </>
          ) : (
            <div className="grid gap-2">
              <label htmlFor="molblock-input" className="pc-lbl">
                Molblock (MDL MOL format)
              </label>
              <textarea
                id="molblock-input"
                className="pc-well resize-y"
                value={molblock}
                onChange={(e) => setMolblock(e.target.value)}
                placeholder={MOLBLOCK_PLACEHOLDER}
                rows={12}
                spellCheck={false}
                required
                aria-required="true"
              />
              <p className="pc-note">
                One complete molecule in MDL MOL format, ending with &quot;M&nbsp;&nbsp;END&quot;.
              </p>
              <input
                ref={fileRef}
                type="file"
                accept=".mol,.sdf,chemical/x-mdl-molfile,chemical/x-mdl-sdfile"
                onChange={handleFileUpload}
                className="sr-only"
                tabIndex={-1}
                aria-hidden="true"
              />
            </div>
          )}

          <div className="pc-row pc-actions">
            <Button
              type="submit"
              disabled={(inputMethod === "smiles" ? !smiles.trim() : !molblock.trim()) || loading}
            >
              {loading ? <Loader2 className="animate-spin" /> : <FileBarChart />}
              {loading ? "Standardizing…" : "Standardize molecule"}
            </Button>
            {inputMethod === "molblock" && (
              <Button type="button" variant="outline" onClick={() => fileRef.current?.click()}>
                <Upload />
                Upload MOL/SDF file
              </Button>
            )}
          </div>
        </form>
      </Cell>

      <Cell id="sd-about" title="About standardization" className="pc-t-side">
        <p className="pc-tool-note">
          Standardizes a molecule with the ChEMBL curation pipeline for a consistent molecular
          representation.
        </p>
        <dl className="pc-defs">
          <div>
            <dt>One molecule at a time</dt>
            <dd>
              Several molecules (dot-separated SMILES or several records in a molblock) cause an
              error.
            </dd>
          </div>
        </dl>
      </Cell>

      <Cell
        id="sd-results"
        title={done ? "Standardized structure" : "Results"}
        className="pc-t-wide"
      >
        {loading && <ToolSkeleton variant="molecule" />}
        {error && !loading && (
          <GlassErrorCard
            message={error}
            onRetry={() => {
              setError(null);
              document
                .getElementById(inputMethod === "smiles" ? "smiles-input" : "molblock-input")
                ?.focus();
            }}
          />
        )}
        {!standardizedData && !loading && !error && (
          <p className="pc-empty">Enter a molblock or a SMILES string to standardize it.</p>
        )}
        {done && standardizedData && (
          <div className="pc-tool-form">
            <div className="grid gap-4 md:grid-cols-2">
              <MoleculeCard
                smiles={standardizedData.canonical_smiles}
                title="Standardized Structure"
                showActions
                size="md"
              />
              <dl className="pc-ids">
                {ids.map(([label, value]) => (
                  <div className="pc-id" key={label}>
                    <dt>{label}</dt>
                    <dd>
                      {value}
                      <CopyButton icon text={value} label={`Copy ${label}`} />
                    </dd>
                  </div>
                ))}
              </dl>
            </div>
            <div className="grid gap-2">
              <span className="pc-lbl">Standardized molblock</span>
              <pre className="max-h-60 overflow-auto rounded-xl border border-border bg-muted p-3 font-mono text-xs whitespace-pre text-foreground">
                {standardizedData.standardized_mol}
              </pre>
              <div className="pc-row">
                <CopyButton text={standardizedData.standardized_mol} label="Copy molblock" />
              </div>
            </div>
          </div>
        )}
      </Cell>

      {done && extra.length > 0 && (
        <Cell id="sd-extra" title="Additional properties" className="pc-t-wide">
          <dl className="pc-defs md:grid-cols-2">
            {extra.map(([key, value]) => (
              <div key={key}>
                <dt className="capitalize">{key.replace(/_/g, " ")}</dt>
                <dd className="font-mono break-all">
                  {typeof value === "object" ? JSON.stringify(value) : String(value)}
                </dd>
              </div>
            ))}
          </dl>
        </Cell>
      )}
    </div>
  );
};

export default StandardizeView;
